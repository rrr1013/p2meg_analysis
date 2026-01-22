// src/MakeACCGridPdf.cc
#include "p2meg/MakeACCGridPdf.h"

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <limits>

#include "TFile.h"
#include "TNamed.h"
#include "TParameter.h"
#include "TH2.h"
#include "THn.h"
#include "TAxis.h"

#include "p2meg/AnalysisWindow.h"
#include "p2meg/AnalysisWindowUtils.h"
#include "p2meg/AngleUtils.h"
#include "p2meg/DetectorResolution.h"
#include "p2meg/HistUtils.h"
#include "p2meg/MathUtils.h"

//============================================================
// 内部設定
//============================================================

// ---- 4D格子ビニング（Ee, Eg, phi_e, phi_g）----
static constexpr int kNBins_Ee = 40;
static constexpr int kNBins_Eg = 40;

// ---- TSB の全時間範囲（基本案）----
// 実データの取得レンジに合わせて変更すること。
static constexpr double kTAllMin = -10.0; // [ns]
static constexpr double kTAllMax =  10.0; // [ns]

// ---- 因子化後の平滑化設定（E軸方向のみ）----
// 物理カットではなく、TSB 統計の疎さによる「穴」を埋めるための数値的平滑化。
static constexpr int    kSmoothRadiusBins = 2;   // 近傍の半径（ビン単位）
static constexpr double kSmoothSigmaBins  = 1.0; // ガウス重みの幅（ビン単位）

//============================================================
// 内部補助
//============================================================

// 2D( E, phi ) の全ビン総和（カウントの総和）
static double SumAllBins2(const TH2D& h) {
  const int nx = h.GetXaxis()->GetNbins();
  const int ny = h.GetYaxis()->GetNbins();
  double sum = 0.0;
  for (int ix = 1; ix <= nx; ++ix) {
    for (int iy = 1; iy <= ny; ++iy) {
      const double v = h.GetBinContent(ix, iy);
      if (v > 0.0 && Math_IsFinite(v)) sum += v;
    }
  }
  return sum;
}

// 2D( E, phi ) を「密度」に変換し、総積分 Σ_phi ∫ dE p(E,phi) = 1 に正規化する。
//  - phi は離散変数として扱うため、phi のビン幅は正規化に含めない
//  - 密度の単位は [1/MeV]
static int ConvertToDensityAndNormalize2E(TH2D& h, double total_mass) {
  if (!(total_mass > 0.0) || !Math_IsFinite(total_mass)) return 1;

  const TAxis* axE = h.GetXaxis();
  const int nx = axE->GetNbins();
  const int ny = h.GetYaxis()->GetNbins();

  for (int ix = 1; ix <= nx; ++ix) {
    const double wE = axE->GetBinWidth(ix);
    if (!(wE > 0.0) || !Math_IsFinite(wE)) continue;

    for (int iy = 1; iy <= ny; ++iy) {
      const double C = h.GetBinContent(ix, iy);
      double density = 0.0;
      if (C > 0.0 && Math_IsFinite(C)) {
        density = (C / total_mass) / wE;
      }
      h.SetBinContent(ix, iy, density);
    }
  }
  return 0;
}

// 2D( E, phi ) 密度の正規化確認： Σ_phi ∫ dE p(E,phi) = 1
static double CheckNormalizationE2(const TH2D& h) {
  const TAxis* axE = h.GetXaxis();
  const int nx = axE->GetNbins();
  const int ny = h.GetYaxis()->GetNbins();

  double sum = 0.0;
  for (int ix = 1; ix <= nx; ++ix) {
    const double wE = axE->GetBinWidth(ix);
    if (!(wE > 0.0) || !Math_IsFinite(wE)) continue;

    for (int iy = 1; iy <= ny; ++iy) {
      const double p = h.GetBinContent(ix, iy);
      if (p > 0.0 && Math_IsFinite(p)) sum += p * wE;
    }
  }
  return sum;
}

// E軸方向のみ平滑化（各 phi ビンごとに 1D ガウス重み平均）
//  - 入力は「密度」(1/MeV)
//  - 端では重みを再正規化して連続性を保つ
static void SmoothAlongEPerPhi(TH2D& h, int radius_bins, double sigma_bins) {
  if (radius_bins <= 0) return;
  if (!(sigma_bins > 0.0) || !Math_IsFinite(sigma_bins)) return;

  const int nx = h.GetXaxis()->GetNbins();
  const int ny = h.GetYaxis()->GetNbins();

  std::vector<double> oldv(nx + 1, 0.0);
  std::vector<double> newv(nx + 1, 0.0);

  for (int iy = 1; iy <= ny; ++iy) {
    for (int ix = 1; ix <= nx; ++ix) {
      oldv[ix] = h.GetBinContent(ix, iy);
    }

    for (int ix = 1; ix <= nx; ++ix) {
      double sw = 0.0;
      double sv = 0.0;

      const int k0 = (ix - radius_bins >= 1) ? (ix - radius_bins) : 1;
      const int k1 = (ix + radius_bins <= nx) ? (ix + radius_bins) : nx;

      for (int k = k0; k <= k1; ++k) {
        const double dk = static_cast<double>(k - ix);
        const double w = std::exp(-0.5 * (dk * dk) / (sigma_bins * sigma_bins));
        const double v = oldv[k];
        if (w > 0.0 && Math_IsFinite(w) && v > 0.0 && Math_IsFinite(v)) {
          sw += w;
          sv += w * v;
        }
      }

      newv[ix] = (sw > 0.0 && Math_IsFinite(sw)) ? (sv / sw) : 0.0;
    }

    for (int ix = 1; ix <= nx; ++ix) {
      h.SetBinContent(ix, iy, newv[ix]);
    }
  }

  // 平滑化で数値積分がズレるのを抑えるため、全体で再正規化
  const double norm = CheckNormalizationE2(h);
  if (norm > 0.0 && Math_IsFinite(norm)) {
    for (int ix = 1; ix <= nx; ++ix) {
      for (int iy = 1; iy <= ny; ++iy) {
        const double v = h.GetBinContent(ix, iy);
        h.SetBinContent(ix, iy, (v > 0.0 && Math_IsFinite(v)) ? (v / norm) : 0.0);
      }
    }
  }
}

// 4Dヒストを「密度」に変換し、指定した total_mass で正規化する。
// density4 = (C / total_mass) / (dEe * dEg)
//  - phi は離散変数として扱うため、phi のビン幅は正規化に含めない
static int ConvertToDensityAndNormalize4EeEg(THnD& h, double total_mass) {
  const TAxis* ax0 = h.GetAxis(0);
  const TAxis* ax1 = h.GetAxis(1);
  const TAxis* ax2 = h.GetAxis(2);
  const TAxis* ax3 = h.GetAxis(3);

  const int n0 = ax0->GetNbins();
  const int n1 = ax1->GetNbins();
  const int n2 = ax2->GetNbins();
  const int n3 = ax3->GetNbins();

  if (!(total_mass > 0.0) || !Math_IsFinite(total_mass)) return 1;

  std::vector<int> idx(4, 1);

  for (int i0 = 1; i0 <= n0; ++i0) {
    idx[0] = i0;
    const double w0 = ax0->GetBinWidth(i0);
    for (int i1 = 1; i1 <= n1; ++i1) {
      idx[1] = i1;
      const double w1 = ax1->GetBinWidth(i1);
      for (int i2 = 1; i2 <= n2; ++i2) {
        idx[2] = i2;
        for (int i3 = 1; i3 <= n3; ++i3) {
          idx[3] = i3;

          const double vol = w0 * w1; // [MeV^2]
          if (!(vol > 0.0) || !Math_IsFinite(vol)) continue;

          const Long64_t bin = h.GetBin(idx.data());
          const double C = h.GetBinContent(bin);

          double density = 0.0;
          if (C > 0.0 && Math_IsFinite(C)) {
            density = (C / total_mass) / vol;
          }
          h.SetBinContent(bin, density);
        }
      }
    }
  }
  return 0;
}

// 密度の規格化確認（Ee/Eg のみ積分）
static double CheckNormalizationEeEg(const THnD& h) {
  const TAxis* ax0 = h.GetAxis(0);
  const TAxis* ax1 = h.GetAxis(1);

  const int n0 = ax0->GetNbins();
  const int n1 = ax1->GetNbins();
  const int n2 = h.GetAxis(2)->GetNbins();
  const int n3 = h.GetAxis(3)->GetNbins();

  std::vector<int> idx(4, 1);
  double sum = 0.0;

  for (int i0 = 1; i0 <= n0; ++i0) {
    idx[0] = i0;
    const double w0 = ax0->GetBinWidth(i0);
    for (int i1 = 1; i1 <= n1; ++i1) {
      idx[1] = i1;
      const double w1 = ax1->GetBinWidth(i1);
      const double vol = w0 * w1;
      if (!(vol > 0.0) || !Math_IsFinite(vol)) continue;
      for (int i2 = 1; i2 <= n2; ++i2) {
        idx[2] = i2;
        for (int i3 = 1; i3 <= n3; ++i3) {
          idx[3] = i3;
          const Long64_t bin = h.GetBin(idx.data());
          const double v = h.GetBinContent(bin);
          if (v > 0.0 && Math_IsFinite(v)) sum += v * vol;
        }
      }
    }
  }
  return sum;
}

static std::string BuildMetaString(long n_total, long n_finite,
                                   long n_in_window, long n_tsb,
                                   long n_fill, double total_mass,
                                   double norm_check, int N_phi_e, int N_phi_g,
                                   double t_all_min, double t_all_max,
                                   double w_tsb) {
  const double dphi_e = Detector_PhiStep(detres.phi_e_min, detres.phi_e_max, N_phi_e);
  const double dphi_g = Detector_PhiStep(detres.phi_g_min, detres.phi_g_max, N_phi_g);

  std::ostringstream oss;
  oss << "MakeACCGridPdf meta (4D grid, time sideband)\n";
  oss << "model: factorized p4(Ee,Eg,phi_e,phi_g)=pE(Ee,phi_e)*pG(Eg,phi_g)\n";
  oss << "smoothing: along E only (radius_bins=" << kSmoothRadiusBins
      << ", sigma_bins=" << kSmoothSigmaBins << ")\n";
  oss << "bins: Ee=" << kNBins_Ee << ", Eg=" << kNBins_Eg
      << ", phi_e=" << (N_phi_e + 1) << ", phi_g=" << (N_phi_g + 1) << "\n";
  oss << "phi_e range: [" << detres.phi_e_min << "," << detres.phi_e_max
      << "] rad, dphi_e=" << dphi_e << " rad\n";
  oss << "phi_g range: [" << detres.phi_g_min << "," << detres.phi_g_max
      << "] rad, dphi_g=" << dphi_g << " rad\n";
  oss << "window: Ee=[" << analysis_window.Ee_min << "," << analysis_window.Ee_max << "] MeV\n";
  oss << "window: Eg=[" << analysis_window.Eg_min << "," << analysis_window.Eg_max << "] MeV\n";
  oss << "window: t=[" << analysis_window.t_min << "," << analysis_window.t_max << "] ns\n";
  oss << "t_all: [" << t_all_min << "," << t_all_max << "] ns\n";
  oss << "TSB: t in [t_all] and outside window (width=" << w_tsb << " ns)\n";
  oss << "window: theta=[" << analysis_window.theta_min << "," << analysis_window.theta_max
      << "] rad (theta_eg=|phi_e-phi_g|)\n";
  oss << "normalization: sum_{phi_e,phi_g} integral dEe dEg p4 = 1\n";
  oss << "events_total=" << n_total << "\n";
  oss << "events_finite=" << n_finite << "\n";
  oss << "events_in_window=" << n_in_window << "\n";
  oss << "events_time_sideband=" << n_tsb << "\n";
  oss << "filled_entries=" << n_fill << "\n";
  oss << "raw_mass=" << total_mass << "\n";
  oss << "norm_check_EeEg=" << norm_check << "\n";
  oss << "saved keys: <key>\n";
  return oss.str();
}

//============================================================
// 本体
//============================================================

int MakeACCGridPdf(const std::vector<Event>& events,
                   const char* out_filepath,
                   const char* key) {
  if (!out_filepath || !key) {
    std::cerr << "[MakeACCGridPdf] invalid arguments\n";
    return 1;
  }

  const int N_phi_e = Math_GetNPhiE(detres);
  const int N_phi_g = Math_GetNPhiG(detres);
  if (!Detector_IsPhiRangeValid(detres.phi_e_min, detres.phi_e_max, N_phi_e) ||
      !Detector_IsPhiRangeValid(detres.phi_g_min, detres.phi_g_max, N_phi_g)) {
    std::cerr << "[MakeACCGridPdf] phi range is invalid.\n";
    return 2;
  }
  const double phi_e_max_plus =
      std::nextafter(detres.phi_e_max, std::numeric_limits<double>::infinity());
  const double phi_g_max_plus =
      std::nextafter(detres.phi_g_max, std::numeric_limits<double>::infinity());

  // ---- 4Dヒスト（Ee, Eg, phi_e, phi_g） ※最終出力 ----
  const int ndim = 4;
  int nbins[ndim] = {kNBins_Ee, kNBins_Eg, N_phi_e + 1, N_phi_g + 1};
  double xmin[ndim] = {analysis_window.Ee_min, analysis_window.Eg_min,
                       detres.phi_e_min, detres.phi_g_min};
  double xmax[ndim] = {analysis_window.Ee_max, analysis_window.Eg_max,
                       phi_e_max_plus, phi_g_max_plus};

  THnD h("acc_grid_tmp", "ACC grid (phi);Ee;Eg;phi_e;phi_g", ndim, nbins, xmin, xmax);
  h.Sumw2();

  h.GetAxis(0)->SetTitle("Ee [MeV]");
  h.GetAxis(1)->SetTitle("Eg [MeV]");
  h.GetAxis(2)->SetTitle("phi_detector_e [rad]");
  h.GetAxis(3)->SetTitle("phi_detector_g [rad]");

  const std::vector<double> phi_edges_e =
      Detector_PhiEdgesFromGrid(detres.phi_e_min, detres.phi_e_max, N_phi_e);
  const std::vector<double> phi_edges_g =
      Detector_PhiEdgesFromGrid(detres.phi_g_min, detres.phi_g_max, N_phi_g);
  h.GetAxis(2)->Set(N_phi_e + 1, phi_edges_e.data());
  h.GetAxis(3)->Set(N_phi_g + 1, phi_edges_g.data());

  // ---- 因子化用 2D ヒスト（Ee,phi_e）と（Eg,phi_g） ----
  // phi の軸定義は 4D と揃える（離散点が FindBin で一致するように）。
  TH2D hE("acc_e_tmp", "ACC factor pE;Ee [MeV];phi_detector_e [rad]",
          kNBins_Ee, analysis_window.Ee_min, analysis_window.Ee_max,
          N_phi_e + 1, phi_edges_e.data());
  TH2D hG("acc_g_tmp", "ACC factor pG;Eg [MeV];phi_detector_g [rad]",
          kNBins_Eg, analysis_window.Eg_min, analysis_window.Eg_max,
          N_phi_g + 1, phi_edges_g.data());
  hE.Sumw2();
  hG.Sumw2();

  long n_total = 0;
  long n_finite = 0;
  long n_in_window = 0;
  long n_tsb = 0;
  long n_fill = 0;

  const double w_tsb = AnalysisWindow_TimeSidebandWidth(analysis_window, kTAllMin, kTAllMax);

  for (const auto& ev : events) {
    ++n_total;
    const double Ee = ev.Ee;
    const double Eg = ev.Eg;
    const double t = ev.t;
    const double phi_e = ev.phi_detector_e;
    const double phi_g = ev.phi_detector_g;

    if (!Math_IsFinite(Ee) || !Math_IsFinite(Eg) || !Math_IsFinite(t) ||
        !Math_IsFinite(phi_e) || !Math_IsFinite(phi_g)) {
      continue;
    }
    ++n_finite;

    int idx_e = Detector_PhiIndexFromValue(phi_e, detres.phi_e_min, detres.phi_e_max, N_phi_e);
    int idx_g = Detector_PhiIndexFromValue(phi_g, detres.phi_g_min, detres.phi_g_max, N_phi_g);
    if (idx_e < 0 || idx_g < 0) continue;
    if (!Detector_IsAllowedPhiPairIndex(idx_e, idx_g, detres)) continue;

    const double phi_e_disc = Detector_PhiGridPoint(idx_e, detres.phi_e_min, detres.phi_e_max, N_phi_e);
    const double phi_g_disc = Detector_PhiGridPoint(idx_g, detres.phi_g_min, detres.phi_g_max, N_phi_g);
    const double theta_eg = std::fabs(phi_e_disc - phi_g_disc);

    if (!AnalysisWindow_In3D(analysis_window, Ee, Eg, theta_eg)) continue;
    ++n_in_window;

    if (!AnalysisWindow_InTSB(analysis_window, Ee, Eg, theta_eg, t, kTAllMin, kTAllMax)) {
      continue;
    }
    ++n_tsb;

    // 因子化：TSB から (Ee,phi_e) と (Eg,phi_g) を別々に詰める
    hE.Fill(Ee, phi_e_disc, 1.0);
    hG.Fill(Eg, phi_g_disc, 1.0);
    ++n_fill;
  }

  std::cout << "[MakeACCGridPdf] events_total=" << n_total
            << " finite=" << n_finite
            << " in_window=" << n_in_window
            << " time_sideband=" << n_tsb
            << " filled=" << n_fill << "\n";

  // ---- 2D 因子の密度化・正規化・平滑化 ----
  const double massE = SumAllBins2(hE);
  const double massG = SumAllBins2(hG);
  if (!(massE > 0.0) || !Math_IsFinite(massE) ||
      !(massG > 0.0) || !Math_IsFinite(massG)) {
    std::cerr << "[MakeACCGridPdf] factor hist mass is not positive (massE="
              << massE << ", massG=" << massG << ")\n";
    return 2;
  }

  if (ConvertToDensityAndNormalize2E(hE, massE) != 0 ||
      ConvertToDensityAndNormalize2E(hG, massG) != 0) {
    std::cerr << "[MakeACCGridPdf] normalization failed for factor hists\n";
    return 3;
  }

  // E方向のみ平滑化（phi ごと）
  SmoothAlongEPerPhi(hE, kSmoothRadiusBins, kSmoothSigmaBins);
  SmoothAlongEPerPhi(hG, kSmoothRadiusBins, kSmoothSigmaBins);

  const double normE = CheckNormalizationE2(hE);
  const double normG = CheckNormalizationE2(hG);
  std::cout << "[MakeACCGridPdf] norm_check_factor_E=" << normE
            << " norm_check_factor_G=" << normG << "\n";

  // ---- 4D 格子を因子化で構成（C=確率質量）----
  // 4Dのビン内容 C は、最終的に ConvertToDensityAndNormalize4EeEg で
  // density4 = (C/total_mass)/ (dEe*dEg) に変換される。
  // ここでは C = pE(Ee,phi_e) * pG(Eg,phi_g) * (dEe*dEg) としておく。
  const TAxis* axEe = h.GetAxis(0);
  const TAxis* axEg = h.GetAxis(1);
  const TAxis* axPe = h.GetAxis(2);
  const TAxis* axPg = h.GetAxis(3);

  const int nEe = axEe->GetNbins();
  const int nEg = axEg->GetNbins();
  const int nPe = axPe->GetNbins();
  const int nPg = axPg->GetNbins();

  std::vector<int> idx(4, 1);
  for (int iPe = 1; iPe <= nPe; ++iPe) {
    const int idx_e = iPe - 1;
    idx[2] = iPe;
    for (int iPg = 1; iPg <= nPg; ++iPg) {
      const int idx_g = iPg - 1;
      if (!Detector_IsAllowedPhiPairIndex(idx_e, idx_g, detres)) continue;
      idx[3] = iPg;
      for (int iEe = 1; iEe <= nEe; ++iEe) {
        idx[0] = iEe;
        const double wEe = axEe->GetBinWidth(iEe);
        const double pE = hE.GetBinContent(iEe, iPe); // [1/MeV]
        if (!(wEe > 0.0) || !Math_IsFinite(wEe) || !(pE > 0.0) || !Math_IsFinite(pE)) continue;

        for (int iEg = 1; iEg <= nEg; ++iEg) {
          idx[1] = iEg;
          const double wEg = axEg->GetBinWidth(iEg);
          const double pG = hG.GetBinContent(iEg, iPg); // [1/MeV]
          if (!(wEg > 0.0) || !Math_IsFinite(wEg) || !(pG > 0.0) || !Math_IsFinite(pG)) continue;

          const double C = pE * pG * (wEe * wEg); // 確率質量（無次元）
          if (!(C > 0.0) || !Math_IsFinite(C)) continue;

          const Long64_t bin = h.GetBin(idx.data());
          h.SetBinContent(bin, C);
        }
      }
    }
  }

  const double total_mass = Hist_SumAllBins4(h);
  if (!(total_mass > 0.0) || !Math_IsFinite(total_mass)) {
    std::cerr << "[MakeACCGridPdf] total mass is not positive (mass="
              << total_mass << ")\n";
    return 2;
  }

  if (ConvertToDensityAndNormalize4EeEg(h, total_mass) != 0) {
    std::cerr << "[MakeACCGridPdf] normalization failed for grid\n";
    return 3;
  }

  const double norm_check = CheckNormalizationEeEg(h);
  std::cout << "[MakeACCGridPdf] norm_check_EeEg=" << norm_check << "\n";

  // ---- 保存 ----
  TFile fout(out_filepath, "RECREATE");
  if (fout.IsZombie()) {
    std::cerr << "[MakeACCGridPdf] cannot open output file: " << out_filepath << "\n";
    return 4;
  }

  h.SetName(key);
  h.Write(key);

  const std::string meta = BuildMetaString(n_total, n_finite, n_in_window,
                                           n_tsb, n_fill, total_mass,
                                           norm_check, N_phi_e, N_phi_g,
                                           kTAllMin, kTAllMax, w_tsb);
  TNamed meta_obj((std::string(key) + "_meta").c_str(), meta.c_str());
  meta_obj.Write();

  TParameter<int> par_NphiE((std::string(key) + "_N_phi_e").c_str(), N_phi_e);
  par_NphiE.Write();

  TParameter<int> par_NphiG((std::string(key) + "_N_phi_g").c_str(), N_phi_g);
  par_NphiG.Write();

  TParameter<double> par_phi_e_min((std::string(key) + "_phi_e_min").c_str(), detres.phi_e_min);
  par_phi_e_min.Write();

  TParameter<double> par_phi_e_max((std::string(key) + "_phi_e_max").c_str(), detres.phi_e_max);
  par_phi_e_max.Write();

  TParameter<double> par_phi_g_min((std::string(key) + "_phi_g_min").c_str(), detres.phi_g_min);
  par_phi_g_min.Write();

  TParameter<double> par_phi_g_max((std::string(key) + "_phi_g_max").c_str(), detres.phi_g_max);
  par_phi_g_max.Write();

  TParameter<int> par_bins_Ee((std::string(key) + "_bins_Ee").c_str(), kNBins_Ee);
  par_bins_Ee.Write();

  TParameter<int> par_bins_Eg((std::string(key) + "_bins_Eg").c_str(), kNBins_Eg);
  par_bins_Eg.Write();

  TParameter<double> par_Ee_min((std::string(key) + "_Ee_min").c_str(), analysis_window.Ee_min);
  par_Ee_min.Write();

  TParameter<double> par_Ee_max((std::string(key) + "_Ee_max").c_str(), analysis_window.Ee_max);
  par_Ee_max.Write();

  TParameter<double> par_Eg_min((std::string(key) + "_Eg_min").c_str(), analysis_window.Eg_min);
  par_Eg_min.Write();

  TParameter<double> par_Eg_max((std::string(key) + "_Eg_max").c_str(), analysis_window.Eg_max);
  par_Eg_max.Write();

  TH2I hmask((std::string(key) + "_phi_mask").c_str(),
             "phi mask;phi_e index;phi_g index",
             N_phi_e + 1, -0.5, N_phi_e + 0.5,
             N_phi_g + 1, -0.5, N_phi_g + 0.5);
  for (int ie = 0; ie <= N_phi_e; ++ie) {
    for (int ig = 0; ig <= N_phi_g; ++ig) {
      const int allowed = Detector_IsAllowedPhiPairIndex(ie, ig, detres) ? 1 : 0;
      hmask.SetBinContent(ie + 1, ig + 1, allowed);
    }
  }
  hmask.Write();

  fout.Close();

  std::cout << "[MakeACCGridPdf] saved (4D): " << out_filepath
            << " (key=" << key << ")\n";
  return 0;
}
