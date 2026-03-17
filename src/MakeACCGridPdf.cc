// src/MakeACCGridPdf.cc
#include "p2meg/MakeACCGridPdf.h"

#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <limits>

#include "TFile.h"
#include "TNamed.h"
#include "TParameter.h"
#include "TH1.h"
#include "TH2.h"
#include "THn.h"
#include "TAxis.h"
#include "TF1.h"
#include "TFitResultPtr.h"

#include "p2meg/AnalysisWindow.h"
#include "p2meg/AccTimeFit.h"
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
static constexpr int kNBins_tShape = 400;
static constexpr double kBlindCoreNs = 20.0; // [ns]

// ---- TSB の全時間範囲（基本案）----
// 実データの取得レンジに合わせて変更すること。
static constexpr double kTAllMin = -500.0; // [ns]
static constexpr double kTAllMax =  500.0; // [ns]

// ---- 最終解析用 Eg 条件付き time template のカテゴリ ----
// 低エネルギー側 [0,20) は診断用とみなし、ACC PDF の時間項には使わない。
static constexpr int kNAccTimeEgBins = 2;
static constexpr double kAccTimeEgEdges[kNAccTimeEgBins + 1] = {
    20.0, 30.0, 80.0
};

// ---- 因子化後の平滑化設定（E軸方向のみ）----
// 物理カットではなく、TSB 統計の疎さによる「穴」を埋めるための数値的平滑化。
static constexpr int    kSmoothRadiusBins = 2;   // 近傍の半径（ビン単位）
static constexpr double kSmoothSigmaBins  = 1.0; // ガウス重みの幅（ビン単位）
//============================================================
// 内部補助
//============================================================

static double GetBlindCoreNsFromEnv() {
  const char* s = std::getenv("P2MEG_ACC_BLIND_CORE_NS");
  if (!s) return kBlindCoreNs;
  char* endptr = nullptr;
  const double v = std::strtod(s, &endptr);
  if (endptr == s || !Math_IsFinite(v) || !(v > 0.0)) return kBlindCoreNs;
  return v;
}

static double GetSigmaMinNsFromEnv() {
  const char* s = std::getenv("P2MEG_ACC_SIGMA_MIN_NS");
  if (!s) return 20.0;
  char* endptr = nullptr;
  const double v = std::strtod(s, &endptr);
  if (endptr == s || !Math_IsFinite(v) || !(v > 0.0)) return 20.0;
  return v;
}

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

// 1Dヒストを count/bin から density[count/ns] へ変換する
static void ConvertCountsToDensity1D(TH1D& h) {
  const int nb = h.GetXaxis()->GetNbins();
  for (int ib = 1; ib <= nb; ++ib) {
    const double w = h.GetXaxis()->GetBinWidth(ib);
    if (!(w > 0.0) || !Math_IsFinite(w)) {
      h.SetBinContent(ib, 0.0);
      h.SetBinError(ib, 0.0);
      continue;
    }
    h.SetBinContent(ib, h.GetBinContent(ib) / w);
    h.SetBinError(ib, h.GetBinError(ib) / w);
  }
}

// piecewise-constant 密度を区間積分する
static double IntegrateDensityRange1D(const TH1D& h, double x_min, double x_max) {
  if (!Math_IsFinite(x_min) || !Math_IsFinite(x_max)) return 0.0;
  if (!(x_max > x_min)) return 0.0;

  const int nb = h.GetXaxis()->GetNbins();
  double sum = 0.0;
  for (int ib = 1; ib <= nb; ++ib) {
    const double lo = h.GetXaxis()->GetBinLowEdge(ib);
    const double hi = lo + h.GetXaxis()->GetBinWidth(ib);
    const double ov = std::min(hi, x_max) - std::max(lo, x_min);
    if (!(ov > 0.0)) continue;
    const double dens = h.GetBinContent(ib);
    if (dens > 0.0 && Math_IsFinite(dens)) sum += dens * ov;
  }
  return (sum > 0.0 && Math_IsFinite(sum)) ? sum : 0.0;
}

// TSB 積分（t_all 内かつ解析窓外）
static double IntegrateTimeSidebandDensity1D(const TH1D& h,
                                             double t_all_min,
                                             double t_all_max) {
  const double left = IntegrateDensityRange1D(h, t_all_min, analysis_window.t_min);
  const double right = IntegrateDensityRange1D(h, analysis_window.t_max, t_all_max);
  const double s = left + right;
  return (s > 0.0 && Math_IsFinite(s)) ? s : 0.0;
}

// TSB 由来の raw ヒストグラムを fit し、その関数から全時間域の density ヒストを作る
static int BuildFitBasedTimeTemplate(const TH1D& h_raw,
                                     TH1D& h_fit_out,
                                     int& fit_status_out,
                                     double& chi2_ndf_out,
                                     double pars_out[3]) {
  fit_status_out = -999;
  chi2_ndf_out = 0.0;
  pars_out[0] = 0.0;
  pars_out[1] = 0.0;
  pars_out[2] = 0.0;

  TH1D* h_tmp = dynamic_cast<TH1D*>(h_raw.Clone("acc_time_fit_input_tmp"));
  if (!h_tmp) return 1;

  AccTimeFitResult fitres{};
  const double blind_core_ns = GetBlindCoreNsFromEnv();
  const double sigma_min_ns = GetSigmaMinNsFromEnv();
  const int fit_ret =
      AccTimeFit_Run(*h_tmp, -blind_core_ns, blind_core_ns,
                     sigma_min_ns, 400.0, fitres);
  delete h_tmp;

  fit_status_out = fitres.fit_status;
  chi2_ndf_out = fitres.chi2_ndf;
  pars_out[0] = fitres.A;
  pars_out[1] = fitres.sigma;
  pars_out[2] = fitres.C;
  if (fit_ret != 0) return fit_ret;
  return AccTimeFit_FillDensityHistogramFromResult(h_fit_out, fitres);
}

// 4D 解析窓 (Ee, Eg, theta_eg, t) の内側かどうか
//  - 物理的な signal-like 領域を time template 学習から除外するために使う
static bool IsInsideFullAnalysisWindow(double Ee, double Eg,
                                       double theta_eg, double t) {
  if (!Math_IsFinite(Ee) || !Math_IsFinite(Eg) ||
      !Math_IsFinite(theta_eg) || !Math_IsFinite(t)) {
    return false;
  }
  if (Ee < analysis_window.Ee_min || Ee > analysis_window.Ee_max) return false;
  if (Eg < analysis_window.Eg_min || Eg > analysis_window.Eg_max) return false;
  if (theta_eg < analysis_window.theta_min || theta_eg > analysis_window.theta_max) return false;
  if (t < analysis_window.t_min || t > analysis_window.t_max) return false;
  return true;
}

static int FindAccTimeEgBin(double Eg) {
  if (!Math_IsFinite(Eg)) return -1;
  for (int i = 0; i < kNAccTimeEgBins; ++i) {
    const double lo = kAccTimeEgEdges[i];
    const double hi = kAccTimeEgEdges[i + 1];
    if ((Eg >= lo && Eg < hi) || (i == kNAccTimeEgBins - 1 && Eg >= lo && Eg <= hi)) {
      return i;
    }
  }
  return -1;
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
                                   long n_fill, long n_tshape_fill,
                                   long n_tshape_fit_success,
                                   double blind_core_ns,
                                   double total_mass,
                                   double norm_check, int N_phi_e, int N_phi_g,
                                   double t_all_min, double t_all_max,
                                   double w_tsb,
                                   double tshape_aw_mass,
                                   double tshape_tsb_mass,
                                   double s_tshape) {
  const double dphi_e = Detector_PhiStep(detres.phi_e_min, detres.phi_e_max, N_phi_e);
  const double dphi_g = Detector_PhiStep(detres.phi_g_min, detres.phi_g_max, N_phi_g);

  std::ostringstream oss;
  oss << "MakeACCGridPdf meta (4D grid, time sideband)\n";
  oss << "model: factorized p4(Ee,Eg,phi_e,phi_g)=pE(Ee,phi_e)*pG(Eg,phi_g)\n";
  oss << "smoothing: along E only (radius_bins=" << kSmoothRadiusBins
      << ", sigma_bins=" << kSmoothSigmaBins << ")\n";
  oss << "bins: Ee=" << kNBins_Ee << ", Eg=" << kNBins_Eg
      << ", phi_e=" << (N_phi_e + 1) << ", phi_g=" << (N_phi_g + 1) << "\n";
  oss << "phi axis: grid points i=0..N_phi (phi_i=phi_min + i*dphi)\n";
  oss << "phi integration: discrete grid (equal weight)\n";
  oss << "phi_e: range=[" << detres.phi_e_min << "," << detres.phi_e_max
      << "], dphi=" << dphi_e << " rad (upper edge phi_max+eps)\n";
  oss << "phi_g: range=[" << detres.phi_g_min << "," << detres.phi_g_max
      << "], dphi=" << dphi_g << " rad (upper edge phi_max+eps)\n";
  oss << "window: Ee=[" << analysis_window.Ee_min << "," << analysis_window.Ee_max << "] MeV\n";
  oss << "window: Eg=[" << analysis_window.Eg_min << "," << analysis_window.Eg_max << "] MeV\n";
  oss << "window: t=[" << analysis_window.t_min << "," << analysis_window.t_max << "] ns\n";
  oss << "t_all: [" << t_all_min << "," << t_all_max << "] ns\n";
  oss << "TSB: t in [t_all] and outside window (width=" << w_tsb << " ns)\n";
  oss << "tshape: Eg-conditional templates for final analysis bins\n";
  oss << "tshape source: fit function (1 Gaussian + constant) from TSB-only data\n";
  oss << "tshape blind core: |t|<" << blind_core_ns << " ns\n";
  oss << "tshape Eg bins: [20,30], [30,80] MeV\n";
  oss << "tshape_AW_mass=" << tshape_aw_mass << "\n";
  oss << "tshape_TSB_mass=" << tshape_tsb_mass << "\n";
  oss << "tshape_scale_AW_over_TSB=" << s_tshape << "\n";
  oss << "window: theta=[" << analysis_window.theta_min << "," << analysis_window.theta_max
      << "] rad (theta_eg=|phi_e-phi_g|)\n";
  oss << "normalization: sum_{phi_e,phi_g} integral dEe dEg p4 = 1\n";
  oss << "events_total=" << n_total << "\n";
  oss << "events_finite=" << n_finite << "\n";
  oss << "events_in_window=" << n_in_window << "\n";
  oss << "events_time_sideband=" << n_tsb << "\n";
  oss << "filled_entries=" << n_fill << "\n";
  oss << "tshape_filled_entries=" << n_tshape_fill << "\n";
  oss << "tshape_fit_success=" << n_tshape_fit_success << "\n";
  oss << "raw_mass=" << total_mass << "\n";
  oss << "norm_check_EeEg=" << norm_check << "\n";
  oss << "saved keys: <key>, <key>_tshape, <key>_tshape_egbin*\n";
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

  const int N_phi_e = detres.N_phi_e;
  const int N_phi_g = detres.N_phi_g;
  const double phi_e_min = detres.phi_e_min;
  const double phi_e_max = detres.phi_e_max;
  const double phi_g_min = detres.phi_g_min;
  const double phi_g_max = detres.phi_g_max;

  const std::vector<double> phi_edges_e =
      Detector_PhiEdgesFromGrid(phi_e_min, phi_e_max, N_phi_e);
  const std::vector<double> phi_edges_g =
      Detector_PhiEdgesFromGrid(phi_g_min, phi_g_max, N_phi_g);

  // ---- 4Dヒスト（Ee, Eg, phi_e, phi_g） ※最終出力 ----
  const int ndim = 4;
  int nbins[ndim] = {kNBins_Ee, kNBins_Eg, N_phi_e + 1, N_phi_g + 1};
  double xmin[ndim] = {analysis_window.Ee_min, analysis_window.Eg_min,
                       phi_e_min, phi_g_min};
  // ROOT のビンは [low, high) なので、phi=phi_max を in-range に入れるため上端を僅かに広げる
  double xmax[ndim] = {analysis_window.Ee_max, analysis_window.Eg_max,
                       phi_edges_e.back(), phi_edges_g.back()};

  THnD h("acc_grid_tmp", "ACC grid (phi);Ee;Eg;phi_e;phi_g", ndim, nbins, xmin, xmax);
  h.Sumw2();

  h.GetAxis(0)->SetTitle("Ee [MeV]");
  h.GetAxis(1)->SetTitle("Eg [MeV]");
  h.GetAxis(2)->SetTitle("phi_detector_e [rad]");
  h.GetAxis(3)->SetTitle("phi_detector_g [rad]");

  // phi 軸は「最近傍格子」になる可変ビンを設定
  h.GetAxis(2)->Set(N_phi_e + 1, phi_edges_e.data());
  h.GetAxis(3)->Set(N_phi_g + 1, phi_edges_g.data());

  // ---- 因子化用 2D ヒスト（Ee,phi_e）と（Eg,phi_g） ----
  // phi の軸定義は 4D と揃える（離散点が FindBin で一致するように）。
  TH2D hE("acc_e_tmp", "ACC factor pE;Ee [MeV];phi_detector_e [rad]",
          kNBins_Ee, analysis_window.Ee_min, analysis_window.Ee_max,
          N_phi_e + 1, phi_e_min, phi_edges_e.back());
  TH2D hG("acc_g_tmp", "ACC factor pG;Eg [MeV];phi_detector_g [rad]",
          kNBins_Eg, analysis_window.Eg_min, analysis_window.Eg_max,
          N_phi_g + 1, phi_g_min, phi_edges_g.back());
  hE.Sumw2();
  hG.Sumw2();
  hE.GetYaxis()->Set(N_phi_e + 1, phi_edges_e.data());
  hG.GetYaxis()->Set(N_phi_g + 1, phi_edges_g.data());

  // 後方互換用の全 Eg 時間テンプレート（Eg in analysis window）
  TH1D hTShapeRaw("acc_tshape_raw_tmp", "ACC raw time-shape (all Eg in analysis window);t [ns];density [arb/ns]",
               kNBins_tShape, kTAllMin, kTAllMax);
  hTShapeRaw.Sumw2();
  TH1D hTShapeEgRaw[kNAccTimeEgBins] = {
      TH1D("acc_tshape_egbin0_raw_tmp", "ACC raw time-shape Eg[20,30];t [ns];density [arb/ns]",
           kNBins_tShape, kTAllMin, kTAllMax),
      TH1D("acc_tshape_egbin1_raw_tmp", "ACC raw time-shape Eg[30,80];t [ns];density [arb/ns]",
           kNBins_tShape, kTAllMin, kTAllMax)
  };
  TH1D hTShape("acc_tshape_tmp", "ACC fit-based time-shape (all Eg in analysis window);t [ns];density [arb/ns]",
               kNBins_tShape, kTAllMin, kTAllMax);
  TH1D hTShapeEg[kNAccTimeEgBins] = {
      TH1D("acc_tshape_egbin0_tmp", "ACC fit-based time-shape Eg[20,30];t [ns];density [arb/ns]",
           kNBins_tShape, kTAllMin, kTAllMax),
      TH1D("acc_tshape_egbin1_tmp", "ACC fit-based time-shape Eg[30,80];t [ns];density [arb/ns]",
           kNBins_tShape, kTAllMin, kTAllMax)
  };
  for (int i = 0; i < kNAccTimeEgBins; ++i) {
    hTShapeEgRaw[i].Sumw2();
    hTShapeEg[i].Sumw2();
  }
  hTShape.Sumw2();

  long n_total = 0;
  long n_finite = 0;
  long n_in_window = 0;
  long n_tsb = 0;
  long n_fill = 0;
  long n_tshape_fill = 0;
  long n_tshape_aw_excluded = 0;
  long n_tshape_fit_success = 0;

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

    if (t >= kTAllMin && t <= kTAllMax &&
        Ee >= analysis_window.Ee_min && Ee <= analysis_window.Ee_max &&
        theta_eg >= analysis_window.theta_min && theta_eg <= analysis_window.theta_max) {
      const int ieg = FindAccTimeEgBin(Eg);
      if (ieg >= 0) {
        if (!IsInsideFullAnalysisWindow(Ee, Eg, theta_eg, t)) {
          hTShapeEgRaw[ieg].Fill(t);
          hTShapeRaw.Fill(t);
          ++n_tshape_fill;
        } else {
          ++n_tshape_aw_excluded;
        }
      }
    }

    if (!AnalysisWindow_In3D(analysis_window, Ee, Eg, theta_eg)) continue;
    ++n_in_window;

    if (!AnalysisWindow_InTSB(analysis_window, Ee, Eg, theta_eg, t, kTAllMin, kTAllMax)) {
      continue;
    }
    ++n_tsb;

    // 因子化：TSB から (Ee,phi_e) と (Eg,phi_g) を別々に詰める
    // phi は離散点として等重みで扱う
    hE.Fill(Ee, phi_e_disc);
    hG.Fill(Eg, phi_g_disc);
    ++n_fill;
  }

  std::cout << "[MakeACCGridPdf] events_total=" << n_total
            << " finite=" << n_finite
            << " in_window=" << n_in_window
            << " time_sideband=" << n_tsb
            << " filled=" << n_fill
            << " tshape_filled=" << n_tshape_fill
            << " tshape_aw_excluded=" << n_tshape_aw_excluded << "\n";

  ConvertCountsToDensity1D(hTShapeRaw);
  for (int i = 0; i < kNAccTimeEgBins; ++i) {
    ConvertCountsToDensity1D(hTShapeEgRaw[i]);
  }

  int fit_status_all = -999;
  double chi2_ndf_all = 0.0;
  double pars_all[3] = {0.0, 0.0, 0.0};
  if (BuildFitBasedTimeTemplate(hTShapeRaw, hTShape, fit_status_all, chi2_ndf_all, pars_all) == 0) {
    ++n_tshape_fit_success;
  } else {
    std::cerr << "[MakeACCGridPdf] fit-based global tshape build failed\n";
    return 2;
  }
  for (int i = 0; i < kNAccTimeEgBins; ++i) {
    int fit_status_eg = -999;
    double chi2_ndf_eg = 0.0;
    double pars_eg[3] = {0.0, 0.0, 0.0};
    if (BuildFitBasedTimeTemplate(hTShapeEgRaw[i], hTShapeEg[i], fit_status_eg, chi2_ndf_eg, pars_eg) == 0) {
      ++n_tshape_fit_success;
    } else {
      std::cerr << "[MakeACCGridPdf] fit-based Eg tshape build failed for bin " << i << "\n";
      return 2;
    }
  }

  const double tshape_aw_mass =
      IntegrateDensityRange1D(hTShape, analysis_window.t_min, analysis_window.t_max);
  const double tshape_tsb_mass =
      IntegrateTimeSidebandDensity1D(hTShape, kTAllMin, kTAllMax);
  const double s_tshape =
      (tshape_tsb_mass > 0.0) ? (tshape_aw_mass / tshape_tsb_mass) : 0.0;
  const double blind_core_ns = GetBlindCoreNsFromEnv();

  if (!(tshape_tsb_mass > 0.0) || !Math_IsFinite(tshape_tsb_mass)) {
    std::cerr << "[MakeACCGridPdf] t-shape TSB mass is not positive (AW="
              << tshape_aw_mass << ", TSB=" << tshape_tsb_mass << ")\n";
    return 2;
  }
  if (!(tshape_aw_mass > 0.0) || !Math_IsFinite(tshape_aw_mass)) {
    std::cerr << "[MakeACCGridPdf] t-shape AW mass is not positive even after fit (AW="
              << tshape_aw_mass << ", TSB=" << tshape_tsb_mass << ")\n";
    return 2;
  }

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

  TH1D hTShapeOut = hTShape;
  hTShapeOut.SetName((std::string(key) + "_tshape").c_str());
  hTShapeOut.Write();
  for (int i = 0; i < kNAccTimeEgBins; ++i) {
    TH1D hTShapeEgOut = hTShapeEg[i];
    hTShapeEgOut.SetName((std::string(key) + "_tshape_egbin" + std::to_string(i)).c_str());
    hTShapeEgOut.Write();
  }

  const std::string meta = BuildMetaString(n_total, n_finite, n_in_window,
                                           n_tsb, n_fill, n_tshape_fill, n_tshape_fit_success, blind_core_ns, total_mass,
                                           norm_check, N_phi_e, N_phi_g,
                                           kTAllMin, kTAllMax, w_tsb,
                                           tshape_aw_mass, tshape_tsb_mass, s_tshape);
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

  TParameter<double> par_tshape_aw((std::string(key) + "_tshape_aw_mass").c_str(), tshape_aw_mass);
  par_tshape_aw.Write();

  TParameter<double> par_tshape_tsb((std::string(key) + "_tshape_tsb_mass").c_str(), tshape_tsb_mass);
  par_tshape_tsb.Write();

  TParameter<double> par_tshape_scale((std::string(key) + "_tshape_scale").c_str(), s_tshape);
  par_tshape_scale.Write();

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

  std::cout << "[MakeACCGridPdf] saved (4D+tshape): " << out_filepath
            << " (key=" << key << ", " << key << "_tshape)\n";
  return 0;
}
