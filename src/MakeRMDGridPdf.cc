// src/MakeRMDGridPdf.cc
#include "p2meg/MakeRMDGridPdf.h"

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TNamed.h"
#include "TParameter.h"
#include "TRandom3.h"
#include "THn.h"
#include "TAxis.h"

#include "p2meg/AnalysisWindow.h"
#include "p2meg/DetectorResolution.h"
#include "p2meg/RMDSpectrum.h"
#include "p2meg/Constants.h"

//============================================================
// 内部設定
//
// 単位：Ee, Eg は MeV、t は ns
// 角度は DetectorResolutionConst::N_theta による離散化（角度スメアなし）
// cosΔφ = +1 固定（Δφ=0）
//============================================================

// ---- エネルギー格子ビニング（Ee, Eg）----
static constexpr int kNBins_Ee = 40;
static constexpr int kNBins_Eg = 40;

// ---- 生成統計 ----
static constexpr long kNTruthSamples  = 1000000L; // 真値サンプル数（Ee,Eg を一様）
static constexpr int  kNSmearPerTruth = 50;       // 1真値あたりのスメア回数
static constexpr unsigned long kSeed  = 20251216UL;

// ---- RMD 理論関数の d_min ----
static constexpr double kDMin = 1e-6;

// ---- 平面内：cosΔφ=+1 固定 ----
static constexpr double kCosDeltaPhi = 1.0;

//============================================================
// 内部補助関数
//============================================================

static bool IsFinite(double x) {
  return std::isfinite(x);
}

static double Clamp(double x, double lo, double hi) {
  if (x < lo) return lo;
  if (x > hi) return hi;
  return x;
}

static bool IsInsideWindow2D(double Ee, double Eg) {
  if (Ee < analysis_window.Ee_min || Ee > analysis_window.Ee_max) return false;
  if (Eg < analysis_window.Eg_min || Eg > analysis_window.Eg_max) return false;
  return true;
}

// 離散角度：theta_i = i*pi/N_theta (i=0..N_theta)
static double ThetaFromIndex(int i, int N_theta) {
  if (!(N_theta >= 1)) return 0.0;
  if (i < 0) i = 0;
  if (i > N_theta) i = N_theta;
  return pi * static_cast<double>(i) / static_cast<double>(N_theta);
}

// (ie, ig) から cosThetaE, cosThetaG, cosThetaEG を計算（cosΔφ=+1 固定）
// 平面内で Δφ=0 とすると
//   cosThetaEG = cos(theta_e - theta_g)
//              = cosE*cosG + sinE*sinG
static bool AnglesFromIndices(int ie, int ig, int N_theta,
                              double& cosThetaE, double& cosThetaG,
                              double& cosThetaEG, double& thetaEG) {
  if (!(N_theta >= 1)) return false;
  if (ie < 0 || ie > N_theta) return false;
  if (ig < 0 || ig > N_theta) return false;

  const double thetaE = ThetaFromIndex(ie, N_theta);
  const double thetaG = ThetaFromIndex(ig, N_theta);

  const double cE = std::cos(thetaE);
  const double sE = std::sin(thetaE);
  const double cG = std::cos(thetaG);
  const double sG = std::sin(thetaG);

  // cosΔφ=+1（Δφ=0）
  const double cEG = cE * cG + sE * sG * kCosDeltaPhi;

  cosThetaE  = Clamp(cE,  -1.0, 1.0);
  cosThetaG  = Clamp(cG,  -1.0, 1.0);
  cosThetaEG = Clamp(cEG, -1.0, 1.0);
  thetaEG    = std::acos(cosThetaEG);
  return IsFinite(thetaEG);
}

// THnD（4D）の全ビン総和
static double SumAllBins4(const THnD& h) {
  const int n0 = h.GetAxis(0)->GetNbins();
  const int n1 = h.GetAxis(1)->GetNbins();
  const int n2 = h.GetAxis(2)->GetNbins();
  const int n3 = h.GetAxis(3)->GetNbins();

  std::vector<int> idx(4, 1);
  double sum = 0.0;

  for (int i0 = 1; i0 <= n0; ++i0) {
    idx[0] = i0;
    for (int i1 = 1; i1 <= n1; ++i1) {
      idx[1] = i1;
      for (int i2 = 1; i2 <= n2; ++i2) {
        idx[2] = i2;
        for (int i3 = 1; i3 <= n3; ++i3) {
          idx[3] = i3;
          const Long64_t bin = h.GetBin(idx.data());
          sum += h.GetBinContent(bin);
        }
      }
    }
  }
  return sum;
}

// THnD（4D）を「Ee,Eg の密度」に変換して、
//  Σ_{ie,ig} ∫∫ pdf(Ee,Eg,ie,ig) dEe dEg = 1 に正規化する
// 離散軸 (ie, ig) は幅1のビンとして扱う（測度は和）
// density4 = (C / total_mass) / (ΔEe * ΔEg)
static int ConvertToDensityAndNormalize4(THnD& h) {
  const TAxis* ax0 = h.GetAxis(0);
  const TAxis* ax1 = h.GetAxis(1);
  const TAxis* ax2 = h.GetAxis(2);
  const TAxis* ax3 = h.GetAxis(3);

  const int n0 = ax0->GetNbins();
  const int n1 = ax1->GetNbins();
  const int n2 = ax2->GetNbins();
  const int n3 = ax3->GetNbins();

  const double total_mass = SumAllBins4(h);
  if (!(total_mass > 0.0) || !IsFinite(total_mass)) return 1;

  std::vector<int> idx(4, 1);

  for (int i0 = 1; i0 <= n0; ++i0) {
    idx[0] = i0;
    const double w0 = ax0->GetBinWidth(i0);
    for (int i1 = 1; i1 <= n1; ++i1) {
      idx[1] = i1;
      const double w1 = ax1->GetBinWidth(i1);

      const double vol = w0 * w1; // [MeV^2]
      if (!(vol > 0.0) || !IsFinite(vol)) continue;

      for (int i2 = 1; i2 <= n2; ++i2) {
        idx[2] = i2;
        for (int i3 = 1; i3 <= n3; ++i3) {
          idx[3] = i3;

          const Long64_t bin = h.GetBin(idx.data());
          const double C = h.GetBinContent(bin);

          double density = 0.0;
          if (C > 0.0 && IsFinite(C)) {
            density = (C / total_mass) / vol;
          }
          h.SetBinContent(bin, density);
        }
      }
    }
  }
  return 0;
}

static std::string BuildMetaString(long n_setting_ok, long n_setting_ng,
                                   long n_wpos_ok, long n_filled4,
                                   int N_theta, double P_mu) {
  std::ostringstream oss;
  oss << "MakeRMDGridPdf meta (4D: Ee,Eg,ie,ig)\n";
  oss << "bins: Ee=" << kNBins_Ee << ", Eg=" << kNBins_Eg
      << ", ie=" << (N_theta + 1) << ", ig=" << (N_theta + 1) << "\n";
  oss << "truth_samples=" << kNTruthSamples << "\n";
  oss << "smear_per_truth=" << kNSmearPerTruth << "\n";
  oss << "seed=" << kSeed << "\n";
  oss << "d_min=" << kDMin << "\n";
  oss << "cosDeltaPhi=" << kCosDeltaPhi << "\n";
  oss << "N_theta=" << N_theta << "\n";
  oss << "P_mu=" << P_mu << "\n";

  oss << "setting_theta_window_ok=" << n_setting_ok << "\n";
  oss << "setting_theta_window_ng=" << n_setting_ng << "\n";
  oss << "wpos_ok=" << n_wpos_ok << "\n";
  oss << "filled_entries_4d=" << n_filled4 << "\n";

  oss << "window: Ee=[" << analysis_window.Ee_min << "," << analysis_window.Ee_max << "] MeV\n";
  oss << "window: Eg=[" << analysis_window.Eg_min << "," << analysis_window.Eg_max << "] MeV\n";
  oss << "window: t=[" << analysis_window.t_min << "," << analysis_window.t_max
      << "] ns (applied at evaluation)\n";
  oss << "window: theta=[" << analysis_window.theta_min << "," << analysis_window.theta_max
      << "] rad (applied via (ie,ig) -> thetaEG)\n";

  oss << "res: sigma_Ee=" << detres.sigma_Ee << " MeV\n";
  oss << "res: sigma_Eg=" << detres.sigma_Eg << " MeV\n";
  oss << "res: sigma_t=" << detres.sigma_t << " ns (used analytically at evaluation)\n";
  oss << "res: t_mean=" << detres.t_mean << " ns (used analytically at evaluation)\n";

  return oss.str();
}

//============================================================
// 本体
//============================================================

int MakeRMDGridPdf(const char* out_filepath, const char* key) {
  if (!out_filepath || !key) {
    std::cerr << "[MakeRMDGridPdf] invalid arguments\n";
    return 1;
  }

  // soft photon 発散対策：Eg > Eg_min
  if (!(analysis_window.Eg_min > 0.0)) {
    std::cerr << "[MakeRMDGridPdf] Eg_min must be > 0 to avoid soft photon divergence.\n";
    return 2;
  }

  // 分割数
  const int N_theta = detres.N_theta;
  if (!(N_theta >= 1)) {
    std::cerr << "[MakeRMDGridPdf] detres.N_theta must be >= 1\n";
    return 3;
  }

  // ---- 4Dヒスト（Ee, Eg, ie, ig）----
  const int ndim4 = 4;
  int nbins4[ndim4] = {kNBins_Ee, kNBins_Eg, N_theta + 1, N_theta + 1};

  // ie, ig は整数 0..N_theta をビン中央に入れるため [-0.5, N_theta+0.5]
  double xmin4[ndim4] = {
    analysis_window.Ee_min,
    analysis_window.Eg_min,
    -0.5,
    -0.5
  };
  double xmax4[ndim4] = {
    analysis_window.Ee_max,
    analysis_window.Eg_max,
    static_cast<double>(N_theta) + 0.5,
    static_cast<double>(N_theta) + 0.5
  };

  THnD h4("rmd_grid_tmp",
          "RMD smeared grid (4D);Ee;Eg;ie;ig",
          ndim4, nbins4, xmin4, xmax4);
  h4.Sumw2();
  h4.GetAxis(0)->SetTitle("Ee [MeV]");
  h4.GetAxis(1)->SetTitle("Eg [MeV]");
  h4.GetAxis(2)->SetTitle("i_e (theta_e=i*pi/N_theta)");
  h4.GetAxis(3)->SetTitle("i_g (theta_g=i*pi/N_theta)");

  TRandom3 rng(kSeed);

  const double Ee_min = analysis_window.Ee_min;
  const double Ee_max = analysis_window.Ee_max;
  const double Eg_min = analysis_window.Eg_min;
  const double Eg_max = analysis_window.Eg_max;

  const double sigma_Ee = detres.sigma_Ee;
  const double sigma_Eg = detres.sigma_Eg;

  if (!(sigma_Ee > 0.0) || !(sigma_Eg > 0.0)) {
    std::cerr << "[MakeRMDGridPdf] sigma_Ee and sigma_Eg must be > 0\n";
    return 4;
  }

  const double P_mu = detres.P_mu;

  long n_setting_ok = 0; // (ie,ig) で thetaEG が解析窓に入る回数
  long n_setting_ng = 0; // thetaEG が解析窓外だった回数
  long n_wpos_ok    = 0; // w>0 だった回数
  long n_filled4    = 0; // 4D に fill した回数

  for (long it = 0; it < kNTruthSamples; ++it) {
    const double Ee_true = rng.Uniform(Ee_min, Ee_max);
    const double Eg_true = rng.Uniform(Eg_min, Eg_max);

    // 設定（全部あり）：ie, ig を一様に選ぶ
    const int ie = static_cast<int>(rng.Integer(static_cast<ULong64_t>(N_theta + 1)));
    const int ig = static_cast<int>(rng.Integer(static_cast<ULong64_t>(N_theta + 1)));

    // 角度（離散）→ cos を計算し、thetaEG が解析窓に入る設定のみ使う
    double cosThetaE = 0.0, cosThetaG = 0.0, cosThetaEG = 0.0, thetaEG = 0.0;
    if (!AnglesFromIndices(ie, ig, N_theta, cosThetaE, cosThetaG, cosThetaEG, thetaEG)) {
      ++n_setting_ng;
      continue;
    }

    if (thetaEG < analysis_window.theta_min || thetaEG > analysis_window.theta_max) {
      ++n_setting_ng;
      continue;
    }
    ++n_setting_ok;

    // 理論重み（運動学的に許されない領域では 0）
    const double w0 = RMD_d6B_dEe_dEg_dOmegae_dOmegag(
      Ee_true, Eg_true,
      cosThetaEG, cosThetaE, cosThetaG,
      P_mu, kDMin
    );
    if (!(w0 > 0.0) || !IsFinite(w0)) continue;
    ++n_wpos_ok;

    // 重要度補正は不要（Ee,Eg,設定は一様サンプルで提案密度は定数）
    const double w = w0;

    // 同じ真値点から複数回スメアして観測分布を埋める（Ee,Eg のみ）
    for (int is = 0; is < kNSmearPerTruth; ++is) {
      const double Ee_obs = Ee_true + rng.Gaus(0.0, sigma_Ee);
      const double Eg_obs = Eg_true + rng.Gaus(0.0, sigma_Eg);

      if (!IsInsideWindow2D(Ee_obs, Eg_obs)) continue;

      double x4[4] = {Ee_obs, Eg_obs, static_cast<double>(ie), static_cast<double>(ig)};
      h4.Fill(x4, w);
      ++n_filled4;
    }
  }

  std::cout << "[MakeRMDGridPdf] setting theta-window ok: " << n_setting_ok
            << " / " << kNTruthSamples << "\n";
  std::cout << "[MakeRMDGridPdf] setting theta-window ng: " << n_setting_ng
            << " / " << kNTruthSamples << "\n";
  std::cout << "[MakeRMDGridPdf] w>0 ok: " << n_wpos_ok
            << " / " << n_setting_ok << "\n";
  std::cout << "[MakeRMDGridPdf] filled entries (4D): " << n_filled4 << "\n";

  // 正規化（Σ_{ie,ig} ∫ pdf dEe dEg = 1）
  const int norm4 = ConvertToDensityAndNormalize4(h4);
  if (norm4 != 0) {
    std::cerr << "[MakeRMDGridPdf] normalization (4D) failed\n";
    return 5;
  }

  // ---- 保存 ----
  TFile fout(out_filepath, "RECREATE");
  if (fout.IsZombie()) {
    std::cerr << "[MakeRMDGridPdf] cannot open output file: " << out_filepath << "\n";
    return 6;
  }

  h4.SetName(key);
  h4.Write(key);

  const std::string meta = BuildMetaString(n_setting_ok, n_setting_ng, n_wpos_ok, n_filled4, N_theta, P_mu);
  TNamed meta_obj((std::string(key) + "_meta").c_str(), meta.c_str());
  meta_obj.Write();

  TParameter<double> par_dmin((std::string(key) + "_d_min").c_str(), kDMin);
  par_dmin.Write();

  TParameter<Long64_t> par_seed((std::string(key) + "_seed").c_str(),
                                static_cast<Long64_t>(kSeed));
  par_seed.Write();

  TParameter<int> par_Ntheta((std::string(key) + "_N_theta").c_str(), N_theta);
  par_Ntheta.Write();

  TParameter<double> par_Pmu((std::string(key) + "_P_mu").c_str(), P_mu);
  par_Pmu.Write();

  TParameter<double> par_cosdphi((std::string(key) + "_cosDeltaPhi").c_str(), kCosDeltaPhi);
  par_cosdphi.Write();

  fout.Close();

  std::cout << "[MakeRMDGridPdf] saved (4D): " << out_filepath << " (key=" << key << ")\n";
  return 0;
}
