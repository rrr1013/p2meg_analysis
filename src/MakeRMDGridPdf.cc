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
#include "p2meg/Constants.h"
#include "p2meg/RMDSpectrum.h"

//============================================================
// 内部設定
//
// 単位：Ee, Eg は MeV、t は ns、角度は rad
//
// 出力格子（4D）:
//   (Ee, Eg, phi_detector_e, phi_detector_g)
//
// ※ phi_detector_e/g は「偏極軸と各検出器方向のなす角 φ（0..π）」で、
//   離散化は φ を 0..π を N_theta 等分した格子点で行う。
//   ここでは、phi 軸を「最近傍θ格子」に対応する可変ビンにしておき、
//   later の評価側は "最近傍θ" に丸めることで整合を取る。
//============================================================

// ---- 4D格子ビニング（Ee, Eg, phi_e, phi_g）----
static constexpr int kNBins_Ee = 40;
static constexpr int kNBins_Eg = 40;

// ---- 生成統計 ----
static constexpr long kNTruthSamples  = 5000000L; // 真値サンプル数（Ee,Eg を一様）
static constexpr int  kNSmearPerTruth = 500;       // 1真値あたりのエネルギースメア回数
static constexpr unsigned long kSeed  = 20260109UL;

// ---- RMD 理論関数の d_min ----
static constexpr double kDMin = 1e-6;

//============================================================
// 内部補助
//============================================================

static bool IsFinite(double x) { return std::isfinite(x); }

static double Clamp(double x, double lo, double hi) {
  if (x < lo) return lo;
  if (x > hi) return hi;
  return x;
}

static bool IsInsideWindow_EeEg(double Ee, double Eg) {
  if (Ee < analysis_window.Ee_min || Ee > analysis_window.Ee_max) return false;
  if (Eg < analysis_window.Eg_min || Eg > analysis_window.Eg_max) return false;
  return true;
}

static bool IsInsideWindowTheta(double theta_eg) {
  if (theta_eg < analysis_window.theta_min || theta_eg > analysis_window.theta_max) return false;
  return true;
}

// 真値サンプル窓のチェック（Ee, Eg のみ使用）
static bool IsValidTruthWindowEeEg(const AnalysisWindow4D& win) {
  if (!IsFinite(win.Ee_min) || !IsFinite(win.Ee_max)) return false;
  if (!IsFinite(win.Eg_min) || !IsFinite(win.Eg_max)) return false;
  if (!(win.Ee_min < win.Ee_max)) return false;
  if (!(win.Eg_min < win.Eg_max)) return false;
  // soft photon 発散のため Eg_min > 0 を要求
  if (!(win.Eg_min > 0.0)) return false;
  return true;
}

// 解析窓から +/-2sigma だけ広げた真値窓を作り、非物理領域（負のエネルギーなど）を除外する
static AnalysisWindow4D ExpandTruthWindow2Sigma(const AnalysisWindow4D& base_win) {
  AnalysisWindow4D out = base_win;

  const double n_sigma = 2.0;
  const double sigma_Ee = detres.sigma_Ee;
  const double sigma_Eg = detres.sigma_Eg;

  out.Ee_min = base_win.Ee_min - n_sigma * sigma_Ee;
  out.Ee_max = base_win.Ee_max + n_sigma * sigma_Ee;
  out.Eg_min = base_win.Eg_min - n_sigma * sigma_Eg;
  out.Eg_max = base_win.Eg_max + n_sigma * sigma_Eg;

  // 物理的な上限（RMD の運動学: x<=1+r, y<=1-r）
  const double mmu = kMassesPDG.m_mu;
  const double me  = kMassesPDG.m_e;
  const double r = (me * me) / (mmu * mmu);

  const double Ee_phys_max = 0.5 * mmu * (1.0 + r);
  const double Eg_phys_max = 0.5 * mmu * (1.0 - r);

  if (out.Ee_min < 0.0) out.Ee_min = 0.0;
  if (out.Eg_min < 0.0) out.Eg_min = 0.0;
  if (out.Ee_max > Ee_phys_max) out.Ee_max = Ee_phys_max;
  if (out.Eg_max > Eg_phys_max) out.Eg_max = Eg_phys_max;

  // soft photon 発散のため Eg_min は正に保つ
  if (out.Eg_min <= 0.0) out.Eg_min = 1e-6;

  return out;
}

// detres.N_theta を int として使う（型が double でもここで丸める）
static int GetNTheta() {
  // N_theta >= 1 を想定（0 は不正）
  const double x = static_cast<double>(detres.N_theta);
  int N = static_cast<int>(std::lround(x));
  if (N < 1) N = 1;
  return N;
}

// φ= i*pi/N (i=0..N) の格子に対応する phi 軸の可変ビン境界を作る。
// 「最近傍θ格子」になるように、φ の中点 (i+0.5)*Δθ を境界にする。
// 出力: edges サイズ = (N+1)+1 = N+2、x は昇順（0..π）
static std::vector<double> BuildPhiEdgesFromThetaMidpoints(int N_theta) {
  const double dth = pi / static_cast<double>(N_theta);

  std::vector<double> edges;
  edges.resize(static_cast<size_t>(N_theta + 2));

  edges[0] = 0.0;
  for (int i = 1; i <= N_theta; ++i) {
    edges[i] = (static_cast<double>(i) - 0.5) * dth;
  }
  edges[N_theta + 1] = pi;

  for (int k = 1; k < static_cast<int>(edges.size()); ++k) {
    if (edges[k] < edges[k - 1]) edges[k] = edges[k - 1];
    edges[k] = Clamp(edges[k], 0.0, pi);
  }

  return edges;
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

// 4Dヒストを「密度」に変換し、指定した total_mass で正規化する。
// density4 = (C / total_mass) / volume4(bin)
static int ConvertToDensityAndNormalize4(THnD& h, double total_mass) {
  const TAxis* ax0 = h.GetAxis(0);
  const TAxis* ax1 = h.GetAxis(1);
  const TAxis* ax2 = h.GetAxis(2);
  const TAxis* ax3 = h.GetAxis(3);

  const int n0 = ax0->GetNbins();
  const int n1 = ax1->GetNbins();
  const int n2 = ax2->GetNbins();
  const int n3 = ax3->GetNbins();

  if (!(total_mass > 0.0) || !IsFinite(total_mass)) return 1;

  std::vector<int> idx(4, 1);

  for (int i0 = 1; i0 <= n0; ++i0) {
    idx[0] = i0;
    const double w0 = ax0->GetBinWidth(i0);
    for (int i1 = 1; i1 <= n1; ++i1) {
      idx[1] = i1;
      const double w1 = ax1->GetBinWidth(i1);
      for (int i2 = 1; i2 <= n2; ++i2) {
        idx[2] = i2;
        const double w2 = ax2->GetBinWidth(i2);
        for (int i3 = 1; i3 <= n3; ++i3) {
          idx[3] = i3;
          const double w3 = ax3->GetBinWidth(i3);

          const double vol = w0 * w1 * w2 * w3; // [MeV^2 * rad^2]
          if (!(vol > 0.0) || !IsFinite(vol)) continue;

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

static std::string BuildMetaString(long n_theta_ok, long n_wpos,
                                   long n_fill, double total_mass,
                                   const AnalysisWindow4D& truth_win) {
  std::ostringstream oss;
  oss << "MakeRMDGridPdf meta (4D grid, phi only)\n";
  oss << "bins: Ee=" << kNBins_Ee << ", Eg=" << kNBins_Eg
      << ", phi_e=" << (GetNTheta() + 1) << ", phi_g=" << (GetNTheta() + 1) << "\n";
  oss << "truth_samples=" << kNTruthSamples << "\n";
  oss << "smear_per_truth=" << kNSmearPerTruth << "\n";
  oss << "seed=" << kSeed << "\n";
  oss << "d_min=" << kDMin << "\n";

  oss << "truth window: Ee=[" << truth_win.Ee_min << "," << truth_win.Ee_max
      << "] MeV (analysis +/-2sigma, physical clip)\n";
  oss << "truth window: Eg=[" << truth_win.Eg_min << "," << truth_win.Eg_max
      << "] MeV (analysis +/-2sigma, physical clip)\n";
  oss << "window: Ee=[" << analysis_window.Ee_min << "," << analysis_window.Ee_max << "] MeV (observed)\n";
  oss << "window: Eg=[" << analysis_window.Eg_min << "," << analysis_window.Eg_max << "] MeV (observed)\n";
  oss << "window: t=[" << analysis_window.t_min << "," << analysis_window.t_max << "] ns (applied at evaluation)\n";
  oss << "window: theta=[" << analysis_window.theta_min << "," << analysis_window.theta_max << "] rad (theta_eg=|phi_e-phi_g| cut at generation)\n";

  oss << "res: sigma_Ee=" << detres.sigma_Ee << " MeV\n";
  oss << "res: sigma_Eg=" << detres.sigma_Eg << " MeV\n";
  oss << "res: sigma_t=" << detres.sigma_t << " ns (used analytically at evaluation)\n";
  oss << "res: t_mean=" << detres.t_mean << " ns (used analytically at evaluation)\n";
  oss << "res: N_theta=" << GetNTheta() << "\n";
  oss << "res: P_mu=" << detres.P_mu << "\n";

  oss << "theta_window_ok=" << n_theta_ok << "\n";
  oss << "w>0 ok=" << n_wpos << "\n";
  oss << "filled entries=" << n_fill << "\n";
  oss << "raw mass=" << total_mass << "\n";

  oss << "saved keys: <key>\n";
  return oss.str();
}

//============================================================
// 本体
//============================================================

int MakeRMDGridPdf(const char* out_filepath, const char* key) {
  return MakeRMDGridPdfWithTruthWindow(out_filepath, key, analysis_window);
}

int MakeRMDGridPdfWithTruthWindow(const char* out_filepath, const char* key,
                                  const AnalysisWindow4D& analysis_win_for_truth) {
  if (!out_filepath || !key) {
    std::cerr << "[MakeRMDGridPdf] invalid arguments\n";
    return 1;
  }

  // 解析窓と真値窓の基本チェック
  const AnalysisWindow4D truth_win = ExpandTruthWindow2Sigma(analysis_win_for_truth);
  if (!IsValidTruthWindowEeEg(truth_win)) {
    std::cerr << "[MakeRMDGridPdf] truth window (Ee/Eg) is invalid.\n";
    return 2;
  }
  if (!(analysis_window.Eg_min > 0.0)) {
    std::cerr << "[MakeRMDGridPdf] analysis Eg_min must be > 0 to avoid soft photon divergence.\n";
    return 2;
  }
  if (truth_win.Ee_min > analysis_window.Ee_min ||
      truth_win.Ee_max < analysis_window.Ee_max ||
      truth_win.Eg_min > analysis_window.Eg_min ||
      truth_win.Eg_max < analysis_window.Eg_max) {
    std::cerr << "[MakeRMDGridPdf] WARNING: truth window (analysis +/-2sigma, physical clip) is narrower than analysis window.\n";
    std::cerr << "[MakeRMDGridPdf] WARNING: smearing-in contributions may be lost.\n";
  }

  const int N_theta = GetNTheta();
  const double dth = pi / static_cast<double>(N_theta);

  // ---- 4Dヒスト（Ee, Eg, phi_e, phi_g） ----
  const int ndim = 4;
  int nbins[ndim] = {kNBins_Ee, kNBins_Eg, N_theta + 1, N_theta + 1};
  double xmin[ndim] = {analysis_window.Ee_min, analysis_window.Eg_min, 0.0, 0.0};
  double xmax[ndim] = {analysis_window.Ee_max, analysis_window.Eg_max, pi, pi};

  THnD h("rmd_grid_tmp", "RMD grid (phi);Ee;Eg;phi_e;phi_g", ndim, nbins, xmin, xmax);
  h.Sumw2();

  h.GetAxis(0)->SetTitle("Ee [MeV]");
  h.GetAxis(1)->SetTitle("Eg [MeV]");
  h.GetAxis(2)->SetTitle("phi_detector_e [rad]");
  h.GetAxis(3)->SetTitle("phi_detector_g [rad]");

  // phi 軸を「最近傍θ格子」になる可変ビンに変更
  const std::vector<double> phi_edges = BuildPhiEdgesFromThetaMidpoints(N_theta);
  h.GetAxis(2)->Set(N_theta + 1, phi_edges.data());
  h.GetAxis(3)->Set(N_theta + 1, phi_edges.data());

  TRandom3 rng(kSeed);

  const double Ee_true_min = truth_win.Ee_min;
  const double Ee_true_max = truth_win.Ee_max;
  const double Eg_true_min = truth_win.Eg_min;
  const double Eg_true_max = truth_win.Eg_max;

  const double sigma_Ee = detres.sigma_Ee;
  const double sigma_Eg = detres.sigma_Eg;

  if (!(sigma_Ee > 0.0) || !(sigma_Eg > 0.0)) {
    std::cerr << "[MakeRMDGridPdf] sigma_Ee/sigma_Eg must be > 0\n";
    return 5;
  }

  const double Pmu = detres.P_mu;

  long n_theta_ok = 0;
  long n_wpos = 0;
  long n_fill = 0;

  for (long it = 0; it < kNTruthSamples; ++it) {
    const double Ee_true = rng.Uniform(Ee_true_min, Ee_true_max);
    const double Eg_true = rng.Uniform(Eg_true_min, Eg_true_max);

    // 検出器角（離散）：θ = i*pi/N (i=0..N)
    const int ie = static_cast<int>(rng.Integer(static_cast<ULong_t>(N_theta + 1)));
    const int ig = static_cast<int>(rng.Integer(static_cast<ULong_t>(N_theta + 1)));

    const double phi_e = static_cast<double>(ie) * dth;
    const double phi_g = static_cast<double>(ig) * dth;

    const double cosThetaE = std::cos(phi_e);
    const double cosThetaG = std::cos(phi_g);
    const double sinE = std::sin(phi_e);
    const double sinG = std::sin(phi_g);

    const double cosThetaEG = Clamp(std::cos(phi_e - phi_g), -1.0, 1.0);
    const double theta_eg = std::fabs(phi_e - phi_g);

    const bool theta_ok = IsInsideWindowTheta(theta_eg);
    if (theta_ok) ++n_theta_ok;

    // 変数変換の Jacobian（dcos = -sinφ dφ）
    const double w_ang = sinE * sinG;
    if (!(w_ang > 0.0) || !IsFinite(w_ang)) continue;

    double w0 = 0.0;
    if (theta_ok) {
      w0 = RMD_d6B_dEe_dEg_dOmegae_dOmegag(
          Ee_true, Eg_true, cosThetaEG, cosThetaE, cosThetaG, Pmu, kDMin);
      if (w0 > 0.0 && IsFinite(w0)) ++n_wpos;
      else w0 = 0.0;
    }

    if (!(w0 > 0.0)) continue;

    // proposal 補正（角度離散の補正）
    const double w = w0 * w_ang;

    // 同じ真値点から複数回（Ee,Eg のみ）スメアして観測分布を埋める
    for (int is = 0; is < kNSmearPerTruth; ++is) {
      const double Ee_obs = Ee_true + rng.Gaus(0.0, sigma_Ee);
      const double Eg_obs = Eg_true + rng.Gaus(0.0, sigma_Eg);

      if (!IsInsideWindow_EeEg(Ee_obs, Eg_obs)) continue;

      double x[4] = {Ee_obs, Eg_obs, phi_e, phi_g};

      if (w > 0.0) {
        h.Fill(x, w);
        ++n_fill;
      }
    }
  }

  std::cout << "[MakeRMDGridPdf] setting theta-window ok: " << n_theta_ok
            << " / " << kNTruthSamples << "\n";
  std::cout << "[MakeRMDGridPdf] w>0 ok=" << n_wpos << "\n";
  std::cout << "[MakeRMDGridPdf] filled entries (4D)=" << n_fill << "\n";

  // ---- 正規化（∫ pdf dEe dEg dphi_e dphi_g = 1 になるよう密度化）----
  const double total_mass = SumAllBins4(h);

  if (!(total_mass > 0.0) || !IsFinite(total_mass)) {
    std::cerr << "[MakeRMDGridPdf] total mass is not positive (mass=" << total_mass
              << ")\n";
    return 3;
  }

  if (ConvertToDensityAndNormalize4(h, total_mass) != 0) {
    std::cerr << "[MakeRMDGridPdf] normalization failed for grid\n";
    return 3;
  }

  // ---- 保存 ----
  TFile fout(out_filepath, "RECREATE");
  if (fout.IsZombie()) {
    std::cerr << "[MakeRMDGridPdf] cannot open output file: " << out_filepath << "\n";
    return 4;
  }

  h.SetName(key);
  h.Write(key);

  const std::string meta = BuildMetaString(n_theta_ok, n_wpos, n_fill, total_mass, truth_win);
  TNamed meta_obj((std::string(key) + "_meta").c_str(), meta.c_str());
  meta_obj.Write();

  TParameter<double> par_dmin((std::string(key) + "_d_min").c_str(), kDMin);
  par_dmin.Write();

  TParameter<Long64_t> par_seed((std::string(key) + "_seed").c_str(),
                                static_cast<Long64_t>(kSeed));
  par_seed.Write();

  TParameter<int> par_Ntheta((std::string(key) + "_N_theta").c_str(), N_theta);
  par_Ntheta.Write();

  TParameter<double> par_Pmu((std::string(key) + "_P_mu").c_str(), Pmu);
  par_Pmu.Write();

  fout.Close();

  std::cout << "[MakeRMDGridPdf] saved (4D): " << out_filepath
            << " (key=" << key << ")\n";
  return 0;
}
