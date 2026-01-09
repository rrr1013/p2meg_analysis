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
//   (Ee, Eg, cos_detector_e, cos_detector_g)
//
// ※ cos_detector_e/g は「偏極軸と各検出器方向のなす角の cos」だが、
//   離散化は θ を 0..π を N_theta 等分した格子点で行う。
//   ここでは、cos 軸を「最近傍θ格子」に対応する可変ビンにしておき、
//   later の評価側は "最近傍θ" に丸めることで整合を取る。
//============================================================

// ---- 4D格子ビニング（Ee, Eg, cos_e, cos_g）----
static constexpr int kNBins_Ee = 40;
static constexpr int kNBins_Eg = 40;

// ---- 生成統計 ----
static constexpr long kNTruthSamples  = 1000000L; // 真値サンプル数（Ee,Eg を一様）
static constexpr int  kNSmearPerTruth = 50;       // 1真値あたりのエネルギースメア回数
static constexpr unsigned long kSeed  = 20260109UL;

// ---- RMD 理論関数の d_min ----
static constexpr double kDMin = 1e-6;

// ---- 正規化モード切替（手で書き換える） ----
// false: plus-only（+枝だけ fill & 正規化、-枝は 0 のまま保存）
// true : both（+枝と-枝を fill、(p+ + p-) で正規化）
static constexpr bool kNormalizeBothBranches = true;

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

static bool IsInsideWindowTheta(double theta) {
  if (theta < analysis_window.theta_min || theta > analysis_window.theta_max) return false;
  return true;
}

// detres.N_theta を int として使う（型が double でもここで丸める）
static int GetNTheta() {
  // N_theta >= 1 を想定（0 は不正）
  const double x = static_cast<double>(detres.N_theta);
  int N = static_cast<int>(std::lround(x));
  if (N < 1) N = 1;
  return N;
}

// θ= i*pi/N (i=0..N) の格子に対応する cos 軸の可変ビン境界を作る。
// 「最近傍θ格子」になるように、θ の中点 (i+0.5)*Δθ を境界にする。
// 出力: edges サイズ = (N+1)+1 = N+2、x は昇順（-1..+1）
static std::vector<double> BuildCosEdgesFromThetaMidpoints(int N_theta) {
  const double dth = pi / static_cast<double>(N_theta);

  std::vector<double> edges;
  edges.resize(static_cast<size_t>(N_theta + 2));

  // 最小端（θ=π）
  edges[0] = -1.0;

  // 中間境界：θ = (i+0.5)*Δθ, i=0..N-1
  // cos は θ 増加で減少するので、昇順にするため θ を大→小の順で詰める
  for (int j = 1; j <= N_theta; ++j) {
    const int i = N_theta - j; // i = N-1, N-2, ..., 0
    const double th = (static_cast<double>(i) + 0.5) * dth;
    edges[j] = std::cos(th);
  }

  // 最大端（θ=0）
  edges[N_theta + 1] = 1.0;

  // 数値誤差ガード（単調性を壊さないようにクリップ）
  for (int k = 1; k < static_cast<int>(edges.size()); ++k) {
    if (edges[k] < edges[k - 1]) edges[k] = edges[k - 1];
    edges[k] = Clamp(edges[k], -1.0, 1.0);
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

          const double vol = w0 * w1 * w2 * w3; // [MeV^2 * (cos)^2]
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

static std::string BuildMetaString(long n_theta_ok, long n_wpos_p, long n_wpos_m,
                                   long n_fill_p, long n_fill_m,
                                   double total_mass_p, double total_mass_m) {
  std::ostringstream oss;
  oss << "MakeRMDGridPdf meta (4D grid, two branches)\n";
  oss << "bins: Ee=" << kNBins_Ee << ", Eg=" << kNBins_Eg
      << ", cos_e=" << (GetNTheta() + 1) << ", cos_g=" << (GetNTheta() + 1) << "\n";
  oss << "truth_samples=" << kNTruthSamples << "\n";
  oss << "smear_per_truth=" << kNSmearPerTruth << "\n";
  oss << "seed=" << kSeed << "\n";
  oss << "d_min=" << kDMin << "\n";
  oss << "normalize_both_branches=" << (kNormalizeBothBranches ? "true" : "false") << "\n";

  oss << "window: Ee=[" << analysis_window.Ee_min << "," << analysis_window.Ee_max << "] MeV\n";
  oss << "window: Eg=[" << analysis_window.Eg_min << "," << analysis_window.Eg_max << "] MeV\n";
  oss << "window: t=[" << analysis_window.t_min << "," << analysis_window.t_max << "] ns (applied at evaluation)\n";
  oss << "window: theta=[" << analysis_window.theta_min << "," << analysis_window.theta_max << "] rad (branch-dependent cut at generation)\n";

  oss << "res: sigma_Ee=" << detres.sigma_Ee << " MeV\n";
  oss << "res: sigma_Eg=" << detres.sigma_Eg << " MeV\n";
  oss << "res: sigma_t=" << detres.sigma_t << " ns (used analytically at evaluation)\n";
  oss << "res: t_mean=" << detres.t_mean << " ns (used analytically at evaluation)\n";
  oss << "res: N_theta=" << GetNTheta() << "\n";
  oss << "res: P_mu=" << detres.P_mu << "\n";

  oss << "theta_window_ok (any branch)=" << n_theta_ok << "\n";
  oss << "w>0 ok: plus=" << n_wpos_p << ", minus=" << n_wpos_m << "\n";
  oss << "filled entries: plus=" << n_fill_p << ", minus=" << n_fill_m << "\n";
  oss << "raw mass: plus=" << total_mass_p << ", minus=" << total_mass_m << "\n";

  oss << "saved keys: <key>_p and <key>_m\n";
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

  // soft photon 発散対策：Eg_min > 0
  if (!(analysis_window.Eg_min > 0.0)) {
    std::cerr << "[MakeRMDGridPdf] Eg_min must be > 0 to avoid soft photon divergence.\n";
    return 2;
  }

  const int N_theta = GetNTheta();
  const double dth = pi / static_cast<double>(N_theta);

  // ---- 4Dヒスト（Ee, Eg, cos_e, cos_g）: +枝 と -枝 ----
  const int ndim = 4;
  int nbins[ndim] = {kNBins_Ee, kNBins_Eg, N_theta + 1, N_theta + 1};
  double xmin[ndim] = {analysis_window.Ee_min, analysis_window.Eg_min, -1.0, -1.0};
  double xmax[ndim] = {analysis_window.Ee_max, analysis_window.Eg_max, +1.0, +1.0};

  THnD hP("rmd_grid_plus_tmp",  "RMD grid (+ branch);Ee;Eg;cos_e;cos_g",  ndim, nbins, xmin, xmax);
  THnD hM("rmd_grid_minus_tmp", "RMD grid (- branch);Ee;Eg;cos_e;cos_g",  ndim, nbins, xmin, xmax);
  hP.Sumw2();
  hM.Sumw2();

  hP.GetAxis(0)->SetTitle("Ee [MeV]");
  hP.GetAxis(1)->SetTitle("Eg [MeV]");
  hP.GetAxis(2)->SetTitle("cos_detector_e");
  hP.GetAxis(3)->SetTitle("cos_detector_g");

  hM.GetAxis(0)->SetTitle("Ee [MeV]");
  hM.GetAxis(1)->SetTitle("Eg [MeV]");
  hM.GetAxis(2)->SetTitle("cos_detector_e");
  hM.GetAxis(3)->SetTitle("cos_detector_g");

  // cos 軸を「最近傍θ格子」になる可変ビンに変更
  const std::vector<double> cos_edges = BuildCosEdgesFromThetaMidpoints(N_theta);
  hP.GetAxis(2)->Set(N_theta + 1, cos_edges.data());
  hP.GetAxis(3)->Set(N_theta + 1, cos_edges.data());
  hM.GetAxis(2)->Set(N_theta + 1, cos_edges.data());
  hM.GetAxis(3)->Set(N_theta + 1, cos_edges.data());

  TRandom3 rng(kSeed);

  const double Ee_min = analysis_window.Ee_min;
  const double Ee_max = analysis_window.Ee_max;
  const double Eg_min = analysis_window.Eg_min;
  const double Eg_max = analysis_window.Eg_max;

  const double sigma_Ee = detres.sigma_Ee;
  const double sigma_Eg = detres.sigma_Eg;

  if (!(sigma_Ee > 0.0) || !(sigma_Eg > 0.0)) {
    std::cerr << "[MakeRMDGridPdf] sigma_Ee/sigma_Eg must be > 0\n";
    return 5;
  }

  const double Pmu = detres.P_mu;

  long n_theta_ok = 0;
  long n_wpos_p = 0;
  long n_wpos_m = 0;
  long n_fill_p = 0;
  long n_fill_m = 0;

  for (long it = 0; it < kNTruthSamples; ++it) {
    const double Ee_true = rng.Uniform(Ee_min, Ee_max);
    const double Eg_true = rng.Uniform(Eg_min, Eg_max);

    // 検出器角（離散）：θ = i*pi/N (i=0..N)
    const int ie = static_cast<int>(rng.Integer(static_cast<ULong_t>(N_theta + 1)));
    const int ig = static_cast<int>(rng.Integer(static_cast<ULong_t>(N_theta + 1)));

    const double th_e = static_cast<double>(ie) * dth;
    const double th_g = static_cast<double>(ig) * dth;

    const double cosE = std::cos(th_e);
    const double cosG = std::cos(th_g);
    const double sinE = std::sin(th_e);
    const double sinG = std::sin(th_g);

    // 枝ごとの cos(theta_eγ)
    // cosΔφ=+1 -> cos(θe-θg), cosΔφ=-1 -> cos(θe+θg)
    const double cosEG_p = Clamp(cosE * cosG + sinE * sinG, -1.0, 1.0);
    const double cosEG_m = Clamp(cosE * cosG - sinE * sinG, -1.0, 1.0);

    const double thEG_p = std::acos(cosEG_p);
    const double thEG_m = std::acos(cosEG_m);

    bool any_theta_ok = false;
    const bool theta_ok_p = IsInsideWindowTheta(thEG_p);
    const bool theta_ok_m = IsInsideWindowTheta(thEG_m);
    if (theta_ok_p || theta_ok_m) any_theta_ok = true;
    if (any_theta_ok) ++n_theta_ok;

    // bin 幅（cos 軸の可変ビン）を proposal 補正として使う
    // ここでは「θ を一様にサンプル」しているため、cos 測度での積分を近似するには
    // bin 幅 Δcos を掛ける（旧コードの cos proposal 補正と同じ発想）。
    const int bin_ce = hP.GetAxis(2)->FindBin(cosE);
    const int bin_cg = hP.GetAxis(3)->FindBin(cosG);
    const double dcos_e = hP.GetAxis(2)->GetBinWidth(bin_ce);
    const double dcos_g = hP.GetAxis(3)->GetBinWidth(bin_cg);

    if (!(dcos_e > 0.0) || !(dcos_g > 0.0) || !IsFinite(dcos_e) || !IsFinite(dcos_g)) continue;
    const double w_ang = dcos_e * dcos_g;

    // ---- +枝の重み ----
    double wP0 = 0.0;
    if (theta_ok_p) {
      wP0 = RMD_d6B_dEe_dEg_dOmegae_dOmegag(Ee_true, Eg_true, cosEG_p, cosE, cosG, Pmu, kDMin);
      if (wP0 > 0.0 && IsFinite(wP0)) ++n_wpos_p;
      else wP0 = 0.0;
    }

    // ---- -枝の重み ----
    double wM0 = 0.0;
    if (kNormalizeBothBranches && theta_ok_m) {
      wM0 = RMD_d6B_dEe_dEg_dOmegae_dOmegag(Ee_true, Eg_true, cosEG_m, cosE, cosG, Pmu, kDMin);
      if (wM0 > 0.0 && IsFinite(wM0)) ++n_wpos_m;
      else wM0 = 0.0;
    }

    if (!(wP0 > 0.0) && !(wM0 > 0.0)) continue;

    // proposal 補正（角度離散の補正）
    const double wP = wP0 * w_ang;
    const double wM = wM0 * w_ang;

    // 同じ真値点から複数回（Ee,Eg のみ）スメアして観測分布を埋める
    for (int is = 0; is < kNSmearPerTruth; ++is) {
      const double Ee_obs = Ee_true + rng.Gaus(0.0, sigma_Ee);
      const double Eg_obs = Eg_true + rng.Gaus(0.0, sigma_Eg);

      if (!IsInsideWindow_EeEg(Ee_obs, Eg_obs)) continue;

      double x[4] = {Ee_obs, Eg_obs, cosE, cosG};

      if (wP > 0.0) {
        hP.Fill(x, wP);
        ++n_fill_p;
      }
      if (kNormalizeBothBranches && wM > 0.0) {
        hM.Fill(x, wM);
        ++n_fill_m;
      }
    }
  }

  std::cout << "[MakeRMDGridPdf] setting theta-window ok: " << n_theta_ok
            << " / " << kNTruthSamples << "\n";
  std::cout << "[MakeRMDGridPdf] w>0 ok: plus=" << n_wpos_p
            << " minus=" << n_wpos_m << "\n";
  std::cout << "[MakeRMDGridPdf] filled entries (4D): plus=" << n_fill_p
            << " minus=" << n_fill_m << "\n";

  // ---- 正規化（∫ pdf dEe dEg dcos_e dcos_g = 1 になるよう密度化）----
  const double mass_p = SumAllBins4(hP);
  const double mass_m = SumAllBins4(hM);

  double total_mass = 0.0;
  if (kNormalizeBothBranches) {
    total_mass = mass_p + mass_m;
  } else {
    total_mass = mass_p; // plus-only
  }

  if (!(total_mass > 0.0) || !IsFinite(total_mass)) {
    std::cerr << "[MakeRMDGridPdf] total mass is not positive (mass_p=" << mass_p
              << ", mass_m=" << mass_m << ")\n";
    return 3;
  }

  if (ConvertToDensityAndNormalize4(hP, total_mass) != 0) {
    std::cerr << "[MakeRMDGridPdf] normalization failed for plus grid\n";
    return 3;
  }
  if (kNormalizeBothBranches) {
    if (ConvertToDensityAndNormalize4(hM, total_mass) != 0) {
      std::cerr << "[MakeRMDGridPdf] normalization failed for minus grid\n";
      return 3;
    }
  } else {
    // plus-only のとき -枝は空のまま保存（評価側で読めるようにキーは用意）
    // （hM は未正規化の 0 のままでも問題なし）
  }

  // ---- 保存 ----
  TFile fout(out_filepath, "RECREATE");
  if (fout.IsZombie()) {
    std::cerr << "[MakeRMDGridPdf] cannot open output file: " << out_filepath << "\n";
    return 4;
  }

  const std::string key_p = std::string(key) + "_p";
  const std::string key_m = std::string(key) + "_m";

  hP.SetName(key_p.c_str());
  hP.Write(key_p.c_str());

  hM.SetName(key_m.c_str());
  hM.Write(key_m.c_str());

  const std::string meta = BuildMetaString(n_theta_ok, n_wpos_p, n_wpos_m,
                                           n_fill_p, n_fill_m,
                                           mass_p, mass_m);
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
            << " (keys=" << key_p << "," << key_m << ")\n";
  return 0;
}
