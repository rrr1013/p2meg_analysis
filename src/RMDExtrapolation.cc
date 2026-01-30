// src/RMDExtrapolation.cc
//
// ============================================================
// p2MEG RMD: ESB(Eg) -> AW 外挿係数 k の推定（Monte Carlo）
//
// k = N_AW / N_ESB
//  - AW : (Ee, Eg, t, theta_eg) が解析窓内
//  - ESB: (Ee, t, theta_eg) は解析窓内、Eg は「Eg のサイドバンド」
//         (全 Eg 範囲 [Eg_all_min, Eg_all_max] 内で解析窓外)
//
// MakeRMDGridPdf と整合させるための重要点:
//  - 角度は detres.phi_e/g_min/max, N_phi_e/g の格子点のみを許す。
//  - phi_detector_e/g は「極角 θ」と同等に扱い、
//      cosThetaE = cos(phi_e), cosThetaG = cos(phi_g)
//    を RMD 式に入力する。
//  - e-γ 相対角は
//      cosThetaEG = cos(phi_e - phi_g)
//      theta_eg   = |phi_e - phi_g|
//    とする。
//  - 許可領域（マスク）は Detector_IsAllowedPhiPairIndex を適用する。
//  - Jacobian（sinθ など）は掛けない（離散角度を等重みで扱う）。
//  - エネルギーは smear_energy_trandom3_e/g でスメア。
//  - time は窓内代表点 t0 固定（比では本質でない）。
// ============================================================

#include "p2meg/RMDExtrapolation.h"

#include <cmath>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

#include "TRandom3.h"

#include "p2meg/Constants.h"
#include "p2meg/AnalysisWindowUtils.h"
#include "p2meg/DetectorResolution.h"
#include "p2meg/RMDSpectrum.h"

// RMD 理論の d_min（共線領域の数値ガード。物理カットではない）
static constexpr double kDMin = 1e-6;

static inline double SafeTimeInWindow(const AnalysisWindow4D& win) {
  const double t0 = 0.5 * (win.t_min + win.t_max);
  if (!std::isfinite(t0)) return 0.0;
  if (t0 < win.t_min) return win.t_min;
  if (t0 > win.t_max) return win.t_max;
  return t0;
}

// detres の定義どおりの phi グリッド点を作る（i=0..N_phi）
static inline std::vector<double> BuildPhiGrid(double phi_min, double phi_max, int N_phi) {
  std::vector<double> phis;
  if (!Detector_IsPhiRangeValid(phi_min, phi_max, N_phi)) return phis;
  phis.reserve((size_t)N_phi + 1);
  for (int i = 0; i <= N_phi; ++i) {
    phis.push_back(Detector_PhiGridPoint(i, phi_min, phi_max, N_phi));
  }
  return phis;
}

// theta_eg=|phi_e-phi_g| が win に入る (ie, ig) だけ列挙（分散低減）
//  - マスク Detector_IsAllowedPhiPairIndex を適用
static inline std::vector<std::pair<int,int>>
BuildAllowedPairsThetaInWin(const std::vector<double>& phis_e,
                            const std::vector<double>& phis_g,
                            const AnalysisWindow4D& win,
                            const DetectorResolutionConst& detres)
{
  std::vector<std::pair<int,int>> pairs;
  const int ne = (int)phis_e.size();
  const int ng = (int)phis_g.size();
  if (ne <= 0 || ng <= 0) return pairs;

  for (int ie = 0; ie < ne; ++ie) {
    for (int ig = 0; ig < ng; ++ig) {
      if (!Detector_IsAllowedPhiPairIndex(ie, ig, detres)) continue;

      const double pe = phis_e[(size_t)ie];
      const double pg = phis_g[(size_t)ig];
      const double theta_eg = std::fabs(pe - pg);
      if (!std::isfinite(theta_eg)) continue;

      // theta_eg は解析窓内が必須（AW も ESB も共通で要求）
      if (theta_eg < win.theta_min || theta_eg > win.theta_max) continue;

      pairs.emplace_back(ie, ig);
    }
  }
  return pairs;
}

KFactorResult RMD_EstimateK_AW_over_ESB_MC(const AnalysisWindow4D& win,
                                           const DetectorResolutionConst& detres,
                                           double Eg_all_min,
                                           double Eg_all_max,
                                           double t_all_min,
                                           double t_all_max,
                                           long long n_gen,
                                           std::uint64_t seed)
{
  KFactorResult out{};
  out.k = 0.0;
  out.sigma_k = 0.0;
  out.n_gen = n_gen;
  out.n_aw = 0;
  out.n_esb = 0;

  if (!(n_gen > 0)) return out;
  if (!std::isfinite(Eg_all_min) || !std::isfinite(Eg_all_max) || !(Eg_all_max > Eg_all_min)) return out;
  if (!std::isfinite(t_all_min) || !std::isfinite(t_all_max) || !(t_all_max > t_all_min)) return out;

  // time: AW に確実に入る代表点
  const double t0 = SafeTimeInWindow(win);
  if (!AnalysisWindow_InTime(win, t0)) return out;
  if (t0 < t_all_min || t0 > t_all_max) return out;

  // phi グリッド（MakeRMDGridPdf と同じ）
  const std::vector<double> phis_e = BuildPhiGrid(detres.phi_e_min, detres.phi_e_max, detres.N_phi_e);
  const std::vector<double> phis_g = BuildPhiGrid(detres.phi_g_min, detres.phi_g_max, detres.N_phi_g);
  if (phis_e.empty() || phis_g.empty()) return out;

  // 候補ペア（マスク + theta 窓）
  const auto pairs = BuildAllowedPairsThetaInWin(phis_e, phis_g, win, detres);
  if (pairs.empty()) return out;

  // truth エネルギーの提案分布：一様（重要度サンプリングの規格化因子は比で相殺）
  const double Emax_phys = 0.5 * kMassesPDG.m_mu;
  const double Ee_true_min = 0.0;
  const double Ee_true_max = Emax_phys;
  const double Eg_true_min = 0.0;
  const double Eg_true_max = Emax_phys;

  TRandom3 rng((ULong_t)seed);

  double sum_aw = 0.0, sum_esb = 0.0;
  double sum2_aw = 0.0, sum2_esb = 0.0;

  for (long long i = 0; i < n_gen; ++i) {
    // (ie, ig) を候補集合から一様に選ぶ（MakeRMDGridPdf の「離散角度を等重み」に対応）
    const int ip = (int)rng.Integer((ULong_t)pairs.size());
    const int ie = pairs[(size_t)ip].first;
    const int ig = pairs[(size_t)ip].second;

    const double phi_e = phis_e[(size_t)ie];
    const double phi_g = phis_g[(size_t)ig];

    const double theta_eg = std::fabs(phi_e - phi_g);
    if (!std::isfinite(theta_eg)) continue;

    // 理論側に渡す cos（設計どおり）
    const double cos_theta_e  = std::cos(phi_e);
    const double cos_theta_g  = std::cos(phi_g);
    const double cos_theta_eg = std::cos(phi_e - phi_g);

    // truth energies
    const double Ee_true = Ee_true_min + (Ee_true_max - Ee_true_min) * rng.Uniform();
    const double Eg_true = Eg_true_min + (Eg_true_max - Eg_true_min) * rng.Uniform();

    // 理論重み
    const double w = RMD_d6B_dEe_dEg_dOmegae_dOmegag(
        Ee_true, Eg_true,
        cos_theta_eg, cos_theta_e, cos_theta_g,
        detres.P_mu,
        kDMin
    );
    if (!(w > 0.0) || !std::isfinite(w)) continue;

    // 観測エネルギー（DetectorResolution の応答でスメア）
    const double Ee_obs = smear_energy_trandom3_e(rng, Ee_true);
    const double Eg_obs = smear_energy_trandom3_g(rng, Eg_true);

    // AW 判定（4D）
    if (AnalysisWindow_In4D(win, Ee_obs, Eg_obs, t0, theta_eg)) {
      ++out.n_aw;
      sum_aw  += w;
      sum2_aw += w * w;
      continue;
    }

    // ESB 判定（Eg だけサイドバンド、他は AW 内）
    if (Ee_obs < win.Ee_min || Ee_obs > win.Ee_max) continue;
    if (theta_eg < win.theta_min || theta_eg > win.theta_max) continue;
    if (!AnalysisWindow_InEnergySidebandEg(win, Eg_obs, Eg_all_min, Eg_all_max)) continue;

    ++out.n_esb;
    sum_esb  += w;
    sum2_esb += w * w;
  }

  if (!(sum_esb > 0.0) || !std::isfinite(sum_esb)) return out;

  const double k_hat = sum_aw / sum_esb;
  if (!std::isfinite(k_hat) || k_hat < 0.0) return out;

  // 誤差（デルタ法）。排他性：各試行で AW と ESB は同時に立たない。
  const double N = (double)n_gen;
  const double mx  = sum_aw / N;
  const double my  = sum_esb / N;
  const double ex2 = sum2_aw / N;
  const double ey2 = sum2_esb / N;

  double varX = ex2 - mx * mx;
  double varY = ey2 - my * my;
  if (varX < 0.0) varX = 0.0;
  if (varY < 0.0) varY = 0.0;

  const double varS1 = N * varX;
  const double varS2 = N * varY;
  const double covS1S2 = N * (0.0 - mx * my); // E[XY]=0 -> Cov=-E[X]E[Y]

  const double s1 = sum_aw;
  const double s2 = sum_esb;

  double vark = 0.0;
  if (s2 > 0.0) {
    const double s2_2 = s2 * s2;
    const double s2_3 = s2_2 * s2;
    const double s2_4 = s2_2 * s2_2;
    vark = varS1 / s2_2 + (s1 * s1) * varS2 / s2_4 - 2.0 * s1 * covS1S2 / s2_3;
  }
  if (vark < 0.0) vark = 0.0;

  out.k = k_hat;
  out.sigma_k = std::sqrt(vark);
  return out;
}
