// include/p2meg/RMDExtrapolation.h

#ifndef P2MEG_RMD_EXTRAPOLATION_H
#define P2MEG_RMD_EXTRAPOLATION_H

#include <cstdint>

// 前方宣言（ヘッダの依存を減らす）
struct AnalysisWindow4D;
struct DetectorResolutionConst;

// ============================================================
// ESB(Eg) -> AW 外挿係数 k の推定結果
//   k = N_AW / N_ESB
// ここでの N_AW, N_ESB は MC 試行数（カウント）で、
// k と sigma_k は「重み付きカウント（理論重み）」から推定する。
// ============================================================

struct KFactorResult {
  double k = 0.0;
  double sigma_k = 0.0;
  long long n_gen = 0;
  long long n_aw = 0;
  long long n_esb = 0;
};

// ============================================================
// RMD の理論カーネル（RMDSpectrum）を用いた k 推定（MC）
//
// MakeRMDGridPdf の取り扱いに合わせること（重要）:
//  - phi は detres.phi_*_min/max, N_phi_* の格子点
//  - マスク Detector_IsAllowedPhiPairIndex を適用
//  - theta_eg = |phi_e - phi_g|, cosThetaEG = cos(phi_e - phi_g)
//  - Jacobian は掛けない（離散角度を等重みとして扱う）
// ============================================================

KFactorResult RMD_EstimateK_AW_over_ESB_MC(const AnalysisWindow4D& win,
                                           const DetectorResolutionConst& detres,
                                           double Eg_all_min,
                                           double Eg_all_max,
                                           double t_all_min,
                                           double t_all_max,
                                           long long n_gen,
                                           std::uint64_t seed);

#endif // P2MEG_RMD_EXTRAPOLATION_H
