// include/p2meg/MichelPolFit.h
#ifndef P2MEG_MICHEL_POL_FIT_H
#define P2MEG_MICHEL_POL_FIT_H

#include <string>

#include "p2meg/MichelPolConfig.h"

// ============================================================
// Michel 偏極測定（別実験）: P_mu 推定 API
//
// 概要（案A）:
//  - 4データセット（入替え運用）から energy-binned cross-ratio 非対称 A_i を作る
//  - 理論テンプレートから K_i = |cosθ| (V_i/U_i) を作る
//  - A_i ≈ P_mu * K_i を重み付き最小二乗で 1パラメータフィット（解析解）
//
// 入力4ファイルの意味（固定）:
//  - path_A_plus : Run#1 で Detector A を「+側（cosθ>0）」に置いたときの Ee-only
//  - path_A_minus: Run#2 で Detector A を「-側（cosθ<0）」に置いたときの Ee-only
//  - path_B_plus : Run#2 で Detector B を「+側（cosθ>0）」に置いたときの Ee-only
//  - path_B_minus: Run#1 で Detector B を「-側（cosθ<0）」に置いたときの Ee-only
//
// 角度設定:
//  - |cosθ| は michel_pol_config.cos_theta_abs を使用（手動で変更）
//
// エラーハンドリング方針:
//  - 不正入力やフィット不能（有効ビン不足等）の場合は status!=0 を返す
// ============================================================

struct MichelPolFitResult {
    int    status;       // 0: success, !=0: failure
    double P_mu_hat;     // 推定 P_mu
    double err_P_mu;     // 1σ 誤差（重み付き最小二乗の解析式）
    double chi2;         // χ^2
    int    ndf;          // 自由度
    int    n_bins_used;  // フィットに使ったビン数
    long long n_read_total;    // 読み込んだ総事象数（4ファイル合計）
    long long n_skipped_total; // スキップした行数（4ファイル合計）
};

// 4ファイルから P_mu を推定する（案A）。
MichelPolFitResult FitMichelPolarizationFrom4Files(const std::string& path_A_plus,
                                                   const std::string& path_A_minus,
                                                   const std::string& path_B_plus,
                                                   const std::string& path_B_minus);

#endif // P2MEG_MICHEL_POL_FIT_H
