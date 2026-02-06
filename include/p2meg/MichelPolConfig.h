// include/p2meg/MichelPolConfig.h
#ifndef P2MEG_MICHEL_POL_CONFIG_H
#define P2MEG_MICHEL_POL_CONFIG_H

// ============================================================
// Michel 偏極測定（別実験）設定
//
// 目的:
//  - Michel 崩壊 e+ のエネルギースペクトル非対称から P_mu を推定するための
//    設定パラメータを一箇所に集約する。
//
// 単位:
//  - Ee: MeV
//
// 角度の扱い:
//  - 本測定では「対向配置」2方向を用い、cosθ は ±|cosθ| として扱う。
//  - |cosθ| はここで手動編集して変更できるようにする（初期値 0.5 など）。
//
// ビニング:
//  - Ee_min..Ee_max を nbins_Ee 等分した等間隔ビンを使用する。
//
// フィット範囲:
//  - fit_Ee_min..fit_Ee_max に入るビンのみを使用する。
//
// ゼロカウント対策:
//  - cross-ratio の各ビンで a,b,c,d のいずれかが 0 の場合は除外するのが基本。
// ============================================================

struct MichelPolConfig {
    // --- 角度設定 ---
    double cos_theta_abs;  // |cosθ|（対向2方向で cosθ = ±|cosθ| とする）

    // --- ヒスト作成範囲 ---
    double Ee_min;         // [MeV]
    double Ee_max;         // [MeV]
    int    nbins_Ee;       // Ee のビン数（等間隔）

    // --- フィットに使う範囲（ヒスト範囲の部分区間） ---
    double fit_Ee_min;     // [MeV]
    double fit_Ee_max;     // [MeV]

    // --- ビン除外条件 ---
    int    min_counts_each; // cross-ratio の a,b,c,d に要求する最小カウント（通常 1）
};

// 初期値（必要に応じて手動で書き換える）
inline constexpr MichelPolConfig michel_pol_config{
    0.86602,   // cos_theta_abs  (初期: |cosθ|=1/2)
    0.0,   // Ee_min [MeV]
    60.0,  // Ee_max [MeV]
    120,   // nbins_Ee
    10.0,  // fit_Ee_min [MeV]
    55.0,  // fit_Ee_max [MeV]
    1      // min_counts_each
};

#endif // P2MEG_MICHEL_POL_CONFIG_H
