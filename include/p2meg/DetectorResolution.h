#ifndef P2MEG_DETECTOR_RESOLUTION_H
#define P2MEG_DETECTOR_RESOLUTION_H

// ============================================================
// p2MEG 分解能モデル（現状の簡略版）
//
// 単位:
//  - Ee, Eg: MeV
//  - t     : ns
//
// 角度 theta の扱い:
//  - 測定設定に合わせて、角度は離散化する
//      theta_i = i * pi / N_theta   (i = 0..N_theta)
//    すなわち、0 から pi までを N_theta 分割した (N_theta+1) 点のみを許す
//  - 崩壊後の e, γ の散乱は無視し、離散化後の角度スメアはかけない
//
// t_mean は Δt の平均値（現状は 0 以外にもなり得るので分離）
// ============================================================

struct DetectorResolutionConst {
    double sigma_Ee;  // [MeV]
    double sigma_Eg;  // [MeV]
    double sigma_t;   // [ns]
    int    N_theta;   // 角度分割数（theta_i = i*pi/N_theta, i=0..N_theta）
    double t_mean;    // [ns]  Δt の平均値
    double P_mu;      // muon polarization (signed, [-1,1])
};

inline constexpr DetectorResolutionConst detres{
    9.264,   // sigma_Ee [MeV]
    9.908,   // sigma_Eg [MeV]
    0.1561,  // sigma_t  [ns]
    18,      // N_theta  （例：0..pi を 18 分割 → 19 点）
    -0.1479, // t_mean [ns]
    -0.8     // P_mu
};

#endif // P2MEG_DETECTOR_RESOLUTION_H
