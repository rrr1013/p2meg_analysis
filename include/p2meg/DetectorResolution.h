#ifndef P2MEG_DETECTOR_RESOLUTION_H
#define P2MEG_DETECTOR_RESOLUTION_H

// ============================================================
// p2MEG 分解能モデル
//
// 単位:
//  - Ee, Eg: MeV
//  - t     : ns
//  - theta : rad
//
// t_mean は Δt の平均値（現状は 0、将来ズレる可能性に備えて分離）
// ============================================================

struct DetectorResolutionConst {
    double sigma_Ee;     // [MeV]
    double sigma_Eg;     // [MeV]
    double sigma_t;      // [ns]
    double sigma_theta;  // [rad]
    double t_mean;       // [ns]  Δt の平均値
};

// 既定値（要調整）
// cf. 10°=0.1745 rad, 20°=0.3491 rad
inline constexpr DetectorResolutionConst detres_rmd{
    5.0,    // sigma_Ee [MeV]（仮）
    5.0,    // sigma_Eg [MeV]（仮）
    0.5,    // sigma_t  [ns]  （仮）
    0.26,   // sigma_theta [rad]（仮：15°相当）
    0.0     // t_mean [ns]（現状は 0）
};

#endif // P2MEG_DETECTOR_RESOLUTION_H
