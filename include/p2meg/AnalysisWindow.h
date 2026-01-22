#ifndef P2MEG_ANALYSIS_WINDOW_H
#define P2MEG_ANALYSIS_WINDOW_H

// ============================================================
// p2MEG 解析窓
//
// 単位:
//  - Ee, Eg: MeV
//  - t     : ns
//  - theta : rad
//
// theta は e と γ のなす角で 0 <= theta <= pi を想定。
// ============================================================

struct AnalysisWindow4D {
    double Ee_min;     // [MeV]
    double Ee_max;     // [MeV]
    double Eg_min;     // [MeV]
    double Eg_max;     // [MeV]
    double t_min;      // [ns]
    double t_max;      // [ns]
    double theta_min;  // [rad]
    double theta_max;  // [rad]
};

inline constexpr AnalysisWindow4D analysis_window{
    30.0, 70.0,    // Ee [MeV]
    30.0, 70.0,    // Eg [MeV]
    -2.0, 2.0,     // t  [ns]
    2, 3.1415926536 // theta [rad]
};

#endif // P2MEG_ANALYSIS_WINDOW_H
