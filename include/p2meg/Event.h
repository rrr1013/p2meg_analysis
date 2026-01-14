#ifndef P2MEG_EVENT_H
#define P2MEG_EVENT_H

// ============================================================
// p2MEG 事象データ構造
//
// 単位:
//  - Ee, Eg : MeV
//  - t      : ns
//  - theta  : rad
//
// phi_detector_e/g の定義:
//  - phi_detector_e = 偏極軸と e+ 側検出器方向のなす角 [rad]
//  - phi_detector_g = 偏極軸と γ 側検出器方向のなす角 [rad]
//  値域は [0, pi] を想定。
//
// d_e-hat, d_g-hat は運用で定義する（例：ターゲット中心→該当検出器中心）。
// ============================================================

struct Event {
    double Ee;    // [MeV]
    double Eg;    // [MeV]
    double t;     // [ns]
    double theta; // [rad]

    double phi_detector_e; // [rad]
    double phi_detector_g; // [rad]
};

#endif // P2MEG_EVENT_H
