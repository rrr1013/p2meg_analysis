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
// cos_detector_e/g の定義:
//  - cos_detector_e = P-hat · d_e-hat（e+ 側検出器の代表方向）
//  - cos_detector_g = P-hat · d_g-hat（gamma 側検出器の代表方向）
//  値域は [-1, 1] を想定（無次元）。
//
// d_e-hat, d_g-hat は運用で定義する（例：ターゲット中心→該当検出器中心）。
// ============================================================

struct Event {
    double Ee;    // [MeV]
    double Eg;    // [MeV]
    double t;     // [ns]
    double theta; // [rad]

    double cos_detector_e; // [-1,1]
    double cos_detector_g; // [-1,1]
};

#endif // P2MEG_EVENT_H
