#ifndef P2MEG_EVENT_H
#define P2MEG_EVENT_H

// ============================================================
// p2MEG 事象データ構造
//
// 単位:
//  - Ee, Eg : MeV
//  - t      : ns
//  - phi_detector_* : rad
//
// phi_detector_e/g の定義:
//  - 偏極軸 P-hat を基準とした検出器方向の方位角
//  - e+/gamma 側検出器の代表方向 d_e-hat / d_g-hat を
//    偏極軸と同一平面内で測った角度（単位 [rad]）
//  - 角度の原点は P-hat 方向とし、取り方は運用で統一する
//
// d_e-hat, d_g-hat は運用で定義する（例：ターゲット中心→該当検出器中心）。
// ============================================================

struct Event {
    double Ee;    // [MeV]
    double Eg;    // [MeV]
    double t;     // [ns]

    double phi_detector_e; // [rad]
    double phi_detector_g; // [rad]
};

#endif // P2MEG_EVENT_H
