#ifndef P2MEG_RMD_SPECTRUM_H
#define P2MEG_RMD_SPECTRUM_H

//============================================================
// RMDSpectrum.h
//
// 停止ミューオン静止系での Radiative Muon Decay (RMD)
//   mu+ -> e+ nu_e nu_mu_bar gamma
// の「三重微分分岐比（shape の核）」を返す関数群。
//
// ここでいう shape の核とは
//   d^3B / (dx dy dcosθ_{eγ})
// をそのまま返す（解析ウィンドウでの正規化はしない）ものです。
//
// 【重要：このコードがやっていないこと】
//  - 偏極項（G, H）は入れていません（F 項のみ）。
//    (x, y, cosθ_{eγ}) に射影して形状だけを見る場合、受容が十分対称なら
//    偏極項は平均で消え得ますが、受容・カットが非対称なら残り得ます。
//    最終的には MC で P=-1 と P=0 を比較して shape への影響が無視できるか確認してください。
//  - QED の高次補正は入れていません（ツリーレベル SM(V-A)）。
//
// 【使用時の注意】
// 1) 赤外発散（soft photon）
//    y -> 0（Eγ -> 0）で 1/y により発散します。
//    必ず Eγ > Eγ_min のカット（解析ウィンドウ）を入れてください。
//    規格化（ウィンドウ積分値）はこのカットに強く依存します。
// 2) 共線領域（collinear）
//    d = 1 - β cosθ_{eγ} -> 0 で 1/d, 1/d^2 が大きくなります。
//    数値安定化のため d_min を設けています（デフォルト 1e-6）。
//    ただし、これはあくまで計算のガードであり、物理カットではありません。
// 3) 解析で PDF として使う場合
//    この関数が返すのは「核」なので、解析ウィンドウ W 内で積分して
//      R = (1/N_W) * (d^3B/dx dy dcos)
//    のように正規化してから尤度の shape-PDF として使ってください。
//============================================================

/// 三重微分分岐比（核）: d^3B / (dx dy dcosθ_{eγ})
/// 入力:
///   x = 2 Ee / m_mu
///   y = 2 Eγ / m_mu
///   cosTheta = cos(theta_{eγ})
///   d_min: 数値安定化のための下限（d=1-βcosθ がこれ未満なら d=d_min に置き換え）
/// 出力:
///   d^3B/dx dy dcosθ（無次元）
/// 注意:
///   y->0 で発散するので、呼び出し側で必ず y>y_min（Eγ>閾値）を課してください。
double RMD_d3B_dxdy_dcos(double x, double y, double cosTheta, double d_min = 1e-6);

/// 便利関数: エネルギー変数での三重微分分岐比（核）
///   d^3B / (dEe dEγ dcosθ_{eγ})
/// 入力は MeV。
double RMD_d3B_dEe_dEg_dcos(double Ee_MeV, double Eg_MeV, double cosTheta, double d_min = 1e-6);

#endif // P2MEG_RMD_SPECTRUM_H
