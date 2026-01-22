// include/p2meg/RMDGridPdf.h
#ifndef P2MEG_RMD_GRID_PDF_H
#define P2MEG_RMD_GRID_PDF_H

// ============================================================
// RMDGridPdf
//
// 目的
//  - オフラインで作成した「RMD 畳み込み済み 格子PDF」を読み込み、
//    任意の観測値 (Ee, Eg, t, phi_detector_e, phi_detector_g) に対して
//    PDF 値を返す。
//  - 格子は 4D（Ee, Eg, phi_detector_e, phi_detector_g）で保存しておき、
//    評価時に theta_eg = |phi_e - phi_g| を作って解析窓カットを行う。
//
// 入力（Eval）
//  - Ee            : 陽電子エネルギー [MeV]
//  - Eg            : ガンマ線エネルギー [MeV]
//  - t             : 到達時間差 Δt [ns]
//  - phi_detector_e: 偏極軸と e 側検出器代表方向の角度 φ_e [rad]
//  - phi_detector_g: 偏極軸と γ 側検出器代表方向の角度 φ_g [rad]
//
// 出力（Eval）
//  - 解析窓内なら PDF 値（密度）を返す
//  - 解析窓外・未ロード・不正入力・phi マスク外なら 0 を返す
//
// 実装方針（重要）
//  - root から読み込む格子は 4D（Ee, Eg, phi_e, phi_g）の THnD を 1つ。
//  - Ee, Eg は 2D（4点）補間、phi 軸は「最近傍θ格子」ビンに落として評価。
//  - 時間因子 p(t) は解析的に掛ける（窓内正規化ガウシアン）。
// ============================================================

// 初期化：root ファイルから格子PDF（4D）を読み込む
//  - filepath : 例 "data/pdf_cache/rmd_grid.root"
//  - key      : 例 "rmd_grid"
// 成功したら true
bool RMDGridPdf_Load(const char* filepath, const char* key);

// 現在ロード済みかどうか（デバッグ・安全用）
bool RMDGridPdf_IsLoaded();

// PDF 評価（ロード済みが前提）
// 解析窓外は 0 を返す
double RMDGridPdf(double Ee, double Eg, double t,
                  double phi_detector_e, double phi_detector_g);

#endif // P2MEG_RMD_GRID_PDF_H
