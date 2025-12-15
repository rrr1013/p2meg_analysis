#ifndef P2MEG_RMD_GRID_PDF_H
#define P2MEG_RMD_GRID_PDF_H

// ============================================================
// RMDGridPdf
//
// 目的
//  - オフラインで作成した「RMD 畳み込み済み 格子PDF」を読み込み、
//    任意の観測値 (Ee, Eg, t, theta) に対して PDF 値を返す。
//
// 入力（Eval）
//  - Ee    : 陽電子エネルギー [MeV]
//  - Eg    : ガンマ線エネルギー [MeV]
//  - t     : 到達時間差 Δt [ns]
//  - theta : e と γ のなす角 θ [rad]（0<=theta<=pi）
//
// 出力（Eval）
//  - 解析窓内なら PDF 値（密度）を返す
//  - 解析窓外なら 0 を返す
//
// 実装方針（重要）
//  - root から読み込む格子は 3D（Ee, Eg, theta）の THnD のみ。
//  - 3D補間（8点）で p3(Ee,Eg,theta) を計算し、時間因子を解析的に掛ける：
//      p(t) = N(t_mean, sigma_t) / ∫_{tmin}^{tmax} N dt
//    ここで t_mean, sigma_t は DetectorResolution.h、tmin/tmax は AnalysisWindow.h から取る。
// ============================================================

// 初期化：root ファイルから格子PDF（3D）を読み込む
//  - filepath : 例 "data/pdf_cache/rmd_grid.root"
//  - key      : 例 "rmd_grid"
// 成功したら true
bool RMDGridPdf_Load(const char* filepath, const char* key);

// 現在ロード済みかどうか（デバッグ・安全用）
bool RMDGridPdf_IsLoaded();

// PDF 評価（ロード済みが前提）
// 解析窓外は 0 を返す
double RMDGridPdf(double Ee, double Eg, double t, double theta);

#endif // P2MEG_RMD_GRID_PDF_H
