#ifndef P2MEG_MAKE_RMD_GRID_PDF_H
#define P2MEG_MAKE_RMD_GRID_PDF_H

// ============================================================
// MakeRMDGridPdf
//
// 目的
//  - RMD 理論式（RMD_d3B_dEe_dEg_dcos）を用い、
//    検出器分解能（独立ガウシアン, 定数）でスメアした
//    畳み込み済み「格子PDF」を作成して root に保存する。
//
// 設計方針（重要）
//  - 出力する格子は 3D（Ee, Eg, theta）のみ。
//    時間 t は RMD のモデルでは独立（t_true = t_mean）なので、
//    PDF 評価側（RMDGridPdf）で解析的に p_t(t) を掛ける。
//  - これにより生成と保存が軽くなり、評価も高速になる。
//
// 外部入力（include/p2meg から参照）
//  - AnalysisWindow.h : 解析窓（Ee, Eg, t, theta）
//  - DetectorResolution.h : 分解能（sigma_Ee, sigma_Eg, sigma_t, sigma_theta）と t_mean
//
// 出力（root）
//  - key 名で 3D の THnD（Ee, Eg, theta）を保存する。
//  - メタ情報（窓・分解能・ビニング・生成統計など）も併せて保存する。
//
// 注意
//  - 格子ビニング、生成統計（サンプル数、シード、スメア回数など）は
//    このモジュール内部に固定で持つ（ヘッダには出さない）。
//  - 角度は theta（e と γ のなす角, 0<=theta<=pi）を格子変数とし、
//    境界 pi 近傍の歪みを避けるため δ=pi-theta を使った折り返しスメアを用いる。
// ============================================================

// 成功: 0、失敗: 非0
int MakeRMDGridPdf(const char* out_filepath, const char* key);

#endif // P2MEG_MAKE_RMD_GRID_PDF_H
