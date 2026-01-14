#ifndef P2MEG_MAKE_RMD_GRID_PDF_H
#define P2MEG_MAKE_RMD_GRID_PDF_H

// ============================================================
// MakeRMDGridPdf
//
// 目的
//  - RMD の完全微分核（偏極込み）
//      RMD_d6B_dEe_dEg_dOmegae_dOmegag
//    を用い、検出器分解能（独立ガウシアン, 定数）で Ee, Eg のみをスメアした
//    畳み込み済み「格子PDF」を作成して ROOT に保存する。
//
// 設計方針（重要）
//  - 角度は「検出器配置の離散化」により扱う（角度スメアなし）。
//    DetectorResolutionConst::N_theta により 0..pi を N_theta 分割し、
//      phi_i = i * pi / N_theta  (i = 0..N_theta)
//    を許される検出器角度とする。
//  - phi_detector_e, phi_detector_g はそれぞれ
//      cosThetaE = cos(phi_ie),  cosThetaG = cos(phi_ig)
//    を使って RMD 式に渡す。
//  - e-γ 相対角は
//      theta_eg = |phi_ie - phi_ig|,
//      cosThetaEG = cos(phi_ie - phi_ig)
//    で与える。
//  - 時間 t は RMD のモデルでは独立（t_true = t_mean）なので、
//    PDF 評価側（RMDGridPdf）で解析的に p_t(t)（窓内正規化ガウシアン）を掛ける。
//    これにより生成・保存を軽くし、評価も高速にする。
//  - 生成側で扱うスメアは Ee, Eg のみ（sigma_Ee, sigma_Eg）。
//
// 外部入力（include/p2meg から参照）
//  - AnalysisWindow.h       : 解析窓（Ee, Eg, t, theta）
//  - DetectorResolution.h   : 分解能（sigma_Ee, sigma_Eg, sigma_t）, t_mean,
//                             および角度分割 N_theta と偏極度 P_mu
//  - RMDSpectrum.h          : RMD_d6B_dEe_dEg_dOmegae_dOmegag
//
// 出力（ROOT）
//  - key 名で「設定込み」の格子PDFを保存する。
//    例：4次元 THnD（Ee, Eg, ie, ig）
//      Ee, Eg : 連続（ビン）
//      ie, ig : 離散（0..N_theta の整数インデックス）
//    ※ 実際の格子表現（THnD/TH2D配列など）は実装側で固定する。
//  - メタ情報（窓・分解能・N_theta・P_mu・ビニング等）も併せて保存する。
//
// 注意
//  - soft photon 発散があるため、呼び出し側（解析窓など）で Eg > Eg_min のカットが必須。
//  - 格子ビニング、生成統計（サンプル数、シード、スメア回数など）は
//    このモジュール内部に固定で持つ（ヘッダには出さない）。
// ============================================================

// 成功: 0、失敗: 非0
int MakeRMDGridPdf(const char* out_filepath, const char* key);

#endif // P2MEG_MAKE_RMD_GRID_PDF_H
