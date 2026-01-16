#ifndef P2MEG_MAKE_RMD_GRID_PDF_H
#define P2MEG_MAKE_RMD_GRID_PDF_H

// ============================================================
// MakeRMDGridPdf
//
// 目的
//  - RMD の完全微分核（偏極込み）
//      RMD_d6B_dEe_dEg_dOmegae_dOmegag
//    を用い、検出器分解能 energy_response_shape_e/g(E_res,E_true) で Ee, Eg のみをスメアした
//    畳み込み済み「格子PDF」を作成して ROOT に保存する。
//
// 設計方針（重要）
//  - 角度は「検出器配置の離散化」により扱う（角度スメアなし）。
//    DetectorResolutionConst::N_theta により 0..pi を N_theta 分割し、
//      phi_i = i * pi / N_theta  (i = 0..N_theta)
//    を許される検出器角度とする。
//  - phi_detector_e, phi_detector_g はそれぞれ
//      cosThetaE = cos(phi_e),  cosThetaG = cos(phi_g)
//    を介して RMD 式に入力する。
//  - e-γ 相対角は
//      cosThetaEG = cos(phi_e - phi_g)
//    で与える（theta_eg = |phi_e - phi_g|）。
//  - 時間 t は RMD のモデルでは独立（t_true = t_mean）なので、
//    PDF 評価側（RMDGridPdf）で解析的に p_t(t)（窓内正規化ガウシアン）を掛ける。
//    これにより生成・保存を軽くし、評価も高速にする。
//  - 生成側で扱うスメアは Ee, Eg のみ（energy_response_shape_e/g に従う）。
//
// 外部入力（include/p2meg から参照）
//  - AnalysisWindow.h       : 解析窓（Ee, Eg, t, theta）
//  - DetectorResolution.h   : 分解能（energy_response_shape_e/g, sigma_t）, t_mean,
//                             および角度分割 N_theta と偏極度 P_mu
//  - RMDSpectrum.h          : RMD_d6B_dEe_dEg_dOmegae_dOmegag
//
// 出力（ROOT）
//  - key 名で「設定込み」の格子PDFを保存する。
//    例：4次元 THnD（Ee, Eg, phi_e, phi_g）
//      Ee, Eg : 連続（ビン）
//      phi_e, phi_g : 離散（0..pi の格子点）
//    ※ 実際の格子表現（THnD/TH2D配列など）は実装側で固定する。
//  - メタ情報（窓・分解能・N_theta・P_mu・ビニング等）も併せて保存する。
//
// 注意
//  - soft photon 発散があるため、呼び出し側（解析窓など）で Eg > Eg_min のカットが必須。
//  - 格子ビニング、生成統計（サンプル数、シード、スメア回数など）は
//    このモジュール内部に固定で持つ（ヘッダには出さない）。
// ============================================================

struct AnalysisWindow4D;

// 成功: 0、失敗: 非0
int MakeRMDGridPdf(const char* out_filepath, const char* key);

// 真値生成窓（Ee, Eg）を解析窓から自動的に作って格子PDFを作成する。
//  - truth_win の Ee/Eg を「解析窓」と見做し、energy_response_shape_e/g の 0.1 倍点まで広げた窓で真値をサンプルする
//  - 非物理領域（負のエネルギーや運動学上限超え）は除外する
//  - t, theta は未使用（観測窓は AnalysisWindow.h の analysis_window を使う）
//  - 単位は MeV（Ee, Eg）
int MakeRMDGridPdfWithTruthWindow(const char* out_filepath,
                                  const char* key,
                                  const AnalysisWindow4D& analysis_win_for_truth);

#endif // P2MEG_MAKE_RMD_GRID_PDF_H
