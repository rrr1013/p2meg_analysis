#ifndef P2MEG_MAKE_ACC_GRID_PDF_H
#define P2MEG_MAKE_ACC_GRID_PDF_H

#include <vector>

#include "p2meg/Event.h"

// ============================================================
// MakeACCGridPdf
//
// 目的
//  - ACC (accidental) 成分の 4D 格子 PDF
//    p_acc(Ee, Eg, phi_detector_e, phi_detector_g)
//    を ROOT に保存する。
//  - 併せて、最終解析で使う Eg 条件付き時間テンプレート
//      f_t(t | Eg category)
//    も保存する。
//
// 設計方針（重要）
//  - 入力イベントは (Ee, Eg, t, phi_detector_e, phi_detector_g)。
//  - TSB（タイミングサイドバンド）は
//      「t が全時間範囲内で、解析窓外」にあること。
//    Ee/Eg/角度は解析窓内にあるものだけを採用する。
//  - phi は DetectorResolution の範囲にクリップして離散点に丸める。
//    phi_i = phi_min + i * (phi_max-phi_min)/N_phi (i=0..N_phi)
//  - theta_eg は離散化後の phi から |phi_e - phi_g| として作る。
//  - 正規化は
//      Σ_{phi_e,phi_g} ∫ dEe dEg p4(Ee,Eg,phi_e,phi_g) = 1
//    となるように行う（Ee/Eg の bin 幅のみを測度に入れる）。
//
// 出力（ROOT）
//  - key 名の 4D THnD（Ee, Eg, phi_e, phi_g）
//  - key_tshape の 1D TH1D（後方互換用の全 Eg 時間テンプレート密度）
//  - key_tshape_egbin0, key_tshape_egbin1, ... の 1D TH1D
//    （Eg category ごとの t_all 上時間テンプレート密度）
//  - メタ情報（N_phi_e/g、ビニング、phi 軸定義、正規化条件）
// ============================================================

// 成功: 0、失敗: 非0
int MakeACCGridPdf(const std::vector<Event>& events,
                   const char* out_filepath,
                   const char* key);

#endif // P2MEG_MAKE_ACC_GRID_PDF_H
