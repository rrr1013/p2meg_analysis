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
//  - 時間 t は解析窓内で一様として評価側で解析的に掛ける。
//
// 設計方針（重要）
//  - 入力イベントは (Ee, Eg, t, phi_detector_e, phi_detector_g)。
//  - TSB（タイミングサイドバンド）は
//      「t が全時間範囲内で、解析窓外」にあること。
//    Ee/Eg/角度は解析窓内にあるものだけを採用する。
//  - phi は [0,pi] にクリップしてから N_theta の離散点に丸める。
//    phi_i = i * pi / N_theta (i=0..N_theta)
//  - theta_eg は離散化後の phi から |phi_e - phi_g| として作る。
//  - 正規化は
//      Σ_{phi_e,phi_g} ∫ dEe dEg p4(Ee,Eg,phi_e,phi_g) = 1
//    となるように行う（Ee/Eg の bin 幅のみを測度に入れる）。
//
// 出力（ROOT）
//  - key 名の 4D THnD（Ee, Eg, phi_e, phi_g）
//  - メタ情報（N_theta、ビニング、phi 軸定義、正規化条件）
// ============================================================

// 成功: 0、失敗: 非0
int MakeACCGridPdf(const std::vector<Event>& events,
                   const char* out_filepath,
                   const char* key);

#endif // P2MEG_MAKE_ACC_GRID_PDF_H
