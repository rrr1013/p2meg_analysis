#ifndef P2MEG_MICHEL_SPECTRUM_H
#define P2MEG_MICHEL_SPECTRUM_H

#include "p2meg/Constants.h"

// ============================================================
// ミシェル崩壊スペクトル（shape = 形状）ユーティリティ
//
// ・返り値は「非正規化の shape」（微分崩壊率に比例）
// ・全体の規格化定数は尤度側（イベント数など）に吸収する想定
//
// 実装上の仮定（このファイルの関数が使うモデル）
// ・μ静止系
// ・端点エネルギーは Emax = m_mu/2 を採用（me による端点補正は入れていない）
// ・me は η 項の x0 = me/Emax にのみ反映（簡略化）
// ============================================================

// ---- 補助関数 ----
double Michel_Emax(const ParticleMasses& ms);          // Emax = m_mu/2
double Michel_x(double Ee, const ParticleMasses& ms);  // x = Ee/Emax
double Michel_x0(const ParticleMasses& ms);            // x0 = me/Emax

// ---- スペクトル（非正規化 shape） ----
// Ee [MeV], costh = cos(theta), P_mu は [-1,1] を想定
double Michel_d2Shape_dE_dCosTheta(double Ee,
                                   double costh,
                                   double P_mu,
                                   const MichelParams& mp = kMichelSM,
                                   const ParticleMasses& ms = kMassesPDG);

double Michel_dShape_dE(double Ee,
                        const MichelParams& mp = kMichelSM,
                        const ParticleMasses& ms = kMassesPDG);

#endif // P2MEG_MICHEL_SPECTRUM_H
