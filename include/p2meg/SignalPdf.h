#ifndef P2MEG_SIGNAL_PDF_H
#define P2MEG_SIGNAL_PDF_H

// ============================================================
// p2MEG Signal PDF (μ+→e+γ)
//
// 単位:
//  - Ee, Eg: MeV
//  - t     : ns
//  - theta : rad
//
// モデル:
//  - 真値: Ee0=Eg0=m_mu/2, theta0=pi, t0=t_mean
//  - Ee, Eg, t は独立なトランケート正規分布（解析窓内で正規化）
//  - 角度は δ = pi - theta (δ>=0) に写像して half-normal（折り返し）
//    を用い、解析窓に対応する δ 範囲で窓内正規化してから θ に戻す。
//    （|dδ/dθ|=1 のためヤコビアン因子は不要）
//
// 戻り値:
//  - 解析窓外、thetaが[0,pi]外、sigma<=0、正規化定数<=0、非数などは 0
// ============================================================

#include "p2meg/AnalysisWindow.h"
#include "p2meg/DetectorResolution.h"
#include "p2meg/Constants.h"

double SignalPdf(double Ee, double Eg, double t, double theta,
                 const AnalysisWindow4D& win,
                 const DetectorResolutionConst& res,
                 const ParticleMasses& ms = kMassesPDG);

#endif // P2MEG_SIGNAL_PDF_H
