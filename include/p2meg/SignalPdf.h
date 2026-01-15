// include/p2meg/SignalPdf.h
#ifndef P2MEG_SIGNAL_PDF_H
#define P2MEG_SIGNAL_PDF_H

// ============================================================
// p2MEG Signal PDF (μ+→e+γ)
//
// 単位:
//  - Ee, Eg         : MeV
//  - t              : ns
//  - phi_detector_* : rad
//
// モデル:
//  - 真値: Ee0=Eg0=m_mu/2, theta0=pi, t0=t_mean
//  - Ee, Eg, t は独立なトランケート正規分布（解析窓内で正規化）
//  - 角度は離散化して扱う（角度スメアなし）
//      theta_i = i * pi / N_theta   (i = 0..N_theta)
//    入力は (phi_detector_e, phi_detector_g) とし、phi ベースでの角度評価を行う。
//      theta_eg = |phi_detector_e - phi_detector_g|
//    を [0,pi] で扱い、最も近い格子点に丸めて評価する。
//    信号は理想化（散乱なし）により theta=pi のみに重みを持つ。
//    したがって角度因子は
//      P(theta=pi)=1/Area_pi, それ以外=0
//    とし、phi 空間での正規化（∫dphi_e dphi_g = 1）を満たす。
//    Δtheta = pi/N_theta は角度離散化ビン幅で、
//    theta=pi に丸め込まれる領域の面積は Area_pi = (Δtheta/2)^2 となる。
//    （丸め後の theta_i が解析窓外なら 0）
//
// 戻り値:
//  - 解析窓外、phi が不正、thetaが[0,pi]外、sigma<=0、N_theta<1、正規化定数<=0、非数などは 0
// ============================================================

#include "p2meg/AnalysisWindow.h"
#include "p2meg/DetectorResolution.h"
#include "p2meg/Constants.h"

double SignalPdf(double Ee, double Eg, double t,
                 double phi_detector_e, double phi_detector_g,
                 const AnalysisWindow4D& win,
                 const DetectorResolutionConst& res,
                 const ParticleMasses& ms = kMassesPDG);

#endif // P2MEG_SIGNAL_PDF_H
