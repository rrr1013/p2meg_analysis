// include/p2meg/RMDSpectrum.h
#ifndef P2MEG_RMD_SPECTRUM_H
#define P2MEG_RMD_SPECTRUM_H

// ============================================================
// Radiative Muon Decay (RMD) : μ+ → e+ ν ν̄ γ
//
// ここでは「停止ミューオン静止系」での tree-level SM(V-A) の式を使い、
// 分岐比（branching ratio）に比例する微分分布（shape）を返す。
// 単位系は x,y が無次元、Ee,Eg は MeV 入力。
// ============================================================

//------------------------------------------------------------
// 「全体回転」を平均した後の三重微分（偏極を落とした形）
//   d^3B/(dx dy dcosθ_eγ)
//------------------------------------------------------------
double RMD_d3B_dxdy_dcos(double x, double y, double cosThetaEG, double d_min = 1e-6);

double RMD_d3B_dEe_dEg_dcos(double Ee_MeV, double Eg_MeV, double cosThetaEG, double d_min = 1e-6);

//------------------------------------------------------------
// 完全式（F,G,H を含む・偏極あり）
//   d^6B/(dx dy dΩ_e dΩ_γ)
//
// 入力として
//   cosThetaEG = p̂_e · p̂_γ
//   cosThetaE  = P̂ · p̂_e
//   cosThetaG  = P̂ · p̂_γ
// を受け取り、Pmu は偏極度（スカラー；符号込み）として扱う。
// （P ベクトルが必要なら、Pmu=|P| とし cosThetaE,G を P̂ との内積で与える）
//------------------------------------------------------------
double RMD_d6B_dxdy_dOmegae_dOmegag(double x, double y,
                                   double cosThetaEG,
                                   double cosThetaE,
                                   double cosThetaG,
                                   double Pmu,
                                   double d_min = 1e-6);

double RMD_d6B_dEe_dEg_dOmegae_dOmegag(double Ee_MeV, double Eg_MeV,
                                      double cosThetaEG,
                                      double cosThetaE,
                                      double cosThetaG,
                                      double Pmu,
                                      double d_min = 1e-6);

#endif // P2MEG_RMD_SPECTRUM_H
