#ifndef P2MEG_HIST_UTILS_H
#define P2MEG_HIST_UTILS_H

#include "THn.h"
#include "TAxis.h"

// ============================================================
// p2MEG ヒストグラム補助関数
// ============================================================

// 軸が等間隔ビンである前提で、座標 x を隣り合う2ビンに落とす（Ee,Eg 用）。
// 値は「ビン中心に定義された格子値」と見做して補間するための (i0, i1, f) を返す。
//  - i0, i1 は [1..nbins] のビン番号
//  - f は 0..1 の補間係数（i0 側が (1-f)、i1 側が f）
int Hist_AxisBracketUniform(const TAxis& ax, double x, int& i0, int& i1, double& f);

// Ee,Eg を 2D（4点）補間して、phi軸は固定ビン（最近傍）で評価する。
double Hist_InterpEeEg4(const THnD& h, double Ee, double Eg, int bin_phi_e, int bin_phi_g);

// THnD（4D）の全ビン総和
double Hist_SumAllBins4(const THnD& h);

#endif // P2MEG_HIST_UTILS_H
