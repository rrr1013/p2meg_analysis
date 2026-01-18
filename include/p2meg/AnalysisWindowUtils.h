#ifndef P2MEG_ANALYSIS_WINDOW_UTILS_H
#define P2MEG_ANALYSIS_WINDOW_UTILS_H

#include "p2meg/AnalysisWindow.h"

// ============================================================
// p2MEG 解析窓ユーティリティ
//
// 単位:
//  - Ee, Eg: MeV
//  - t     : ns
//  - theta : rad
// ============================================================

// 4D 解析窓判定（窓外は false）
static inline bool AnalysisWindow_In4D(const AnalysisWindow4D& win,
                                       double Ee, double Eg, double t, double theta) {
  if (Ee < win.Ee_min || Ee > win.Ee_max) return false;
  if (Eg < win.Eg_min || Eg > win.Eg_max) return false;
  if (t  < win.t_min  || t  > win.t_max ) return false;
  if (theta < win.theta_min || theta > win.theta_max) return false;
  return true;
}

// 3D 解析窓判定（Ee, Eg, theta）
static inline bool AnalysisWindow_In3D(const AnalysisWindow4D& win,
                                       double Ee, double Eg, double theta) {
  if (Ee < win.Ee_min || Ee > win.Ee_max) return false;
  if (Eg < win.Eg_min || Eg > win.Eg_max) return false;
  if (theta < win.theta_min || theta > win.theta_max) return false;
  return true;
}

// 1D 解析窓判定（t のみ）
static inline bool AnalysisWindow_InTime(const AnalysisWindow4D& win, double t) {
  if (t < win.t_min || t > win.t_max) return false;
  return true;
}

#endif // P2MEG_ANALYSIS_WINDOW_UTILS_H
