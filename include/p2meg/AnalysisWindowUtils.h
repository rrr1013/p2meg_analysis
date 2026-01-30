// include/p2meg/AnalysisWindowUtils.h

#ifndef P2MEG_ANALYSIS_WINDOW_UTILS_H
#define P2MEG_ANALYSIS_WINDOW_UTILS_H

#include <cmath>
#include "p2meg/AnalysisWindow.h"

// ============================================================
// p2MEG 解析窓ユーティリティ
//
// 単位:
//  - Ee, Eg: MeV
//  - t     : ns
//  - theta : rad
//
// TSB（タイミングサイドバンド）:
//  - 全時間範囲 [t_all_min, t_all_max] を与え、
//    解析窓 [win.t_min, win.t_max] の外側を TSB とみなす。
//  - Ee/Eg/theta の窓内条件と組み合わせた判定も用意する。
//
// ESB（エネルギーサイドバンド, Eg）:
//  - 全エネルギー範囲 [Eg_all_min, Eg_all_max] を与え、
//    解析窓 [win.Eg_min, win.Eg_max] の外側を ESB とみなす。
//  - Ee/theta の窓内条件と組み合わせた判定も用意する。
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

// ---- TSB 判定：全時間範囲内で、解析窓の外側か？ ----
static inline bool AnalysisWindow_InTimeSideband(const AnalysisWindow4D& win,
                                                 double t,
                                                 double t_all_min,
                                                 double t_all_max) {
  if (!std::isfinite(t) || !std::isfinite(t_all_min) || !std::isfinite(t_all_max)) return false;
  if (t_all_max <= t_all_min) return false;

  if (t < t_all_min || t > t_all_max) return false;

  // 端の重複を避けるため、左右で不等号を非対称にする
  const bool in_left  = (t >= t_all_min && t <  win.t_min);
  const bool in_right = (t >  win.t_max  && t <= t_all_max);
  return in_left || in_right;
}

// ---- TSB 判定： (Ee, Eg, theta) は解析窓内、t は全時間範囲内で解析窓外 ----
static inline bool AnalysisWindow_InTSB(const AnalysisWindow4D& win,
                                       double Ee, double Eg, double theta,
                                       double t,
                                       double t_all_min,
                                       double t_all_max) {
  if (!AnalysisWindow_In3D(win, Ee, Eg, theta)) return false;
  return AnalysisWindow_InTimeSideband(win, t, t_all_min, t_all_max);
}

// ---- TSB 幅（スケール係数用）。不正時は 0 ----
static inline double AnalysisWindow_TimeSidebandWidth(const AnalysisWindow4D& win,
                                                      double t_all_min,
                                                      double t_all_max) {
  if (!std::isfinite(t_all_min) || !std::isfinite(t_all_max)) return 0.0;
  if (t_all_max <= t_all_min) return 0.0;

  const double left  = (win.t_min > t_all_min) ? (win.t_min - t_all_min) : 0.0;
  const double right = (t_all_max > win.t_max) ? (t_all_max - win.t_max) : 0.0;
  const double w = left + right;
  return (w > 0.0 && std::isfinite(w)) ? w : 0.0;
}

// ============================================================
// ESB (Eg) utilities
// ============================================================

// ---- ESB(Eg) 判定：全 Eg 範囲内で、解析窓(Eg) の外側か？ ----
static inline bool AnalysisWindow_InEnergySidebandEg(const AnalysisWindow4D& win,
                                                     double Eg,
                                                     double Eg_all_min,
                                                     double Eg_all_max) {
  if (!std::isfinite(Eg) || !std::isfinite(Eg_all_min) || !std::isfinite(Eg_all_max)) return false;
  if (Eg_all_max <= Eg_all_min) return false;

  if (Eg < Eg_all_min || Eg > Eg_all_max) return false;

  // 端の重複を避けるため、左右で不等号を非対称にする
  const bool in_left  = (Eg >= Eg_all_min && Eg <  win.Eg_min);
  const bool in_right = (Eg >  win.Eg_max  && Eg <= Eg_all_max);
  return in_left || in_right;
}

// ---- ESB(Eg) 判定： (Ee, theta) は解析窓内、Eg は ESB、t は全時間範囲内 ----
// 注意: ESB での t フィット用途を想定し、t は解析窓内に縛らない。
static inline bool AnalysisWindow_InESB_Eg(const AnalysisWindow4D& win,
                                          double Ee, double Eg, double theta,
                                          double t,
                                          double Eg_all_min,
                                          double Eg_all_max,
                                          double t_all_min,
                                          double t_all_max) {
  if (!std::isfinite(t) || !std::isfinite(t_all_min) || !std::isfinite(t_all_max)) return false;
  if (t_all_max <= t_all_min) return false;
  if (t < t_all_min || t > t_all_max) return false;

  if (Ee < win.Ee_min || Ee > win.Ee_max) return false;
  if (theta < win.theta_min || theta > win.theta_max) return false;

  return AnalysisWindow_InEnergySidebandEg(win, Eg, Eg_all_min, Eg_all_max);
}

// ---- ESB(Eg) 幅（メタ情報・sanity check 用）。不正時は 0 ----
static inline double AnalysisWindow_EnergySidebandWidthEg(const AnalysisWindow4D& win,
                                                          double Eg_all_min,
                                                          double Eg_all_max) {
  if (!std::isfinite(Eg_all_min) || !std::isfinite(Eg_all_max)) return 0.0;
  if (Eg_all_max <= Eg_all_min) return 0.0;

  const double left  = (win.Eg_min > Eg_all_min) ? (win.Eg_min - Eg_all_min) : 0.0;
  const double right = (Eg_all_max > win.Eg_max) ? (Eg_all_max - win.Eg_max) : 0.0;
  const double w = left + right;
  return (w > 0.0 && std::isfinite(w)) ? w : 0.0;
}

#endif // P2MEG_ANALYSIS_WINDOW_UTILS_H
