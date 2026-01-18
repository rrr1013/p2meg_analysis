#ifndef P2MEG_MATH_UTILS_H
#define P2MEG_MATH_UTILS_H

#include <cmath>

#include "p2meg/DetectorResolution.h"

// ============================================================
// p2MEG 数値ユーティリティ
//
// 物理カットではない数値ガードをまとめる。
// ============================================================

// 有限値チェック
static inline bool Math_IsFinite(double x) {
  return std::isfinite(x);
}

// 範囲クリップ（物理カットではない）
static inline double Math_Clamp(double x, double lo, double hi) {
  if (x < lo) return lo;
  if (x > hi) return hi;
  return x;
}

// detres.N_theta を int として安全に使う
// - N_theta >= 1 を保証する（物理カットではない）
static inline int Math_GetNTheta(const DetectorResolutionConst& res) {
  const double x = static_cast<double>(res.N_theta);
  int N = static_cast<int>(std::lround(x));
  if (N < 1) N = 1;
  return N;
}

#endif // P2MEG_MATH_UTILS_H
