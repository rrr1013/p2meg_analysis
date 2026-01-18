#ifndef P2MEG_ANGLE_UTILS_H
#define P2MEG_ANGLE_UTILS_H

#include <cmath>

#include "p2meg/Constants.h"
#include "p2meg/MathUtils.h"

// ============================================================
// p2MEG 角度ユーティリティ
//
// 単位: rad
// ============================================================

// phi を [0,pi] にクリップする（不正入力は 0 に落とす）
static inline double Angle_ClipPhi0Pi(double phi) {
  if (!Math_IsFinite(phi)) return 0.0;
  if (phi < 0.0) return 0.0;
  if (phi > pi) return pi;
  return phi;
}

// phi を [0,pi] にクリップして N_theta の離散点に丸める
//  - theta_i = i*pi/N_theta (i=0..N_theta)
static inline double Angle_DiscretizePhi(double phi, int N_theta) {
  if (!(N_theta >= 1)) return 0.0;

  const double p = Angle_ClipPhi0Pi(phi);
  const double step = pi / static_cast<double>(N_theta);
  if (!(step > 0.0) || !Math_IsFinite(step)) return 0.0;

  long long i = std::llround(p / step);
  if (i < 0LL) i = 0LL;
  if (i > static_cast<long long>(N_theta)) i = static_cast<long long>(N_theta);
  return step * static_cast<double>(i);
}

// phi から e-γ の相対角 theta を作る（範囲外は不正扱い）
//  - phi は [0,pi] を想定
//  - 不正なら -1 を返す（物理カットではない）
static inline double Angle_ThetaFromPhiStrict(double phi_e, double phi_g) {
  if (!Math_IsFinite(phi_e) || !Math_IsFinite(phi_g)) return -1.0;
  if (phi_e < 0.0 || phi_e > pi) return -1.0;
  if (phi_g < 0.0 || phi_g > pi) return -1.0;
  return std::fabs(phi_e - phi_g);
}

// phi から e-γ の相対角 theta を作る（クリップ付き）
//  - 数値ガードとして [0,pi] に収める
static inline double Angle_ThetaFromPhiClipped(double phi_e, double phi_g) {
  if (!Math_IsFinite(phi_e) || !Math_IsFinite(phi_g)) return 0.0;
  const double pe = Angle_ClipPhi0Pi(phi_e);
  const double pg = Angle_ClipPhi0Pi(phi_g);
  return std::fabs(pe - pg);
}

#endif // P2MEG_ANGLE_UTILS_H
