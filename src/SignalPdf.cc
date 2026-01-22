// src/SignalPdf.cc
#include "p2meg/SignalPdf.h"

#include <cmath>   // exp, sqrt, erf, isfinite, llround, fmod, fabs

#include "p2meg/AngleUtils.h"
#include "p2meg/MathUtils.h"

// ------------------------------------------------------------
// 内部: 標準正規PDF/CDFとトランケート正規
// ------------------------------------------------------------
static double NormalPdf(double x, double mu, double sigma) {
  if (!(sigma > 0.0)) return 0.0;
  if (!std::isfinite(x) || !std::isfinite(mu) || !std::isfinite(sigma)) return 0.0;

  const double z = (x - mu) / sigma;
  const double norm = 1.0 / (std::sqrt(2.0 * pi) * sigma);
  return norm * std::exp(-0.5 * z * z);
}

static double NormalCdf(double x, double mu, double sigma) {
  if (!(sigma > 0.0)) return 0.0;
  if (!std::isfinite(x) || !std::isfinite(mu) || !std::isfinite(sigma)) return 0.0;

  const double z = (x - mu) / (sigma * std::sqrt(2.0));
  // Φ(x) = 1/2 [1 + erf(z)]
  return 0.5 * (1.0 + std::erf(z));
}

static double TruncNormalPdf(double x, double mu, double sigma,
                             double a, double b) {
  // a<=x<=b の窓内で正規化したトランケート正規
  if (!(sigma > 0.0)) return 0.0;
  if (!(a < b)) return 0.0;
  if (!std::isfinite(x) || !std::isfinite(mu) || !std::isfinite(sigma)) return 0.0;
  if (!std::isfinite(a) || !std::isfinite(b)) return 0.0;

  if (x < a || x > b) return 0.0;

  const double Z = NormalCdf(b, mu, sigma) - NormalCdf(a, mu, sigma);
  if (!(Z > 0.0) || !std::isfinite(Z)) return 0.0;

  const double p = NormalPdf(x, mu, sigma) / Z;
  return std::isfinite(p) ? p : 0.0;
}

// ------------------------------------------------------------
// 内部: 角度の離散格子への丸め
//   theta_i = i*pi/N_theta (i=0..N_theta)
// 入力 theta を最も近い格子点に丸め、インデックス i を返す。
// 不正なら -1 を返す。
// ------------------------------------------------------------
static int ThetaNearestIndex(double theta, int N_theta) {
  if (!(N_theta >= 1)) return -1;
  if (!std::isfinite(theta)) return -1;
  if (theta < 0.0 || theta > pi) return -1;

  const double step = pi / static_cast<double>(N_theta);
  if (!(step > 0.0) || !std::isfinite(step)) return -1;

  long long i_ll = std::llround(theta / step);

  // 数値丸めの端のはみ出しを安全にクリップ（物理カットではない）
  if (i_ll < 0LL) i_ll = 0LL;
  if (i_ll > static_cast<long long>(N_theta)) i_ll = static_cast<long long>(N_theta);

  return static_cast<int>(i_ll);
}

// ------------------------------------------------------------
// 内部: theta=pi に丸め込まれる領域の面積を数値的に評価
//  - phi_e, phi_g の離散格子と許可マスクに従う
// ------------------------------------------------------------
static double ComputeAreaPi(const DetectorResolutionConst& res) {
  const int N_phi_e = Math_GetNPhiE(res);
  const int N_phi_g = Math_GetNPhiG(res);
  if (!(N_phi_e >= 1) || !(N_phi_g >= 1)) return 0.0;

  if (!Detector_IsPhiRangeValid(res.phi_e_min, res.phi_e_max, N_phi_e)) return 0.0;
  if (!Detector_IsPhiRangeValid(res.phi_g_min, res.phi_g_max, N_phi_g)) return 0.0;

  double area = 0.0;
  for (int ie = 0; ie <= N_phi_e; ++ie) {
    const double phi_e = Detector_PhiGridPoint(ie, res.phi_e_min, res.phi_e_max, N_phi_e);
    const double w_e = Detector_PhiBinWidth(ie, res.phi_e_min, res.phi_e_max, N_phi_e);
    if (!(w_e > 0.0)) continue;
    for (int ig = 0; ig <= N_phi_g; ++ig) {
      if (!Detector_IsAllowedPhiPairIndex(ie, ig, res)) continue;
      const double phi_g = Detector_PhiGridPoint(ig, res.phi_g_min, res.phi_g_max, N_phi_g);
      const double w_g = Detector_PhiBinWidth(ig, res.phi_g_min, res.phi_g_max, N_phi_g);
      if (!(w_g > 0.0)) continue;

      const double theta_eg = std::fabs(phi_e - phi_g);
      const int ith = ThetaNearestIndex(theta_eg, res.N_theta);
      if (ith == res.N_theta) area += w_e * w_g;
    }
  }
  return (area > 0.0 && std::isfinite(area)) ? area : 0.0;
}

// ------------------------------------------------------------
// 公開: SignalPdf
// ------------------------------------------------------------
double SignalPdf(double Ee, double Eg, double t,
                 double phi_detector_e, double phi_detector_g,
                 const AnalysisWindow4D& win,
                 const DetectorResolutionConst& res,
                 const ParticleMasses& ms) {
  // 基本チェック
  if (!std::isfinite(Ee) || !std::isfinite(Eg) || !std::isfinite(t)) return 0.0;

  // 分解能チェック（角度は離散なので sigma_theta は不要）
  if (!(res.sigma_t  > 0.0)) return 0.0;
  if (!(res.N_theta  >= 1))  return 0.0;

  // phi の許可領域チェック（範囲 + マスク）
  int idx_e = -1;
  int idx_g = -1;
  if (!Detector_IsAllowedPhiPairValue(phi_detector_e, phi_detector_g, res, idx_e, idx_g)) {
    return 0.0;
  }

  const int N_phi_e = Math_GetNPhiE(res);
  const int N_phi_g = Math_GetNPhiG(res);
  const double phi_e_snap = Detector_PhiGridPoint(idx_e, res.phi_e_min, res.phi_e_max, N_phi_e);
  const double phi_g_snap = Detector_PhiGridPoint(idx_g, res.phi_g_min, res.phi_g_max, N_phi_g);

  // phi から相対角 theta_eg を作る（phi ベースで角度評価）
  const double theta_eg = std::fabs(phi_e_snap - phi_g_snap);
  if (!(theta_eg >= 0.0)) return 0.0;

  // theta を最も近い格子点に丸めてから以後の判定・評価に使う
  const int ith = ThetaNearestIndex(theta_eg, res.N_theta);
  if (ith < 0) return 0.0;

  const double step = pi / static_cast<double>(res.N_theta);
  const double theta_snap = step * static_cast<double>(ith);

  // 解析窓チェック（窓外は0）※角度は丸め後 theta_snap で判定
  if (Ee < win.Ee_min || Ee > win.Ee_max) return 0.0;
  if (Eg < win.Eg_min || Eg > win.Eg_max) return 0.0;
  if (t  < win.t_min  || t  > win.t_max ) return 0.0;
  if (theta_snap < win.theta_min || theta_snap > win.theta_max) return 0.0;

  // 真値（停止μの2体崩壊）
  const double Ee0 = 0.5 * ms.m_mu;
  const double Eg0 = 0.5 * ms.m_mu;
  const double t0  = res.t_mean;

  // Ee, Eg は energy_response_shape_e/g で評価（窓内で正規化）
  const double pEe = energy_response_pdf_window_e(Ee, Ee0, win.Ee_min, win.Ee_max);
  if (!(pEe > 0.0)) return 0.0;

  const double pEg = energy_response_pdf_window_g(Eg, Eg0, win.Eg_min, win.Eg_max);
  if (!(pEg > 0.0)) return 0.0;

  const double pt  = TruncNormalPdf(t,  t0,  res.sigma_t,  win.t_min,  win.t_max);
  if (!(pt > 0.0)) return 0.0;

  // 角度: 離散（散乱なし）→ 信号の角度デルタ条件として theta=pi のみに重み
  // theta_i = i*pi/N_theta なので、pi は i=N_theta に対応
  // phi 空間での正規化を満たすよう、theta=pi に丸め込まれる領域面積を使う。
  const double area_pi = ComputeAreaPi(res);
  if (!(area_pi > 0.0) || !std::isfinite(area_pi)) return 0.0;
  const double pth = (ith == res.N_theta) ? (1.0 / area_pi) : 0.0;
  if (!(pth > 0.0)) return 0.0;

  const double p = pEe * pEg * pt * pth;
  return std::isfinite(p) ? p : 0.0;
}
