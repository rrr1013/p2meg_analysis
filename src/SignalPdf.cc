#include "p2meg/SignalPdf.h"

#include <cmath>   // exp, sqrt, erf, isfinite

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
// 内部: half-normal（折り返し）とそのトランケート版
// ------------------------------------------------------------
static double HalfNormalPdf(double d, double sigma) {
  // d >= 0 で half-normal: sqrt(2/pi)/sigma * exp(-d^2/(2 sigma^2))
  if (!(sigma > 0.0)) return 0.0;
  if (!std::isfinite(d) || !std::isfinite(sigma)) return 0.0;

  if (d < 0.0) return 0.0;

  const double z = d / sigma;
  const double norm = std::sqrt(2.0 / pi) / sigma;
  return norm * std::exp(-0.5 * z * z);
}

static double TruncHalfNormalPdf(double d, double sigma,
                                 double dmin, double dmax) {
  // dmin<=d<=dmax（かつ d>=0）で窓内正規化した half-normal
  if (!(sigma > 0.0)) return 0.0;
  if (!(dmin < dmax)) return 0.0;
  if (!std::isfinite(d) || !std::isfinite(sigma)) return 0.0;
  if (!std::isfinite(dmin) || !std::isfinite(dmax)) return 0.0;

  // δ>=0 の定義なので、窓端も 0 未満なら物理的に扱わない（0を返す）
  if (dmax < 0.0) return 0.0;

  // half-normal は d>=0 のみなので下限を 0 にクリップ
  const double a = (dmin < 0.0) ? 0.0 : dmin;
  const double b = dmax;

  if (!(a < b)) return 0.0;
  if (d < a || d > b) return 0.0;

  // half-normal のCDF: F(d)=erf(d/(sqrt2*sigma))  (d>=0)
  const double u = a / (std::sqrt(2.0) * sigma);
  const double v = b / (std::sqrt(2.0) * sigma);

  const double Z = std::erf(v) - std::erf(u);
  if (!(Z > 0.0) || !std::isfinite(Z)) return 0.0;

  const double p = HalfNormalPdf(d, sigma) / Z;
  return std::isfinite(p) ? p : 0.0;
}

// ------------------------------------------------------------
// 公開: SignalPdf
// ------------------------------------------------------------
double SignalPdf(double Ee, double Eg, double t, double theta,
                 const AnalysisWindow4D& win,
                 const DetectorResolutionConst& res,
                 const ParticleMasses& ms) {
  // 基本チェック
  if (!std::isfinite(Ee) || !std::isfinite(Eg) || !std::isfinite(t) || !std::isfinite(theta)) return 0.0;

  // theta は [0,pi] を想定（解析窓以前の物理的チェック）
  if (theta < 0.0 || theta > pi) return 0.0;

  // 解析窓チェック（窓外は0）
  if (Ee < win.Ee_min || Ee > win.Ee_max) return 0.0;
  if (Eg < win.Eg_min || Eg > win.Eg_max) return 0.0;
  if (t  < win.t_min  || t  > win.t_max ) return 0.0;
  if (theta < win.theta_min || theta > win.theta_max) return 0.0;

  // 分解能チェック
  if (!(res.sigma_Ee > 0.0)) return 0.0;
  if (!(res.sigma_Eg > 0.0)) return 0.0;
  if (!(res.sigma_t  > 0.0)) return 0.0;
  if (!(res.sigma_theta > 0.0)) return 0.0;

  // 真値（停止μの2体崩壊）
  const double Ee0 = 0.5 * ms.m_mu;
  const double Eg0 = 0.5 * ms.m_mu;
  const double t0  = res.t_mean;

  // Ee, Eg, t の窓内正規化トランケート正規
  const double pEe = TruncNormalPdf(Ee, Ee0, res.sigma_Ee, win.Ee_min, win.Ee_max);
  if (!(pEe > 0.0)) return 0.0;

  const double pEg = TruncNormalPdf(Eg, Eg0, res.sigma_Eg, win.Eg_min, win.Eg_max);
  if (!(pEg > 0.0)) return 0.0;

  const double pt  = TruncNormalPdf(t,  t0,  res.sigma_t,  win.t_min,  win.t_max);
  if (!(pt > 0.0)) return 0.0;

  // 角度: δ = pi - theta (δ>=0) で折り返し half-normal
  // 解析窓 [theta_min, theta_max] は δ 範囲 [pi-theta_max, pi-theta_min]
  const double delta  = pi - theta;
  const double dmin   = pi - win.theta_max;
  const double dmax   = pi - win.theta_min;

  const double pth = TruncHalfNormalPdf(delta, res.sigma_theta, dmin, dmax);
  if (!(pth > 0.0)) return 0.0;

  const double p = pEe * pEg * pt * pth;
  return std::isfinite(p) ? p : 0.0;
}
