#ifndef P2MEG_ACC_TIME_FIT_H
#define P2MEG_ACC_TIME_FIT_H

#include <algorithm>
#include <cmath>

#include "TF1.h"
#include "TFitResultPtr.h"
#include "TH1D.h"

// ============================================================
// ACC 時間形状 fit の共通ヘルパ
//
// 目的:
//  - ACC 数見積もり用 macro と ACC PDF 生成器で、
//    同じ時間形状 model / 同じ fit 仕様を使い回す。
//
// model:
//   f(t) = A * exp(-t^2 / (2 sigma^2)) + C
//
// 注意:
//  - t_reject_min <= t <= t_reject_max は fit から除外する。
//  - これは AW を見ずに TSB だけから時間形状を学習するための設定である。
// ============================================================

struct AccTimeFitResult {
  int fit_status;   // ROOT fit status
  double chi2;      // chi2
  int ndf;          // number of degrees of freedom
  double chi2_ndf;  // chi2/ndf
  double A;         // [counts/ns]
  double sigma;     // [ns]
  double C;         // [counts/ns]
};

static double gAccTimeFitRejectMin = 0.0;
static double gAccTimeFitRejectMax = 0.0;

// 共通の密度 model
static inline double AccTimeFit_DensityCoreValue(double t, double A, double sigma, double C) {
  if (!(sigma > 0.0) || !std::isfinite(A) || !std::isfinite(sigma) || !std::isfinite(C)) {
    return 0.0;
  }
  const double g = std::exp(-0.5 * (t * t) / (sigma * sigma));
  const double v = A * g + C;
  return (v > 0.0 && std::isfinite(v)) ? v : 0.0;
}

static inline double AccTimeFit_DensityCore(double* x, double* p) {
  return AccTimeFit_DensityCoreValue(x[0], p[0], p[1], p[2]);
}

// reject window を持つ fit 用
static inline double AccTimeFit_DensityRejectWindow(double* x, double* p) {
  const double t = x[0];
  if (t >= gAccTimeFitRejectMin && t <= gAccTimeFitRejectMax) {
    TF1::RejectPoint();
    return 0.0;
  }
  return AccTimeFit_DensityCore(x, p);
}

// TH1D を同じ仕様で fit する
static inline int AccTimeFit_Run(TH1D& h,
                                 double t_reject_min,
                                 double t_reject_max,
                                 double sigma_min,
                                 double sigma_max,
                                 AccTimeFitResult& out) {
  out.fit_status = -999;
  out.chi2 = 0.0;
  out.ndf = 0;
  out.chi2_ndf = 0.0;
  out.A = 0.0;
  out.sigma = 0.0;
  out.C = 0.0;

  if (h.GetEntries() <= 0.0) return 1;

  const double pedestal0 =
      0.5 * (h.GetBinContent(1) + h.GetBinContent(h.GetNbinsX()));
  const double amp0 = std::max(0.0, h.GetMaximum() - pedestal0);

  gAccTimeFitRejectMin = t_reject_min;
  gAccTimeFitRejectMax = t_reject_max;

  TF1 f("acc_time_fit_shared", AccTimeFit_DensityRejectWindow,
        h.GetXaxis()->GetXmin(), h.GetXaxis()->GetXmax(), 3);
  f.SetParNames("A", "sigma", "C");
  f.SetParameters(amp0, 60.0, std::max(0.0, pedestal0));
  f.SetParLimits(0, 0.0, std::max(10.0 * h.GetMaximum(), 1.0));
  f.SetParLimits(1, sigma_min, sigma_max);
  f.SetParLimits(2, 0.0, std::max(5.0 * h.GetMaximum(), 1.0));

  TFitResultPtr fit_result = h.Fit(&f, "SRQ0");
  out.fit_status = static_cast<int>(fit_result);
  out.chi2 = f.GetChisquare();
  out.ndf = f.GetNDF();
  out.chi2_ndf = (out.ndf > 0) ? (out.chi2 / static_cast<double>(out.ndf)) : 0.0;
  out.A = f.GetParameter(0);
  out.sigma = f.GetParameter(1);
  out.C = f.GetParameter(2);

  if (!(out.sigma > 0.0) || !std::isfinite(out.A) ||
      !std::isfinite(out.sigma) || !std::isfinite(out.C)) {
    return 2;
  }
  return 0;
}

// fit 結果から density ヒストグラムを作る
static inline int AccTimeFit_FillDensityHistogramFromResult(TH1D& h_out,
                                                            const AccTimeFitResult& fit) {
  if (!(fit.sigma > 0.0) || !std::isfinite(fit.A) ||
      !std::isfinite(fit.sigma) || !std::isfinite(fit.C)) {
    return 1;
  }
  for (int ib = 1; ib <= h_out.GetNbinsX(); ++ib) {
    const double lo = h_out.GetXaxis()->GetBinLowEdge(ib);
    const double hi = lo + h_out.GetXaxis()->GetBinWidth(ib);
    const double w = hi - lo;
    double dens = 0.0;
    if (w > 0.0 && std::isfinite(w)) {
      TF1 f_eval("acc_time_eval_shared", AccTimeFit_DensityCore,
                 h_out.GetXaxis()->GetXmin(), h_out.GetXaxis()->GetXmax(), 3);
      f_eval.SetParameters(fit.A, fit.sigma, fit.C);
      dens = f_eval.Integral(lo, hi) / w;
    }
    h_out.SetBinContent(ib, (dens > 0.0 && std::isfinite(dens)) ? dens : 0.0);
    h_out.SetBinError(ib, 0.0);
  }
  return 0;
}

#endif // P2MEG_ACC_TIME_FIT_H
