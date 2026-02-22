// check/study_michel_ps_tagged.C
//
// 目的:
//   PS波高ヒストからミシェル判定しきい値を自動決定し、
//   NaIは関数フィット(差の指数)でパルスを分解・積分して、
//   dt の二峰(prompt / delayed)ごとの ModuleA(A1+A2), ModuleB(B1+B2)
//   エネルギーヒストを比較する。
//
// 出力:
//   doc/mainexp/study_michel_ps_tagged_<tag>.pdf
//
// 実行:
//   root -l -q 'check/study_michel_ps_tagged.C'

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>

#include "TCanvas.h"
#include "TF1.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"

// ============================================================
// 入出力
// ============================================================

static const char* kInputDir = "data/rawdata/mainexp/120°/run8";
static const char* kInputTag = "120°_run8";

static const char* kPSFiles[2] = {
  "wave_PS_A_run8.txt",
  "wave_PS_B_run8.txt"
};

static const char* kNaIFiles[4] = {
  "wave_NaI_A1_run8.txt",
  "wave_NaI_A2_run8.txt",
  "wave_NaI_B1_run8.txt",
  "wave_NaI_B2_run8.txt"
};

static const char* kOutputDir = "doc/mainexp";

// ============================================================
// 共通パラメータ
// ============================================================

static const int    kSamplesPerEvent = 1000;
static const double kDtNs = 4.0;

// PSパルス抽出（check/hist_PS_peak.C と同じ簡易条件）
static const double kPSSeedThr = 30.0;
static const double kPSEndThr  = 15.0;
static const int    kPSMinSep  = 100;
static const int    kPSMinWidth = 1;

// NaI候補抽出（hist_NaI_fitint_clean.C と同系）
static const double kNaISeedThr = 60.0;
static const double kNaIEndThr  = 30.0;
static const int    kNaIMinSep  = 100;

// NaIフィット窓・品質
static const int    kFitPreWindow = 40;
static const int    kFitHalfWindow = 150;
static const double kT0FitRange = 30.0;
static const double kRelShapeSys = 0.02;
static const double kChi2Cut = 60.0;
static const int    kMaxAcceptedPerCh = 6;
static const int    kMaxIterPerCh = 20;
static const double kSubtractTailFactor = 8.0;
static const int    kCandidateLookahead = 220;

// フィット形状（run8で既知の代表値）
// s(t)=A*(exp(-(t-t0)/(tr+d))-exp(-(t-t0)/tr))*Theta + C, ここでは C=0 固定
static const double kTrFix[4] = {11.79, 11.23, 12.16, 11.99};
static const double kDFix[4]  = {36.94, 49.64, 24.69, 38.49};

// moduleペア化
static const int kPairMatchMaxDt = 30; // [samples]

// dt二峰探索
static const int kPromptSearchMin = -60;
static const int kPromptSearchMax = 60;
static const int kDelayedSearchMin = 80;
static const int kDelayedSearchMax = 260;
static const int kDtBandHalfWindow = 20;

// PSしきい値自動決定
static const double kPSPeakSearchMin = 140.0;
static const double kPSPeakSearchMax = 500.0;
static const int    kPSValleyBackBins = 120;
static const int    kPSValleyMarginBins = 6;
static const double kPSThresholdFallback = 80.0;

// ヒスト軸
static const int    kPSAmpBins = 1200;
static const double kPSAmpMin = 0.0;
static const double kPSAmpMax = 1200.0;

static const int    kDtBins = 600;
static const double kDtMin = -600.0;
static const double kDtMax = 600.0;

static const int    kEnergyBins = 2500;
static const double kEnergyMin = 0.0;
static const double kEnergyMax = 2.0e5;

// 描画
static const bool kUseLogY = true;
static const int  kCanvasW = 1400;
static const int  kCanvasH = 1000;

// ============================================================

struct Pulse {
  int start = -1;
  int peak = -1;
  int end = -1;
  double amp = 0.0;
};

struct FitPulse {
  int t0_init = -1;
  double t0_fit = 0.0;
  double A = 0.0;
  double C = 0.0;
  double pre_offset = 0.0;
  double I = 0.0; // = A*d
  double chi2ndf = 0.0;
  bool ok = false;
};

struct ModulePulse {
  double t = 0.0;
  double area_sum = 0.0;
};

// ============================================================
// 基本ユーティリティ
// ============================================================

static double ModeBaselineInEvent(const std::vector<double>& samples)
{
  if (samples.empty()) return 0.0;

  int vmin = (int)std::lround(samples[0]);
  int vmax = vmin;
  for (double x : samples) {
    const int v = (int)std::lround(x);
    if (v < vmin) vmin = v;
    if (v > vmax) vmax = v;
  }

  const int range = vmax - vmin + 1;
  if (range <= 0) return (double)vmin;

  std::vector<int> cnt(range, 0);
  for (double x : samples) {
    const int v = (int)std::lround(x);
    const int i = v - vmin;
    if (0 <= i && i < range) cnt[i]++;
  }

  int best_i = 0;
  int best_c = cnt[0];
  for (int i = 1; i < range; ++i) {
    if (cnt[i] > best_c) {
      best_c = cnt[i];
      best_i = i;
    }
  }
  return (double)(vmin + best_i);
}

static void BuildS(const std::vector<double>& adc, double baseline, std::vector<double>& s)
{
  s.resize(adc.size());
  for (size_t i = 0; i < adc.size(); ++i) s[i] = baseline - adc[i];
}

static double MeanInRange(const std::vector<double>& x, int i0, int i1)
{
  if (x.empty()) return 0.0;
  const int n = (int)x.size();
  i0 = std::max(0, i0);
  i1 = std::min(n - 1, i1);
  if (i1 < i0) return 0.0;

  double sum = 0.0;
  int cnt = 0;
  for (int i = i0; i <= i1; ++i) {
    sum += x[i];
    cnt++;
  }
  if (cnt <= 0) return 0.0;
  return sum / (double)cnt;
}

static double EstimateLocalPreOffset(const std::vector<double>& s, int t0)
{
  return MeanInRange(s, t0 - 80, t0 - 20);
}

static double EstimateNoiseSigma(const std::vector<double>& s, int t0)
{
  const int i0 = std::max(0, t0 - kFitHalfWindow);
  const int i1 = std::max(i0, t0 - 20);
  if (i1 <= i0) return 2.0;

  double mu = 0.0;
  int n = 0;
  for (int i = i0; i <= i1; ++i) {
    mu += s[i];
    n++;
  }
  if (n <= 1) return 2.0;
  mu /= (double)n;

  double v = 0.0;
  for (int i = i0; i <= i1; ++i) {
    const double d = s[i] - mu;
    v += d * d;
  }
  v /= (double)(n - 1);
  if (!std::isfinite(v) || v <= 0.0) return 2.0;
  return std::sqrt(v);
}

// ============================================================
// パルス抽出（PS用）
// ============================================================

static std::vector<Pulse> FindPulsesSimple(const std::vector<double>& s,
                                           double seed_thr,
                                           double end_thr,
                                           int min_sep,
                                           int min_width)
{
  std::vector<Pulse> out;
  const int n = (int)s.size();
  if (n < 3) return out;

  bool in = false;
  int dead = 0;
  Pulse cur;
  double amax = -1.0;

  for (int i = 1; i < n; ++i) {
    if (dead > 0) dead--;
    const double y = s[i];

    if (!in) {
      if (dead == 0 && s[i - 1] < seed_thr && y >= seed_thr) {
        in = true;
        cur = Pulse{};
        cur.start = i;
        cur.peak = i;
        amax = y;
      }
      continue;
    }

    if (y > amax) {
      amax = y;
      cur.peak = i;
    }

    if (y <= end_thr) {
      cur.end = i;
      cur.amp = std::max(0.0, amax);
      const int width = cur.end - cur.start + 1;
      if (width >= min_width && cur.amp > 0.0) out.push_back(cur);
      in = false;
      dead = min_sep;
      amax = -1.0;
    }
  }

  if (in) {
    cur.end = n - 1;
    cur.amp = std::max(0.0, amax);
    const int width = cur.end - cur.start + 1;
    if (width >= min_width && cur.amp > 0.0) out.push_back(cur);
  }

  return out;
}

// ============================================================
// NaI 関数フィット（hist_NaI_fitint_clean.C と同じ関数形）
// ============================================================

static double PulseKernelNoOffset(double x, double t0, double tr, double d)
{
  const double tf = tr + d;
  if (tr <= 0.0 || d <= 0.0 || tf <= 0.0) return 0.0;
  const double t = x - t0;
  if (t < 0.0) return 0.0;
  return std::exp(-t / tf) - std::exp(-t / tr);
}

static bool FitFixedBaselineByGrid(const std::vector<double>& x,
                                   const std::vector<double>& y,
                                   const std::vector<double>& ey,
                                   double tr,
                                   double d,
                                   double t0_center,
                                   double t0_range,
                                   double& t0_best,
                                   double& A_best,
                                   double& chi2ndf_best)
{
  const int n = (int)x.size();
  if (n <= 4 || y.size() != x.size() || ey.size() != x.size()) return false;
  if (tr <= 0.0 || d <= 0.0) return false;

  auto eval_at_t0 = [&](double t0, double& A, double& chi2) -> bool {
    double s1 = 0.0;
    double s2 = 0.0;
    for (int i = 0; i < n; ++i) {
      const double sigma = (ey[i] > 0.0) ? ey[i] : 1.0;
      const double w = 1.0 / (sigma * sigma);
      const double g = PulseKernelNoOffset(x[i], t0, tr, d);
      s1 += y[i] * g * w;
      s2 += g * g * w;
    }
    if (s2 <= 0.0 || !std::isfinite(s2)) return false;

    A = s1 / s2;
    if (!std::isfinite(A)) return false;
    if (A < 0.0) A = 0.0;

    chi2 = 0.0;
    for (int i = 0; i < n; ++i) {
      const double sigma = (ey[i] > 0.0) ? ey[i] : 1.0;
      const double yhat = A * PulseKernelNoOffset(x[i], t0, tr, d);
      const double r = (y[i] - yhat) / sigma;
      chi2 += r * r;
    }
    return std::isfinite(chi2);
  };

  const double tmin = t0_center - t0_range;
  const double tmax = t0_center + t0_range;
  const int ngrid_coarse = 25;
  const int ngrid_fine = 19;

  bool found = false;
  double best_chi2 = 0.0;
  double best_t0 = t0_center;
  double best_A = 0.0;

  for (int ig = 0; ig < ngrid_coarse; ++ig) {
    const double t0 = tmin + (tmax - tmin) * (double)ig / (double)(ngrid_coarse - 1);
    double A = 0.0;
    double chi2 = 0.0;
    if (!eval_at_t0(t0, A, chi2)) continue;
    if (!found || chi2 < best_chi2) {
      found = true;
      best_chi2 = chi2;
      best_t0 = t0;
      best_A = A;
    }
  }
  if (!found) return false;

  const double coarse_step = (tmax - tmin) / (double)(ngrid_coarse - 1);
  const double fmin = std::max(tmin, best_t0 - 2.0 * coarse_step);
  const double fmax = std::min(tmax, best_t0 + 2.0 * coarse_step);
  for (int ig = 0; ig < ngrid_fine; ++ig) {
    const double t0 = fmin + (fmax - fmin) * (double)ig / (double)(ngrid_fine - 1);
    double A = 0.0;
    double chi2 = 0.0;
    if (!eval_at_t0(t0, A, chi2)) continue;
    if (chi2 < best_chi2) {
      best_chi2 = chi2;
      best_t0 = t0;
      best_A = A;
    }
  }

  const double ndf = (double)(n - 2);
  chi2ndf_best = (ndf > 0.0) ? (best_chi2 / ndf) : best_chi2;
  t0_best = best_t0;
  A_best = best_A;
  return std::isfinite(chi2ndf_best);
}

static std::vector<int> FindStartCandidates(const std::vector<double>& s,
                                            double seed_thr,
                                            double end_thr,
                                            int min_sep_samples)
{
  std::vector<int> starts;
  const int n = (int)s.size();
  if (n <= 1) return starts;

  int dead = 0;
  for (int i = 1; i < n; ++i) {
    if (dead > 0) {
      dead--;
      continue;
    }

    const bool crossed = (s[i - 1] < seed_thr && s[i] >= seed_thr);
    if (!crossed) continue;

    starts.push_back(i);

    int j = i + 1;
    while (j < n && s[j] >= end_thr) j++;
    const int width_dead = std::max(0, j - i);
    dead = std::max(min_sep_samples, width_dead);
  }
  return starts;
}

static void AddPeakSupplementCandidates(const std::vector<double>& s,
                                        std::vector<int>& starts)
{
  // 閾値crossだけでは取りこぼす重なりを補うため、局所ピークから追加候補を作る。
  const int n = (int)s.size();
  if (n < 5) return;

  auto has_near = [&](int t0, int sep) {
    for (int x : starts) {
      if (std::abs(x - t0) <= sep) return true;
    }
    return false;
  };

  for (int i = 2; i <= n - 3; ++i) {
    if (!(s[i] >= s[i - 1] && s[i] > s[i + 1])) continue;
    if (s[i] < kNaISeedThr + 10.0) continue;

    int t0 = i;
    // 立ち上がり側へ逆探索
    const int jmin = std::max(1, i - 220);
    bool found = false;
    for (int j = i; j >= jmin; --j) {
      if (s[j - 1] < kNaISeedThr && s[j] >= kNaISeedThr) {
        t0 = j;
        found = true;
        break;
      }
    }
    if (!found) t0 = std::max(1, i - 24);

    if (!has_near(t0, 10)) starts.push_back(t0);
  }
}

static double CandidatePeakScore(const std::vector<double>& s, int t0)
{
  const int n = (int)s.size();
  const int i0 = std::max(0, t0);
  const int i1 = std::min(n - 1, t0 + kCandidateLookahead);
  double peak = 0.0;
  for (int i = i0; i <= i1; ++i) {
    if (s[i] > peak) peak = s[i];
  }
  return peak;
}

static FitPulse FitOnePulse(const std::vector<double>& s,
                            int t0_init,
                            double tr_fix,
                            double d_fix)
{
  FitPulse fp;
  fp.t0_init = t0_init;

  const int n = (int)s.size();
  if (n <= 0) return fp;

  const int i0 = t0_init - kFitPreWindow;
  const int i1 = t0_init + kFitHalfWindow;
  if (i0 < 0 || i1 >= n) return fp;

  const int npts = i1 - i0 + 1;
  std::vector<double> x(npts), y(npts), ey(npts);

  const double pre_offset = EstimateLocalPreOffset(s, t0_init);
  const double sigma_noise = EstimateNoiseSigma(s, t0_init);

  double s_peak = 0.0;
  for (int i = 0; i < npts; ++i) {
    const int idx = i0 + i;
    x[i] = (double)idx;
    y[i] = s[idx] - pre_offset;
    if (y[i] > s_peak) s_peak = y[i];
  }
  if (s_peak <= 0.0) return fp;

  for (int i = 0; i < npts; ++i) {
    const double ys = std::max(0.0, y[i]);
    const double sig_sys = kRelShapeSys * ys;
    ey[i] = std::sqrt(sigma_noise * sigma_noise + sig_sys * sig_sys);
    if (!std::isfinite(ey[i]) || ey[i] < 1.0) ey[i] = 1.0;
  }

  double t0_fit = 0.0;
  double A_fit = 0.0;
  double chi2ndf = 0.0;
  if (!FitFixedBaselineByGrid(x, y, ey, tr_fix, d_fix,
                              (double)t0_init, kT0FitRange,
                              t0_fit, A_fit, chi2ndf)) {
    return fp;
  }

  fp.t0_fit = t0_fit;
  fp.A = A_fit;
  fp.C = 0.0;
  fp.pre_offset = pre_offset;
  fp.chi2ndf = chi2ndf;
  fp.I = fp.A * d_fix;
  fp.ok = (fp.A > 0.0 && std::isfinite(fp.chi2ndf) && fp.chi2ndf <= kChi2Cut);
  return fp;
}

static void SubtractAcceptedPulse(std::vector<double>& s_work,
                                  const FitPulse& fp,
                                  double tr_fix,
                                  double d_fix)
{
  if (!fp.ok) return;
  if (fp.A <= 0.0 || tr_fix <= 0.0 || d_fix <= 0.0) return;

  const int n = (int)s_work.size();
  const double tf = tr_fix + d_fix;
  if (n <= 0 || tf <= 0.0) return;

  const int i0 = std::max(0, (int)std::floor(fp.t0_fit) - 2);
  const int i1 = std::min(n - 1, (int)std::ceil(fp.t0_fit + kSubtractTailFactor * tf));
  for (int i = i0; i <= i1; ++i) {
    const double yhat = fp.A * PulseKernelNoOffset((double)i, fp.t0_fit, tr_fix, d_fix);
    s_work[i] -= yhat;
  }
}

static bool HasNearbyT0(const std::vector<int>& t0s, int t0, int sep)
{
  for (int x : t0s) {
    if (std::abs(x - t0) <= sep) return true;
  }
  return false;
}

static std::vector<FitPulse> FitPulsesWithDecomposition(const std::vector<double>& s_input,
                                                        double tr_fix,
                                                        double d_fix)
{
  std::vector<FitPulse> accepted;
  if (s_input.empty()) return accepted;

  std::vector<double> s_work = s_input;
  std::vector<int> tried;

  for (int iter = 0; iter < kMaxIterPerCh; ++iter) {
    if ((int)accepted.size() >= kMaxAcceptedPerCh) break;

    std::vector<int> starts = FindStartCandidates(s_work, kNaISeedThr, kNaIEndThr, kNaIMinSep);
    AddPeakSupplementCandidates(s_work, starts);
    if (starts.empty()) break;

    // 重複除去
    std::sort(starts.begin(), starts.end());
    starts.erase(std::unique(starts.begin(), starts.end()), starts.end());

    int best_t0 = -1;
    double best_sc = -1.0;
    for (int t0 : starts) {
      if (t0 - kFitPreWindow < 0) continue;
      if (t0 + kFitHalfWindow >= (int)s_work.size()) continue;
      if (HasNearbyT0(tried, t0, 10)) continue;
      if (HasNearbyT0(std::vector<int>{}, t0, 0)) {}
      const double sc = CandidatePeakScore(s_work, t0);
      if (sc > best_sc) {
        best_sc = sc;
        best_t0 = t0;
      }
    }
    if (best_t0 < 0) break;

    tried.push_back(best_t0);

    FitPulse fp = FitOnePulse(s_work, best_t0, tr_fix, d_fix);
    if (!fp.ok) continue;

    bool near_dup = false;
    for (const auto& ap : accepted) {
      if (std::abs((int)std::lround(ap.t0_fit) - (int)std::lround(fp.t0_fit)) < 20) {
        near_dup = true;
        break;
      }
    }
    if (near_dup) continue;

    accepted.push_back(fp);
    SubtractAcceptedPulse(s_work, fp, tr_fix, d_fix);
  }

  std::sort(accepted.begin(), accepted.end(),
            [](const FitPulse& a, const FitPulse& b) { return a.t0_fit < b.t0_fit; });
  return accepted;
}

// ============================================================
// PS閾値自動決定
// ============================================================

static double SmoothedBinContent(const TH1D* h, int b, int hw)
{
  const int n = h->GetNbinsX();
  const int b0 = std::max(1, b - hw);
  const int b1 = std::min(n, b + hw);
  double s = 0.0;
  int c = 0;
  for (int i = b0; i <= b1; ++i) {
    s += h->GetBinContent(i);
    c++;
  }
  if (c <= 0) return 0.0;
  return s / (double)c;
}

static double FindPSValleyThreshold(const TH1D* h,
                                    int& peak_bin_out,
                                    int& valley_bin_out)
{
  peak_bin_out = -1;
  valley_bin_out = -1;
  if (!h) return kPSThresholdFallback;

  const int bpk0 = h->GetXaxis()->FindBin(kPSPeakSearchMin);
  const int bpk1 = h->GetXaxis()->FindBin(kPSPeakSearchMax);
  const int nbin = h->GetNbinsX();

  std::vector<double> ys(nbin + 1, 0.0);
  std::vector<double> ys_pk(nbin + 1, 0.0);
  for (int b = 1; b <= nbin; ++b) {
    ys[b] = SmoothedBinContent(h, b, 3);     // 谷探索用（細め）
    ys_pk[b] = SmoothedBinContent(h, b, 6);  // 山探索用（太め）
  }

  // まず探索窓内の global max を取って、局所最大の下限を決める
  double ymax_win = 0.0;
  for (int b = std::min(bpk0, bpk1); b <= std::max(bpk0, bpk1); ++b) {
    if (ys_pk[b] > ymax_win) ymax_win = ys_pk[b];
  }
  const double local_peak_min = 0.20 * ymax_win;

  // 低エネルギー端の単調減少成分ではなく、Michel側の局所ピークを拾う
  int bpeak = -1;
  double ypeak = -1.0;
  double best_prom = 1.05;
  for (int b = std::min(bpk0, bpk1) + 1; b <= std::max(bpk0, bpk1) - 1; ++b) {
    if (ys_pk[b] < local_peak_min) continue;
    if (ys_pk[b] >= ys_pk[b - 1] && ys_pk[b] > ys_pk[b + 1]) {
      const int bleft0 = std::max(1, b - 48);
      const int bleft1 = std::max(bleft0, b - 20);
      double yleft = 0.0;
      int cnt = 0;
      for (int j = bleft0; j <= bleft1; ++j) { yleft += ys_pk[j]; cnt++; }
      if (cnt > 0) yleft /= (double)cnt;
      const double prom = ys_pk[b] / std::max(1.0, yleft);
      if (prom > best_prom || (std::abs(prom - best_prom) < 1e-9 && ys_pk[b] > ypeak)) {
        best_prom = prom;
        ypeak = ys_pk[b];
        bpeak = b;
      }
    }
  }

  // 局所最大が見つからない場合は窓内最大（ただし低xを避ける）
  if (bpeak < 1) {
    int bf0 = h->GetXaxis()->FindBin(140.0);
    int bf1 = std::max(bf0, std::max(bpk0, bpk1));
    for (int b = std::min(bf0, bf1); b <= std::max(bf0, bf1); ++b) {
      if (ys_pk[b] > ypeak) {
        ypeak = ys_pk[b];
        bpeak = b;
      }
    }
  }
  if (bpeak < 1) return kPSThresholdFallback;

  const int bmin = std::max(h->GetXaxis()->FindBin(40.0), bpeak - kPSValleyBackBins);
  const int bmax = std::max(bmin, bpeak - kPSValleyMarginBins);

  // 谷候補（局所最小）のうち、最も低い点を採用する。
  // これにより「山直前の浅い凹み」ではなく、見た目の谷を選びやすくする。
  int bvalley = -1;
  double yvalley = std::numeric_limits<double>::infinity();
  for (int b = bmin + 1; b <= bmax - 1; ++b) {
    if (ys[b] <= ys[b - 1] && ys[b] < ys[b + 1]) {
      if (ys[b] < yvalley) {
        yvalley = ys[b];
        bvalley = b;
      }
    }
  }

  // 局所谷が無い場合は、peak直前の探索窓内で最小値を採用。
  if (bvalley < 1) {
    for (int b = bmin; b <= bmax; ++b) {
      if (ys[b] < yvalley) {
        yvalley = ys[b];
        bvalley = b;
      }
    }
  }

  if (bvalley < 1) return kPSThresholdFallback;
  peak_bin_out = bpeak;
  valley_bin_out = bvalley;
  const double thr = h->GetBinCenter(bvalley);
  return std::max(60.0, thr);
}

// moduleペア化
// ============================================================

static std::vector<ModulePulse> BuildModulePairs(const std::vector<FitPulse>& c1,
                                                 const std::vector<FitPulse>& c2,
                                                 int pair_dt_max)
{
  std::vector<ModulePulse> out;
  if (c1.empty() || c2.empty()) return out;

  std::vector<int> i1(c1.size()), i2(c2.size());
  for (size_t i = 0; i < c1.size(); ++i) i1[i] = (int)i;
  for (size_t i = 0; i < c2.size(); ++i) i2[i] = (int)i;

  std::sort(i1.begin(), i1.end(), [&](int a, int b) { return c1[a].t0_fit < c1[b].t0_fit; });
  std::sort(i2.begin(), i2.end(), [&](int a, int b) { return c2[a].t0_fit < c2[b].t0_fit; });

  std::vector<char> used2(c2.size(), 0);
  for (int ia : i1) {
    int jbest = -1;
    double dmin = 1e30;
    for (int ib : i2) {
      if (used2[ib]) continue;
      const double d = std::abs(c1[ia].t0_fit - c2[ib].t0_fit);
      if (d <= (double)pair_dt_max && d < dmin) {
        dmin = d;
        jbest = ib;
      }
    }
    if (jbest < 0) continue;
    used2[jbest] = 1;

    ModulePulse m;
    m.t = 0.5 * (c1[ia].t0_fit + c2[jbest].t0_fit);
    m.area_sum = c1[ia].I + c2[jbest].I;
    out.push_back(m);
  }

  std::sort(out.begin(), out.end(), [](const ModulePulse& a, const ModulePulse& b) {
    return a.t < b.t;
  });
  return out;
}

static int FindNearestModuleIndex(const std::vector<ModulePulse>& v, double tps)
{
  if (v.empty()) return -1;
  int ibest = -1;
  double dmin = 1e30;
  for (size_t i = 0; i < v.size(); ++i) {
    const double d = std::abs(v[i].t - tps);
    if (d < dmin) {
      dmin = d;
      ibest = (int)i;
    }
  }
  return ibest;
}

static double FindPeakXInRange(TH1D* h, double xmin, double xmax)
{
  if (!h) return 0.0;
  int b0 = h->GetXaxis()->FindBin(xmin);
  int b1 = h->GetXaxis()->FindBin(xmax);
  if (b1 < b0) std::swap(b0, b1);
  int bbest = b0;
  double ymax = -1.0;
  for (int b = b0; b <= b1; ++b) {
    const double y = h->GetBinContent(b);
    if (y > ymax) {
      ymax = y;
      bbest = b;
    }
  }
  return h->GetBinCenter(bbest);
}

// ============================================================

void study_michel_ps_tagged()
{
  gStyle->SetOptStat(1110);
  gStyle->SetTitleFontSize(0.04);

  gSystem->mkdir(kOutputDir, true);
  const TString out_pdf = Form("%s/study_michel_ps_tagged_%s.pdf", kOutputDir, kInputTag);

  // ---------- Pass 1: PS波高ヒスト作成（しきい値自動決定用） ----------
  TH1D* hPSAmp[2] = {
    new TH1D("hPSAmp_A", "PS_A pulse amplitude (simple, all pulses);s = baseline - ADC [counts];Pulses",
             kPSAmpBins, kPSAmpMin, kPSAmpMax),
    new TH1D("hPSAmp_B", "PS_B pulse amplitude (simple, all pulses);s = baseline - ADC [counts];Pulses",
             kPSAmpBins, kPSAmpMin, kPSAmpMax)
  };
  long n_ps_preselected[2] = {0, 0}; // ここでは「しきい値推定用に詰めたパルス数」

  {
    TString fps[2];
    for (int i = 0; i < 2; ++i) fps[i] = Form("%s/%s", kInputDir, kPSFiles[i]);

    std::ifstream in_ps[2];
    for (int i = 0; i < 2; ++i) {
      in_ps[i].open(fps[i].Data());
      if (!in_ps[i]) {
        std::cerr << "[study] ERROR: cannot open " << fps[i] << std::endl;
        return;
      }
    }

    std::vector<double> buf_ps[2];
    for (int i = 0; i < 2; ++i) buf_ps[i].reserve(kSamplesPerEvent);

    while (true) {
      double xps[2];
      bool ok = true;
      for (int i = 0; i < 2; ++i) {
        if (!(in_ps[i] >> xps[i])) {
          ok = false;
          break;
        }
      }
      if (!ok) break;

      for (int i = 0; i < 2; ++i) buf_ps[i].push_back(xps[i]);
      if ((int)buf_ps[0].size() < kSamplesPerEvent) continue;

      // check/hist_PS_peak.C と同じ簡易条件で、PSパルスを全て積む
      for (int ip = 0; ip < 2; ++ip) {
        const double bps = ModeBaselineInEvent(buf_ps[ip]);
        std::vector<double> sps;
        BuildS(buf_ps[ip], bps, sps);
        const std::vector<Pulse> ps = FindPulsesSimple(sps, kPSSeedThr, kPSEndThr, kPSMinSep, kPSMinWidth);
        for (const auto& q : ps) {
          if (q.amp <= 0.0) continue;
          hPSAmp[ip]->Fill(q.amp);
          n_ps_preselected[ip]++;
        }
      }

      for (int i = 0; i < 2; ++i) buf_ps[i].clear();
    }
  }

  int ps_peak_bin[2] = {-1, -1};
  int ps_valley_bin[2] = {-1, -1};
  double ps_thr_auto[2] = {kPSThresholdFallback, kPSThresholdFallback};
  for (int i = 0; i < 2; ++i) {
    ps_thr_auto[i] = FindPSValleyThreshold(hPSAmp[i], ps_peak_bin[i], ps_valley_bin[i]);
  }
  std::cout << "[study] Pass1 entries A/B = "
            << (long)hPSAmp[0]->GetEntries() << " / " << (long)hPSAmp[1]->GetEntries() << "\n";
  std::cout << "[study] Pass1 counted pulses A/B = "
            << n_ps_preselected[0] << " / " << n_ps_preselected[1] << "\n";
  std::cout << "[study] PS thr auto A/B = "
            << ps_thr_auto[0] << " / " << ps_thr_auto[1] << "\n";

  // ---------- Pass 2: PSしきい値をかけて NaI fit + dt ----------

  TH1D* hDtA = new TH1D("hDtA", "ModuleA nearest pair: dt=t_{mod}-t_{PS};dt [samples];Events",
                        kDtBins, kDtMin, kDtMax);
  TH1D* hDtB = new TH1D("hDtB", "ModuleB nearest pair: dt=t_{mod}-t_{PS};dt [samples];Events",
                        kDtBins, kDtMin, kDtMax);

  std::vector<double> dtA_vals;
  std::vector<double> dtB_vals;
  std::vector<double> EA_vals;
  std::vector<double> EB_vals;
  dtA_vals.reserve(20000);
  dtB_vals.reserve(6000);
  EA_vals.reserve(20000);
  EB_vals.reserve(6000);

  long n_events = 0;
  long n_ps_selected = 0;

  {
    TString fps[2], fna[4];
    for (int i = 0; i < 2; ++i) fps[i] = Form("%s/%s", kInputDir, kPSFiles[i]);
    for (int i = 0; i < 4; ++i) fna[i] = Form("%s/%s", kInputDir, kNaIFiles[i]);

    std::ifstream in_ps[2];
    std::ifstream in_na[4];
    for (int i = 0; i < 2; ++i) {
      in_ps[i].open(fps[i].Data());
      if (!in_ps[i]) {
        std::cerr << "[study] ERROR: cannot open " << fps[i] << std::endl;
        return;
      }
    }
    for (int i = 0; i < 4; ++i) {
      in_na[i].open(fna[i].Data());
      if (!in_na[i]) {
        std::cerr << "[study] ERROR: cannot open " << fna[i] << std::endl;
        return;
      }
    }

    std::vector<double> buf_ps[2], buf_na[4];
    for (int i = 0; i < 2; ++i) buf_ps[i].reserve(kSamplesPerEvent);
    for (int i = 0; i < 4; ++i) buf_na[i].reserve(kSamplesPerEvent);

    while (true) {
      double xps[2], xna[4];
      bool ok = true;
      for (int i = 0; i < 2; ++i) {
        if (!(in_ps[i] >> xps[i])) { ok = false; break; }
      }
      if (!ok) break;
      for (int i = 0; i < 4; ++i) {
        if (!(in_na[i] >> xna[i])) { ok = false; break; }
      }
      if (!ok) break;

      for (int i = 0; i < 2; ++i) buf_ps[i].push_back(xps[i]);
      for (int i = 0; i < 4; ++i) buf_na[i].push_back(xna[i]);
      if ((int)buf_ps[0].size() < kSamplesPerEvent) continue;

      n_events++;

      // PS: しきい値を超える最早パルスを採用
      bool has_ps = false;
      double tps = -1.0;
      double aps = 0.0;
      for (int ip = 0; ip < 2; ++ip) {
        const double b = ModeBaselineInEvent(buf_ps[ip]);
        std::vector<double> s;
        BuildS(buf_ps[ip], b, s);
        const std::vector<Pulse> p = FindPulsesSimple(s, kPSSeedThr, kPSEndThr, kPSMinSep, kPSMinWidth);
        for (const auto& q : p) {
          if (q.amp < ps_thr_auto[ip]) continue;
          if (!has_ps || (double)q.start < tps) {
            has_ps = true;
            tps = (double)q.start;
            aps = q.amp;
          }
        }
      }

      if (has_ps) {
        n_ps_selected++;

        // NaI: 各chで関数フィット分解
        std::vector<FitPulse> fp[4];
        for (int ch = 0; ch < 4; ++ch) {
          const double b = ModeBaselineInEvent(buf_na[ch]);
          std::vector<double> s;
          BuildS(buf_na[ch], b, s);
          fp[ch] = FitPulsesWithDecomposition(s, kTrFix[ch], kDFix[ch]);
        }

        // moduleペア
        const std::vector<ModulePulse> modA = BuildModulePairs(fp[0], fp[1], kPairMatchMaxDt);
        const std::vector<ModulePulse> modB = BuildModulePairs(fp[2], fp[3], kPairMatchMaxDt);

        const int ia = FindNearestModuleIndex(modA, tps);
        if (ia >= 0) {
          const double dt = modA[ia].t - tps;
          hDtA->Fill(dt);
          dtA_vals.push_back(dt);
          EA_vals.push_back(modA[ia].area_sum);
        }

        const int ib = FindNearestModuleIndex(modB, tps);
        if (ib >= 0) {
          const double dt = modB[ib].t - tps;
          hDtB->Fill(dt);
          dtB_vals.push_back(dt);
          EB_vals.push_back(modB[ib].area_sum);
        }
      }

      for (int i = 0; i < 2; ++i) buf_ps[i].clear();
      for (int i = 0; i < 4; ++i) buf_na[i].clear();
    }
  }

  // dt二峰を自動決定
  const double dtA_prompt = FindPeakXInRange(hDtA, kPromptSearchMin, kPromptSearchMax);
  const double dtA_delay  = FindPeakXInRange(hDtA, kDelayedSearchMin, kDelayedSearchMax);
  const double dtB_prompt = FindPeakXInRange(hDtB, kPromptSearchMin, kPromptSearchMax);
  const double dtB_delay  = FindPeakXInRange(hDtB, kDelayedSearchMin, kDelayedSearchMax);

  // A/Bの二峰別エネルギーヒスト（要求に合わせてA,Bの2系統のみ）
  TH1D* hEA_prompt = new TH1D("hEA_prompt",
                              Form("ModuleA energy (dt=%.1f#pm%d);I_{A1}+I_{A2} [ADC*samples];Pairs",
                                   dtA_prompt, kDtBandHalfWindow),
                              kEnergyBins, kEnergyMin, kEnergyMax);
  TH1D* hEA_delay  = new TH1D("hEA_delay",
                              Form("ModuleA energy (dt=%.1f#pm%d);I_{A1}+I_{A2} [ADC*samples];Pairs",
                                   dtA_delay, kDtBandHalfWindow),
                              kEnergyBins, kEnergyMin, kEnergyMax);

  TH1D* hEB_prompt = new TH1D("hEB_prompt",
                              Form("ModuleB energy (dt=%.1f#pm%d);I_{B1}+I_{B2} [ADC*samples];Pairs",
                                   dtB_prompt, kDtBandHalfWindow),
                              kEnergyBins, kEnergyMin, kEnergyMax);
  TH1D* hEB_delay  = new TH1D("hEB_delay",
                              Form("ModuleB energy (dt=%.1f#pm%d);I_{B1}+I_{B2} [ADC*samples];Pairs",
                                   dtB_delay, kDtBandHalfWindow),
                              kEnergyBins, kEnergyMin, kEnergyMax);

  for (size_t i = 0; i < EA_vals.size(); ++i) {
    const double dt = dtA_vals[i];
    if (std::abs(dt - dtA_prompt) <= (double)kDtBandHalfWindow) hEA_prompt->Fill(EA_vals[i]);
    if (std::abs(dt - dtA_delay) <= (double)kDtBandHalfWindow) hEA_delay->Fill(EA_vals[i]);
  }
  for (size_t i = 0; i < EB_vals.size(); ++i) {
    const double dt = dtB_vals[i];
    if (std::abs(dt - dtB_prompt) <= (double)kDtBandHalfWindow) hEB_prompt->Fill(EB_vals[i]);
    if (std::abs(dt - dtB_delay) <= (double)kDtBandHalfWindow) hEB_delay->Fill(EB_vals[i]);
  }

  // ---------- 描画 ----------

  TCanvas* c = new TCanvas("c_study_ps_tagged", "study ps tagged", kCanvasW, kCanvasH);
  c->Print(out_pdf + "[");

  // meta
  c->Clear();
  TPaveText* pt = new TPaveText(0.05, 0.05, 0.95, 0.95, "NDC");
  pt->SetFillStyle(0);
  pt->SetBorderSize(0);
  pt->SetTextAlign(12);
  pt->SetTextSize(0.028);
  pt->AddText("=== p2MEG PS-tagged study (fit integral) ===");
  pt->AddText(Form("InputDir : %s", kInputDir));
  pt->AddText(Form("Tag      : %s", kInputTag));
  pt->AddText(Form("Events   : %ld", n_events));
  pt->AddText(Form("Pass1 simple cond (same as hist_PS_peak.C): seed=%.1f, end=%.1f, min_sep=%d",
                   kPSSeedThr, kPSEndThr, kPSMinSep));
  pt->AddText(Form("Pass1 filled pulses A/B = %ld / %ld",
                   (long)hPSAmp[0]->GetEntries(), (long)hPSAmp[1]->GetEntries()));
  pt->AddText(Form("PS_A auto threshold (valley) = %.1f counts", ps_thr_auto[0]));
  pt->AddText(Form("PS_B auto threshold (valley) = %.1f counts", ps_thr_auto[1]));
  pt->AddText(Form("PS selected events = %ld (%.1f%%)", n_ps_selected,
                   100.0 * (double)n_ps_selected / std::max(1L, n_events)));
  pt->AddText("NaI model: C + A*(exp(-(t-t0)/(tr+d)) - exp(-(t-t0)/tr))*Theta, C fixed to 0");
  pt->AddText("NaI integral: I = A*d (fit area, ADC*samples)");
  pt->AddText(Form("dt peaks A: prompt=%.1f, delayed=%.1f samp", dtA_prompt, dtA_delay));
  pt->AddText(Form("dt peaks B: prompt=%.1f, delayed=%.1f samp", dtB_prompt, dtB_delay));
  pt->AddText(Form("A entries prompt/delayed = %ld / %ld",
                   (long)hEA_prompt->GetEntries(), (long)hEA_delay->GetEntries()));
  pt->AddText(Form("B entries prompt/delayed = %ld / %ld",
                   (long)hEB_prompt->GetEntries(), (long)hEB_delay->GetEntries()));
  pt->Draw();
  c->Print(out_pdf);
  delete pt;

  // PS波高(A/B) + しきい値
  c->Clear();
  c->Divide(2, 1);
  for (int i = 0; i < 2; ++i) {
    c->cd(i + 1);
    if (kUseLogY) gPad->SetLogy();
    gPad->SetGridx();
    gPad->SetGridy();
    hPSAmp[i]->SetLineColor((i == 0) ? (kBlue + 1) : (kRed + 1));
    hPSAmp[i]->SetLineWidth(2);
    hPSAmp[i]->Draw("hist");

    const double ymax = hPSAmp[i]->GetMaximum();
    TLine* lthr = new TLine(ps_thr_auto[i], 0.0, ps_thr_auto[i], ymax);
    lthr->SetLineColor(kBlack);
    lthr->SetLineWidth(2);
    lthr->SetLineStyle(2);
    lthr->Draw();

    if (ps_peak_bin[i] > 0) {
      const double xpk = hPSAmp[i]->GetBinCenter(ps_peak_bin[i]);
      TLine* lpk = new TLine(xpk, 0.0, xpk, ymax);
      lpk->SetLineColor(kMagenta + 2);
      lpk->SetLineStyle(3);
      lpk->Draw();
    }

    TLatex t;
    t.SetNDC(true);
    t.SetTextSize(0.04);
    t.DrawLatex(0.14, 0.86, Form("PS_%c threshold(auto valley)=%.1f", (i == 0 ? 'A' : 'B'), ps_thr_auto[i]));
    t.SetTextSize(0.032);
    t.DrawLatex(0.14, 0.80, Form("entries=%ld (simple all pulses)", (long)hPSAmp[i]->GetEntries()));
  }
  c->Print(out_pdf);

  // dt二峰
  c->Clear();
  c->Divide(2, 1);
  c->cd(1);
  gPad->SetGridx(); gPad->SetGridy();
  hDtA->SetLineColor(kBlue + 1);
  hDtA->SetLineWidth(2);
  hDtA->Draw("hist");
  {
    const double ymax = hDtA->GetMaximum();
    TLine* lp = new TLine(dtA_prompt, 0.0, dtA_prompt, ymax);
    TLine* ld = new TLine(dtA_delay, 0.0, dtA_delay, ymax);
    lp->SetLineColor(kRed + 1); lp->SetLineStyle(2); lp->SetLineWidth(2);
    ld->SetLineColor(kMagenta + 2); ld->SetLineStyle(2); ld->SetLineWidth(2);
    lp->Draw(); ld->Draw();
  }

  c->cd(2);
  gPad->SetGridx(); gPad->SetGridy();
  hDtB->SetLineColor(kBlue + 1);
  hDtB->SetLineWidth(2);
  hDtB->Draw("hist");
  {
    const double ymax = hDtB->GetMaximum();
    TLine* lp = new TLine(dtB_prompt, 0.0, dtB_prompt, ymax);
    TLine* ld = new TLine(dtB_delay, 0.0, dtB_delay, ymax);
    lp->SetLineColor(kRed + 1); lp->SetLineStyle(2); lp->SetLineWidth(2);
    ld->SetLineColor(kMagenta + 2); ld->SetLineStyle(2); ld->SetLineWidth(2);
    lp->Draw(); ld->Draw();
  }
  c->Print(out_pdf);

  // NaIエネルギー（A,Bの2つ）
  c->Clear();
  c->Divide(2, 1);

  c->cd(1);
  if (kUseLogY) gPad->SetLogy();
  gPad->SetGridx(); gPad->SetGridy();
  hEA_prompt->SetLineColor(kRed + 1); hEA_prompt->SetLineWidth(2);
  hEA_delay->SetLineColor(kBlue + 1); hEA_delay->SetLineWidth(2);
  hEA_prompt->Draw("hist");
  hEA_delay->Draw("hist same");
  {
    TLegend* lg = new TLegend(0.56, 0.74, 0.92, 0.90);
    lg->SetBorderSize(0);
    lg->SetFillStyle(0);
    lg->AddEntry(hEA_prompt, Form("A prompt dt=%.1f#pm%d", dtA_prompt, kDtBandHalfWindow), "l");
    lg->AddEntry(hEA_delay,  Form("A delayed dt=%.1f#pm%d", dtA_delay, kDtBandHalfWindow), "l");
    lg->Draw();
  }

  c->cd(2);
  if (kUseLogY) gPad->SetLogy();
  gPad->SetGridx(); gPad->SetGridy();
  hEB_prompt->SetLineColor(kRed + 1); hEB_prompt->SetLineWidth(2);
  hEB_delay->SetLineColor(kBlue + 1); hEB_delay->SetLineWidth(2);
  hEB_prompt->Draw("hist");
  hEB_delay->Draw("hist same");
  {
    TLegend* lg = new TLegend(0.56, 0.74, 0.92, 0.90);
    lg->SetBorderSize(0);
    lg->SetFillStyle(0);
    lg->AddEntry(hEB_prompt, Form("B prompt dt=%.1f#pm%d", dtB_prompt, kDtBandHalfWindow), "l");
    lg->AddEntry(hEB_delay,  Form("B delayed dt=%.1f#pm%d", dtB_delay, kDtBandHalfWindow), "l");
    lg->Draw();
  }

  c->Print(out_pdf);

  c->Print(out_pdf + "]");
  delete c;

  std::cout << "[study] output pdf: " << out_pdf << std::endl;
}
