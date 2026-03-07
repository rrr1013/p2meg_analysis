// check/michel_stepwise_hist_based_run8.C
//
// 目的:
//   「ヒストグラムを順に見て閾値を決める」手順で run8 を再解析し、
//   最後に Michel スペクトルとの一致度を評価する。
//
// 決定順序:
//   A) PS振幅ヒスト -> PS hit閾値
//   B) NaI |dt12|ヒスト -> NaIペア時刻差カット
//   C) dt=t_mod-t_PS ヒスト -> prompt時間窓
//   D) prompt vs sideband のEヒスト -> NaIエネルギー閾値
//   E) sideband差分Eヒストに Michel 形状を重ねて一致度評価
//
// 出力:
//   doc/mainexp/michel_stepwise_hist_based_run8.pdf
//
// 実行:
//   root -l -b -q 'check/michel_stepwise_hist_based_run8.C'

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TF1.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TSystem.h"

// 既存の再構成ロジック・定数を再利用
#include "study_michel_ps_tagged.C"

struct SideCandidate {
  bool ok = false;
  double dt = 0.0;
  double e_sum = 0.0;
  double ps_amp = 0.0;
};

struct GausInfo {
  bool ok = false;
  double mean = 0.0;
  double sigma = 0.0;
  double chi2ndf = 0.0;
};

struct MichelFitInfo {
  bool fit_ok = false;
  double chi2ndf = 0.0;
  double r2 = 0.0;
  double n = 0.0;
  double e0 = 0.0;
  double w = 0.0;
  double c = 0.0;
};

static int FindNearestPulseIndexStepwise(const std::vector<FitPulse>& v, double tref)
{
  if (v.empty()) return -1;
  int ibest = -1;
  double dmin = 1.0e30;
  for (size_t i = 0; i < v.size(); ++i) {
    const double d = std::abs(v[i].t0_fit - tref);
    if (d < dmin) {
      dmin = d;
      ibest = (int)i;
    }
  }
  return ibest;
}

static double QuantileFromSortedStepwise(const std::vector<double>& x_sorted, double q)
{
  if (x_sorted.empty()) return 0.0;
  if (q <= 0.0) return x_sorted.front();
  if (q >= 1.0) return x_sorted.back();
  const double pos = q * (double)(x_sorted.size() - 1);
  const int i0 = (int)std::floor(pos);
  const int i1 = std::min((int)x_sorted.size() - 1, i0 + 1);
  const double t = pos - (double)i0;
  return (1.0 - t) * x_sorted[i0] + t * x_sorted[i1];
}

static double FractionLeqFromSortedStepwise(const std::vector<double>& x_sorted, double cut)
{
  if (x_sorted.empty()) return 0.0;
  const auto it = std::upper_bound(x_sorted.begin(), x_sorted.end(), cut);
  const double n = (double)std::distance(x_sorted.begin(), it);
  return n / (double)x_sorted.size();
}

// NaI clean 条件:
//   fit 時刻の近傍(±kFitHalfWindow)に NaI 候補開始点が 1本だけ存在する事象のみ採用する。
//   条件は hist_NaI_fitint_clean.C の「single candidate in fit window」と同じ意図。
static std::vector<int> BuildNaIStartCandidatesStepwise(const std::vector<double>& s)
{
  std::vector<int> starts = FindStartCandidates(s, kNaISeedThr, kNaIEndThr, kNaIMinSep);
  AddPeakSupplementCandidates(s, starts);
  if (starts.empty()) return starts;
  std::sort(starts.begin(), starts.end());
  starts.erase(std::unique(starts.begin(), starts.end()), starts.end());
  return starts;
}

static bool IsSingleCandidateInWindowStepwise(const std::vector<int>& starts,
                                              int t0,
                                              int halfwin)
{
  int n_in = 0;
  for (int s0 : starts) {
    if (std::abs(s0 - t0) <= halfwin) n_in++;
    if (n_in >= 2) return false;
  }
  return (n_in == 1);
}

static void ApplyNaICleanSelectionStepwise(const std::vector<double>& s,
                                           std::vector<FitPulse>& pulses,
                                           long& n_before,
                                           long& n_after)
{
  n_before += (long)pulses.size();
  if (pulses.empty()) return;

  const std::vector<int> starts = BuildNaIStartCandidatesStepwise(s);
  if (starts.empty()) {
    pulses.clear();
    return;
  }

  std::vector<FitPulse> keep;
  keep.reserve(pulses.size());
  for (const auto& fp : pulses) {
    const int t0i = (int)std::lround(fp.t0_fit);
    if (IsSingleCandidateInWindowStepwise(starts, t0i, kFitHalfWindow)) {
      keep.push_back(fp);
    }
  }
  pulses.swap(keep);
  n_after += (long)pulses.size();
}

static double SmoothedBinContentStepwise(const TH1D* h, int b, int hw)
{
  if (!h) return 0.0;
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

static int FindPeakBinStepwise(const TH1D* h, double xmin, double xmax, int hw)
{
  if (!h) return -1;
  int b0 = h->GetXaxis()->FindBin(xmin);
  int b1 = h->GetXaxis()->FindBin(xmax);
  if (b1 < b0) std::swap(b0, b1);
  double ymax = -1.0;
  int bbest = -1;
  for (int b = b0; b <= b1; ++b) {
    const double y = SmoothedBinContentStepwise(h, b, hw);
    if (y > ymax) {
      ymax = y;
      bbest = b;
    }
  }
  return bbest;
}

static int FindValleyBinStepwise(const TH1D* h, double xmin, int bpeak, int hw)
{
  if (!h || bpeak < 1) return -1;
  int b0 = h->GetXaxis()->FindBin(xmin);
  int b1 = bpeak - 6;
  if (b1 <= b0) return -1;

  double ymin = 1.0e99;
  int bbest = -1;
  for (int b = b0; b <= b1; ++b) {
    const double y = SmoothedBinContentStepwise(h, b, hw);
    if (y < ymin) {
      ymin = y;
      bbest = b;
    }
  }
  return bbest;
}

static int FindMinBinInRangeStepwise(const TH1D* h, double xmin, double xmax, int hw)
{
  if (!h) return -1;
  int b0 = h->GetXaxis()->FindBin(xmin);
  int b1 = h->GetXaxis()->FindBin(xmax);
  if (b1 < b0) std::swap(b0, b1);
  double ymin = 1.0e99;
  int bbest = -1;
  for (int b = b0; b <= b1; ++b) {
    const double y = SmoothedBinContentStepwise(h, b, hw);
    if (y < ymin) {
      ymin = y;
      bbest = b;
    }
  }
  return bbest;
}

static double DetermineEnergyThresholdFromHist(const TH1D* h_prompt,
                                               const TH1D* h_side,
                                               double scale_side)
{
  // ルール:
  //   低エネルギー側から見て
  //   1) diff = prompt - scale*side が連続で正
  //   2) ratio = prompt/(scale*side) が連続で 1.5 を超える
  //   の両方を満たす立ち上がりを、背景優勢から信号優勢への遷移点とみなす。
  if (!h_prompt || !h_side) return 0.0;
  const int nbin = h_prompt->GetNbinsX();
  int consec_diff = 0;
  int consec_ratio = 0;
  int bdiff = -1;
  int bratio = -1;

  for (int b = 1; b <= nbin; ++b) {
    const double x = h_prompt->GetBinCenter(b);
    if (x < 5000.0) continue;

    const double p = h_prompt->GetBinContent(b);
    const double s = scale_side * h_side->GetBinContent(b);
    const double diff = p - s;
    if (diff > 0.0) {
      consec_diff++;
      if (consec_diff >= 6 && bdiff < 1) bdiff = b - 5;
    } else {
      consec_diff = 0;
    }

    if (p >= 10.0) {
      const double ratio = p / std::max(1.0, s);
      if (ratio > 1.5) {
        consec_ratio++;
        if (consec_ratio >= 6 && bratio < 1) bratio = b - 5;
      } else {
        consec_ratio = 0;
      }
    } else {
      consec_ratio = 0;
    }
  }

  int bfound = -1;
  if (bdiff > 0 && bratio > 0) bfound = std::max(bdiff, bratio);
  else if (bratio > 0) bfound = bratio;
  else if (bdiff > 0) bfound = bdiff;
  if (bfound < 1) return 10000.0;

  double e = h_prompt->GetBinLowEdge(bfound);
  if (e < 10000.0) e = 10000.0;
  // 見た目で使いやすい 500 単位へ丸める
  return 500.0 * std::floor(e / 500.0);
}

static double MichelCoreStepwise(double e, double e0, double w)
{
  if (w <= 0.0) return 0.0;
  const double x = (e - e0) / w;
  if (x <= 0.0 || x >= 1.0) return 0.0;
  return x * x * (3.0 - 2.0 * x);
}

static double MichelAffinePdfStepwise(double* xx, double* p)
{
  // p0: norm, p1: e0, p2: width, p3: offset
  return p[0] * MichelCoreStepwise(xx[0], p[1], p[2]) + p[3];
}

static GausInfo FitPromptPeakStepwise(TH1D* h, double xmin, double xmax)
{
  GausInfo g;
  if (!h || h->GetEntries() < 100.0) return g;
  static int fid = 0;
  TF1 f(Form("f_gaus_stepwise_%d", fid++), "gaus", xmin, xmax);
  f.SetParameters(std::max(1.0, h->GetMaximum()), 0.0, 6.0);
  const int st = h->Fit(&f, "Q0R");
  if (st != 0) return g;
  g.ok = true;
  g.mean = f.GetParameter(1);
  g.sigma = std::abs(f.GetParameter(2));
  g.chi2ndf = (f.GetNDF() > 0) ? (f.GetChisquare() / (double)f.GetNDF()) : 0.0;
  return g;
}

static MichelFitInfo FitMichelExcessStepwise(TH1D* h, double e_min, const char* tag)
{
  MichelFitInfo out;
  if (!h) return out;

  static int fid = 0;
  TF1 f(Form("fMichel_stepwise_%s_%d", (tag ? tag : "x"), fid++),
        MichelAffinePdfStepwise, std::max(0.0, e_min), 2.0e5, 4);
  f.SetParNames("N", "E0", "W", "C");
  f.SetParameters(std::max(10.0, h->GetMaximum()), std::max(0.0, e_min - 5000.0), 120000.0, 0.0);
  f.SetParLimits(0, 0.0, 1.0e8);
  f.SetParLimits(1, 0.0, 1.2e5);
  f.SetParLimits(2, 2.0e4, 3.0e5);
  f.SetParLimits(3, -2.0e4, 2.0e4);

  const int fit_status = h->Fit(&f, "Q0R");
  out.fit_ok = (fit_status == 0);
  out.chi2ndf = (f.GetNDF() > 0) ? (f.GetChisquare() / (double)f.GetNDF()) : 0.0;
  out.n = f.GetParameter(0);
  out.e0 = f.GetParameter(1);
  out.w = f.GetParameter(2);
  out.c = f.GetParameter(3);

  // R^2 指標
  double sse = 0.0;
  double sst = 0.0;
  double ymean = 0.0;
  int nbfit = 0;
  const int b0_fit = h->GetXaxis()->FindBin(std::max(0.0, e_min));
  const int b1_fit = h->GetXaxis()->FindBin(2.0e5);
  for (int b = b0_fit; b <= b1_fit; ++b) {
    ymean += h->GetBinContent(b);
    nbfit++;
  }
  if (nbfit > 0) ymean /= (double)nbfit;
  for (int b = b0_fit; b <= b1_fit; ++b) {
    const double x = h->GetBinCenter(b);
    const double y = h->GetBinContent(b);
    const double yf = f.Eval(x);
    const double dy = y - yf;
    sse += dy * dy;
    const double dt = y - ymean;
    sst += dt * dt;
  }
  out.r2 = (sst > 0.0) ? (1.0 - sse / sst) : 0.0;
  return out;
}

void michel_stepwise_hist_based_run8()
{
  gStyle->SetOptStat(1110);
  gSystem->mkdir("doc/mainexp", true);

  const TString out_pdf = "doc/mainexp/michel_stepwise_hist_based_run8.pdf";

  // ------------------------------------------------------------
  // Pass 1: PSヒスト + NaI |dt12| ヒスト（閾値決定用）
  // ------------------------------------------------------------
  TH1D* hPSAmpA = new TH1D("hPSAmpA_stepwise", "PS_A pulse amplitude;baseline-ADC [counts];Pulses",
                           kPSAmpBins, kPSAmpMin, kPSAmpMax);
  TH1D* hPSAmpB = new TH1D("hPSAmpB_stepwise", "PS_B pulse amplitude;baseline-ADC [counts];Pulses",
                           kPSAmpBins, kPSAmpMin, kPSAmpMax);
  TH1D* hAbsDt12A = new TH1D("hAbsDt12A_stepwise", "A nearest |#Delta t_{12}|;|#Delta t_{12}| [samples];Pairs",
                             200, 0, 100);
  TH1D* hAbsDt12B = new TH1D("hAbsDt12B_stepwise", "B nearest |#Delta t_{12}|;|#Delta t_{12}| [samples];Pairs",
                             200, 0, 100);
  TH1D* hDt12A = new TH1D("hDt12A_stepwise", "A nearest #Delta t_{12};#Delta t_{12} [samples];Pairs",
                          240, -120, 120);
  TH1D* hDt12B = new TH1D("hDt12B_stepwise", "B nearest #Delta t_{12};#Delta t_{12} [samples];Pairs",
                          240, -120, 120);

  std::vector<double> absdtA;
  std::vector<double> absdtB;
  absdtA.reserve(50000);
  absdtB.reserve(20000);
  long nai_hits_before_clean_pass1 = 0;
  long nai_hits_after_clean_pass1 = 0;

  long n_events_pass1 = 0;
  {
    TString fps[2], fna[4];
    for (int i = 0; i < 2; ++i) fps[i] = Form("%s/%s", kInputDir, kPSFiles[i]);
    for (int i = 0; i < 4; ++i) fna[i] = Form("%s/%s", kInputDir, kNaIFiles[i]);

    std::ifstream in_ps[2];
    std::ifstream in_na[4];
    for (int i = 0; i < 2; ++i) {
      in_ps[i].open(fps[i].Data());
      if (!in_ps[i]) {
        std::cerr << "[stepwise] ERROR: cannot open " << fps[i] << std::endl;
        return;
      }
    }
    for (int i = 0; i < 4; ++i) {
      in_na[i].open(fna[i].Data());
      if (!in_na[i]) {
        std::cerr << "[stepwise] ERROR: cannot open " << fna[i] << std::endl;
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

      n_events_pass1++;

      // PSヒスト
      for (int ip = 0; ip < 2; ++ip) {
        const double b = ModeBaselineInEvent(buf_ps[ip]);
        std::vector<double> s;
        BuildS(buf_ps[ip], b, s);
        const std::vector<Pulse> ps = FindPulsesSimple(s, kPSSeedThr, kPSEndThr, kPSMinSep, kPSMinWidth);
        for (const auto& q : ps) {
          if (q.amp <= 0.0) continue;
          if (ip == 0) hPSAmpA->Fill(q.amp);
          else hPSAmpB->Fill(q.amp);
        }
      }

      // NaI分解
      std::vector<FitPulse> fp[4];
      for (int ch = 0; ch < 4; ++ch) {
        const double b = ModeBaselineInEvent(buf_na[ch]);
        std::vector<double> s;
        BuildS(buf_na[ch], b, s);
        fp[ch] = FitPulsesWithDecomposition(s, kTrFix[ch], kDFix[ch]);
        ApplyNaICleanSelectionStepwise(s, fp[ch],
                                       nai_hits_before_clean_pass1,
                                       nai_hits_after_clean_pass1);
      }

      // A nearest dt12
      if (!fp[0].empty() && !fp[1].empty()) {
        for (const auto& p1 : fp[0]) {
          const int j = FindNearestPulseIndexStepwise(fp[1], p1.t0_fit);
          if (j < 0) continue;
          const double dt12 = fp[1][j].t0_fit - p1.t0_fit;
          const double adt = std::abs(dt12);
          hDt12A->Fill(dt12);
          hAbsDt12A->Fill(adt);
          absdtA.push_back(adt);
        }
      }
      // B nearest dt12
      if (!fp[2].empty() && !fp[3].empty()) {
        for (const auto& p1 : fp[2]) {
          const int j = FindNearestPulseIndexStepwise(fp[3], p1.t0_fit);
          if (j < 0) continue;
          const double dt12 = fp[3][j].t0_fit - p1.t0_fit;
          const double adt = std::abs(dt12);
          hDt12B->Fill(dt12);
          hAbsDt12B->Fill(adt);
          absdtB.push_back(adt);
        }
      }

      for (int i = 0; i < 2; ++i) buf_ps[i].clear();
      for (int i = 0; i < 4; ++i) buf_na[i].clear();
    }
  }

  // A) PS閾値決定
  // study_michel_ps_tagged.C で実績のある valley 抽出をそのまま使う。
  int bpkA = -1, bvlA = -1;
  int bpkB = -1, bvlB = -1;
  double ps_thr_A = FindPSValleyThreshold(hPSAmpA, bpkA, bvlA);
  double ps_thr_B = FindPSValleyThreshold(hPSAmpB, bpkB, bvlB);

  // B) NaI pair dt12 cut 決定
  //    signed #Delta t12 の中心ピーク幅(4sigma)を採用する。
  //    A 側ピークを取りこぼしにくいよう、3sigma より広く取る。
  std::sort(absdtA.begin(), absdtA.end());
  std::sort(absdtB.begin(), absdtB.end());
  TH1D* hEffDt12A = new TH1D("hEffDt12A_stepwise", "A cut efficiency;|#Delta t_{12}| cut [samples];Accepted fraction",
                             80, 0, 80);
  TH1D* hEffDt12B = new TH1D("hEffDt12B_stepwise", "B cut efficiency;|#Delta t_{12}| cut [samples];Accepted fraction",
                             80, 0, 80);
  for (int b = 1; b <= hEffDt12A->GetNbinsX(); ++b) {
    const double x = hEffDt12A->GetBinCenter(b);
    hEffDt12A->SetBinContent(b, FractionLeqFromSortedStepwise(absdtA, x));
    hEffDt12B->SetBinContent(b, FractionLeqFromSortedStepwise(absdtB, x));
  }
  GausInfo gPairA = FitPromptPeakStepwise(hDt12A, -25.0, 25.0);
  GausInfo gPairB = FitPromptPeakStepwise(hDt12B, -25.0, 25.0);
  const int bpkPairA = FindPeakBinStepwise(hAbsDt12A, 0.0, 12.0, 2);
  const int bpkPairB = FindPeakBinStepwise(hAbsDt12B, 0.0, 12.0, 2);
  const double xpkPairA = (bpkPairA > 0) ? hAbsDt12A->GetBinCenter(bpkPairA) : 6.0;
  const double xpkPairB = (bpkPairB > 0) ? hAbsDt12B->GetBinCenter(bpkPairB) : 6.0;
  const double cutA_4sigma = (gPairA.ok && gPairA.sigma > 1.0 && gPairA.sigma < 30.0) ? (4.0 * gPairA.sigma) : 0.0;
  const double cutB_4sigma = (gPairB.ok && gPairB.sigma > 1.0 && gPairB.sigma < 30.0) ? (4.0 * gPairB.sigma) : 0.0;
  const double cut_by_sigma = std::max(cutA_4sigma, cutB_4sigma);
  const double q90A = QuantileFromSortedStepwise(absdtA, 0.90);
  const double q90B = QuantileFromSortedStepwise(absdtB, 0.90);
  const double q95A = QuantileFromSortedStepwise(absdtA, 0.95);
  const double q95B = QuantileFromSortedStepwise(absdtB, 0.95);
  int pair_cut = 30;
  if (cut_by_sigma > 0.0) pair_cut = (int)(5.0 * std::ceil(cut_by_sigma / 5.0));
  else pair_cut = (int)(5.0 * std::ceil(std::max(q90A, q90B) / 5.0));
  if (pair_cut < 20) pair_cut = 20;
  if (pair_cut > 50) pair_cut = 50;

  // ------------------------------------------------------------
  // Pass 2: A/B閾値を使って dt と E 候補を作る
  // ------------------------------------------------------------
  TH1D* hDtA = new TH1D("hDtA_stepwise", "A side dt=t_{mod}-t_{PS};dt [samples];Candidates", 600, -300, 300);
  TH1D* hDtB = new TH1D("hDtB_stepwise", "B side dt=t_{mod}-t_{PS};dt [samples];Candidates", 600, -300, 300);
  TH1D* hAbsDt = new TH1D("hAbsDt_stepwise", "Combined |dt|;|dt| [samples];Candidates", 300, 0, 300);

  std::vector<SideCandidate> candA;
  std::vector<SideCandidate> candB;
  candA.reserve(50000);
  candB.reserve(50000);
  long nai_hits_before_clean_pass2 = 0;
  long nai_hits_after_clean_pass2 = 0;

  long n_events_pass2 = 0;
  {
    TString fps[2], fna[4];
    for (int i = 0; i < 2; ++i) fps[i] = Form("%s/%s", kInputDir, kPSFiles[i]);
    for (int i = 0; i < 4; ++i) fna[i] = Form("%s/%s", kInputDir, kNaIFiles[i]);

    std::ifstream in_ps[2];
    std::ifstream in_na[4];
    for (int i = 0; i < 2; ++i) {
      in_ps[i].open(fps[i].Data());
      if (!in_ps[i]) {
        std::cerr << "[stepwise] ERROR: cannot open " << fps[i] << std::endl;
        return;
      }
    }
    for (int i = 0; i < 4; ++i) {
      in_na[i].open(fna[i].Data());
      if (!in_na[i]) {
        std::cerr << "[stepwise] ERROR: cannot open " << fna[i] << std::endl;
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
      n_events_pass2++;

      // PS候補
      std::vector<double> tps_sel[2];
      std::vector<double> aps_sel[2];
      for (int ip = 0; ip < 2; ++ip) {
        const double b = ModeBaselineInEvent(buf_ps[ip]);
        std::vector<double> s;
        BuildS(buf_ps[ip], b, s);
        const std::vector<Pulse> ps = FindPulsesSimple(s, kPSSeedThr, kPSEndThr, kPSMinSep, kPSMinWidth);
        for (const auto& q : ps) {
          const double thr = (ip == 0) ? ps_thr_A : ps_thr_B;
          if (q.amp < thr) continue;
          tps_sel[ip].push_back((double)q.start);
          aps_sel[ip].push_back(q.amp);
        }
      }

      // NaI分解 + module
      std::vector<FitPulse> fp[4];
      for (int ch = 0; ch < 4; ++ch) {
        const double b = ModeBaselineInEvent(buf_na[ch]);
        std::vector<double> s;
        BuildS(buf_na[ch], b, s);
        fp[ch] = FitPulsesWithDecomposition(s, kTrFix[ch], kDFix[ch]);
        ApplyNaICleanSelectionStepwise(s, fp[ch],
                                       nai_hits_before_clean_pass2,
                                       nai_hits_after_clean_pass2);
      }
      const std::vector<ModulePulse> modA = BuildModulePairs(fp[0], fp[1], pair_cut);
      const std::vector<ModulePulse> modB = BuildModulePairs(fp[2], fp[3], pair_cut);

      SideCandidate cA;
      SideCandidate cB;

      // A side
      if (!tps_sel[0].empty() && !modA.empty()) {
        int iear = 0;
        for (int i = 1; i < (int)tps_sel[0].size(); ++i) {
          if (tps_sel[0][i] < tps_sel[0][iear]) iear = i;
        }
        const double tps = tps_sel[0][iear];
        const int im = FindNearestModuleIndex(modA, tps);
        if (im >= 0) {
          cA.ok = true;
          cA.dt = modA[im].t - tps;
          cA.e_sum = modA[im].area_sum;
          cA.ps_amp = aps_sel[0][iear];
          hDtA->Fill(cA.dt);
          hAbsDt->Fill(std::abs(cA.dt));
        }
      }
      // B side
      if (!tps_sel[1].empty() && !modB.empty()) {
        int iear = 0;
        for (int i = 1; i < (int)tps_sel[1].size(); ++i) {
          if (tps_sel[1][i] < tps_sel[1][iear]) iear = i;
        }
        const double tps = tps_sel[1][iear];
        const int im = FindNearestModuleIndex(modB, tps);
        if (im >= 0) {
          cB.ok = true;
          cB.dt = modB[im].t - tps;
          cB.e_sum = modB[im].area_sum;
          cB.ps_amp = aps_sel[1][iear];
          hDtB->Fill(cB.dt);
          hAbsDt->Fill(std::abs(cB.dt));
        }
      }

      candA.push_back(cA);
      candB.push_back(cB);

      for (int i = 0; i < 2; ++i) buf_ps[i].clear();
      for (int i = 0; i < 4; ++i) buf_na[i].clear();
    }
  }

  // C) prompt時間窓決定（|dt|ヒストの谷から）
  const int bpkDtA = FindPeakBinStepwise(hDtA, -40.0, 40.0, 2);
  const int bpkDtB = FindPeakBinStepwise(hDtB, -40.0, 40.0, 2);
  const double dtpkA = (bpkDtA > 0) ? hDtA->GetBinCenter(bpkDtA) : 0.0;
  const double dtpkB = (bpkDtB > 0) ? hDtB->GetBinCenter(bpkDtB) : 0.0;
  GausInfo gA = FitPromptPeakStepwise(hDtA, dtpkA - 20.0, dtpkA + 20.0);
  GausInfo gB = FitPromptPeakStepwise(hDtB, dtpkB - 20.0, dtpkB + 20.0);
  const int bpkAbsDt = FindPeakBinStepwise(hAbsDt, 0.0, 12.0, 2);
  const double xpkAbsDt = (bpkAbsDt > 0) ? hAbsDt->GetBinCenter(bpkAbsDt) : 4.0;
  const int bvlAbsDt = FindMinBinInRangeStepwise(hAbsDt, xpkAbsDt + 6.0, 90.0, 4);
  std::vector<double> sigmas_dt;
  if (gA.ok && gA.sigma > 2.0 && gA.sigma < 30.0) sigmas_dt.push_back(gA.sigma);
  if (gB.ok && gB.sigma > 2.0 && gB.sigma < 30.0) sigmas_dt.push_back(gB.sigma);
  double sigma_ref = 12.0;
  if (!sigmas_dt.empty()) sigma_ref = *std::max_element(sigmas_dt.begin(), sigmas_dt.end());
  int prompt_halfwin = (int)std::ceil(1.6 * sigma_ref);
  if (bvlAbsDt > 0) {
    const int valley_w = (int)std::lround(hAbsDt->GetBinCenter(bvlAbsDt));
    if (valley_w < prompt_halfwin) prompt_halfwin = valley_w;
  }
  if (prompt_halfwin < 12) prompt_halfwin = 12;
  if (prompt_halfwin > 30) prompt_halfwin = 30;

  // D) エネルギー閾値決定
  TH1D* hEPromptPre = new TH1D("hEPromptPre_stepwise", "E in prompt window;E_{sum} [ADC*samples];Counts",
                               2500, 0, 2.0e5);
  TH1D* hESidePre = new TH1D("hESidePre_stepwise", "E in sideband;E_{sum} [ADC*samples];Counts",
                             2500, 0, 2.0e5);
  const double sideband_lo = 80.0;
  const double sideband_hi = 160.0;
  hEPromptPre->Sumw2();
  hESidePre->Sumw2();

  for (size_t i = 0; i < candA.size(); ++i) {
    const SideCandidate& a = candA[i];
    if (a.ok) {
      const double adt = std::abs(a.dt);
      if (adt <= (double)prompt_halfwin) hEPromptPre->Fill(a.e_sum);
      if (sideband_lo <= adt && adt <= sideband_hi) hESidePre->Fill(a.e_sum);
    }
    const SideCandidate& b = candB[i];
    if (b.ok) {
      const double adt = std::abs(b.dt);
      if (adt <= (double)prompt_halfwin) hEPromptPre->Fill(b.e_sum);
      if (sideband_lo <= adt && adt <= sideband_hi) hESidePre->Fill(b.e_sum);
    }
  }

  TH1D* hEExcessPre = (TH1D*)hEPromptPre->Clone("hEExcessPre_stepwise");
  const double sf_side = (double)prompt_halfwin / (sideband_hi - sideband_lo); // width(prompt)/width(sideband)
  hEExcessPre->Add(hESidePre, -sf_side);
  double e_min = DetermineEnergyThresholdFromHist(hEPromptPre, hESidePre, sf_side);
  const int bpkLowE = FindPeakBinStepwise(hEExcessPre, 1000.0, 20000.0, 8);
  const int bpkHighE = FindPeakBinStepwise(hEExcessPre, 15000.0, 140000.0, 10);
  double e_valley = 0.0;
  if (bpkLowE > 0 && bpkHighE > bpkLowE + 6) {
    const double x0 = hEExcessPre->GetBinCenter(bpkLowE) + 2000.0;
    const double x1 = hEExcessPre->GetBinCenter(bpkHighE) - 2000.0;
    const int bvlE = FindMinBinInRangeStepwise(hEExcessPre, x0, x1, 10);
    if (bvlE > 0) e_valley = hEExcessPre->GetBinLowEdge(bvlE);
  }
  if (e_valley > 0.0) e_min = std::max(e_min, 500.0 * std::floor(e_valley / 500.0));
  if (e_min < 12000.0) e_min = 12000.0;

  // E) 最終ヒスト + Michel重ね合わせ
  TH1D* hDtFinal = new TH1D("hDtFinal_stepwise", "dt after final cuts;dt [samples];Events", 600, -300, 300);
  TH1D* hAbsDtFinal = new TH1D("hAbsDtFinal_stepwise", "|dt| after final cuts;|dt| [samples];Events", 300, 0, 300);
  TH1D* hEPrompt = new TH1D("hEPrompt_stepwise", "Prompt E (final cuts);E_{sum} [ADC*samples];Counts",
                            2500, 0, 2.0e5);
  TH1D* hESide = new TH1D("hESide_stepwise", "Sideband E (final cuts);E_{sum} [ADC*samples];Counts",
                          2500, 0, 2.0e5);
  TH1D* hEPromptA = new TH1D("hEPromptA_stepwise", "Prompt E A-side;E_{sum} [ADC*samples];Counts",
                             2500, 0, 2.0e5);
  TH1D* hESideA = new TH1D("hESideA_stepwise", "Sideband E A-side;E_{sum} [ADC*samples];Counts",
                           2500, 0, 2.0e5);
  TH1D* hEPromptB = new TH1D("hEPromptB_stepwise", "Prompt E B-side;E_{sum} [ADC*samples];Counts",
                             2500, 0, 2.0e5);
  TH1D* hESideB = new TH1D("hESideB_stepwise", "Sideband E B-side;E_{sum} [ADC*samples];Counts",
                           2500, 0, 2.0e5);
  hEPrompt->Sumw2();
  hESide->Sumw2();
  hEPromptA->Sumw2();
  hESideA->Sumw2();
  hEPromptB->Sumw2();
  hESideB->Sumw2();

  for (size_t i = 0; i < candA.size(); ++i) {
    const SideCandidate& a = candA[i];
    if (a.ok && a.e_sum >= e_min) {
      hDtFinal->Fill(a.dt);
      hAbsDtFinal->Fill(std::abs(a.dt));
      const double adt = std::abs(a.dt);
      if (adt <= (double)prompt_halfwin) {
        hEPrompt->Fill(a.e_sum);
        hEPromptA->Fill(a.e_sum);
      }
      if (sideband_lo <= adt && adt <= sideband_hi) {
        hESide->Fill(a.e_sum);
        hESideA->Fill(a.e_sum);
      }
    }
    const SideCandidate& b = candB[i];
    if (b.ok && b.e_sum >= e_min) {
      hDtFinal->Fill(b.dt);
      hAbsDtFinal->Fill(std::abs(b.dt));
      const double adt = std::abs(b.dt);
      if (adt <= (double)prompt_halfwin) {
        hEPrompt->Fill(b.e_sum);
        hEPromptB->Fill(b.e_sum);
      }
      if (sideband_lo <= adt && adt <= sideband_hi) {
        hESide->Fill(b.e_sum);
        hESideB->Fill(b.e_sum);
      }
    }
  }

  TH1D* hEExcess = (TH1D*)hEPrompt->Clone("hEExcess_stepwise");
  hEExcess->SetTitle("Prompt - scaled sideband;E_{sum} [ADC*samples];Excess");
  hEExcess->Add(hESide, -sf_side);
  TH1D* hEExcessA = (TH1D*)hEPromptA->Clone("hEExcessA_stepwise");
  hEExcessA->SetTitle("A-side: Prompt - scaled sideband;E_{sum} [ADC*samples];Excess");
  hEExcessA->Add(hESideA, -sf_side);
  TH1D* hEExcessB = (TH1D*)hEPromptB->Clone("hEExcessB_stepwise");
  hEExcessB->SetTitle("B-side: Prompt - scaled sideband;E_{sum} [ADC*samples];Excess");
  hEExcessB->Add(hESideB, -sf_side);

  const MichelFitInfo fit_all = FitMichelExcessStepwise(hEExcess, e_min, "all");
  const MichelFitInfo fit_a = FitMichelExcessStepwise(hEExcessA, e_min, "A");
  const MichelFitInfo fit_b = FitMichelExcessStepwise(hEExcessB, e_min, "B");

  // --------------------------
  // 出力 (数値ログ)
  // --------------------------
  std::cout << "\n[stepwise] === Decision Summary ===\n";
  std::cout << "[stepwise] events_pass1=" << n_events_pass1
            << " events_pass2=" << n_events_pass2 << "\n";
  std::cout << "[stepwise] NaI clean (single candidate in fit window: +/-" << kFitHalfWindow
            << "): pass1 kept=" << nai_hits_after_clean_pass1
            << "/" << nai_hits_before_clean_pass1
            << " pass2 kept=" << nai_hits_after_clean_pass2
            << "/" << nai_hits_before_clean_pass2 << "\n";
  std::cout << "[stepwise] A_PS_peak=" << ((bpkA > 0) ? hPSAmpA->GetBinCenter(bpkA) : -1.0)
            << " A_PS_valley=" << ((bvlA > 0) ? hPSAmpA->GetBinCenter(bvlA) : -1.0)
            << " -> PS_thr_A=" << ps_thr_A << "\n";
  std::cout << "[stepwise] B_PS_peak=" << ((bpkB > 0) ? hPSAmpB->GetBinCenter(bpkB) : -1.0)
            << " B_PS_valley=" << ((bvlB > 0) ? hPSAmpB->GetBinCenter(bvlB) : -1.0)
            << " -> PS_thr_B=" << ps_thr_B << "\n";
  std::cout << "[stepwise] pair_peakA=" << xpkPairA
            << " sigmaPairA=" << (gPairA.ok ? gPairA.sigma : -1.0)
            << " pair_peakB=" << xpkPairB
            << " sigmaPairB=" << (gPairB.ok ? gPairB.sigma : -1.0)
            << " cutA_4sigma=" << cutA_4sigma
            << " cutB_4sigma=" << cutB_4sigma
            << " sigma_based_cut=" << cut_by_sigma << "\n";
  std::cout << "[stepwise] q90_absdt12_A=" << q90A
            << " q90_absdt12_B=" << q90B << "\n";
  std::cout << "[stepwise] q95_absdt12_A=" << q95A
            << " q95_absdt12_B=" << q95B
            << " -> pair_cut=" << pair_cut << "\n";
  std::cout << "[stepwise] pair_eff_A(|dt12|<cut)=" << FractionLeqFromSortedStepwise(absdtA, (double)pair_cut)
            << " pair_eff_B(|dt12|<cut)=" << FractionLeqFromSortedStepwise(absdtB, (double)pair_cut) << "\n";
  std::cout << "[stepwise] dt_peak_A=" << dtpkA
            << " sigmaA=" << (gA.ok ? gA.sigma : -1.0)
            << " dt_peak_B=" << dtpkB
            << " sigmaB=" << (gB.ok ? gB.sigma : -1.0)
            << " absdt_peak=" << xpkAbsDt
            << " absdt_valley=" << ((bvlAbsDt > 0) ? hAbsDt->GetBinCenter(bvlAbsDt) : -1.0)
            << " -> prompt_halfwin=" << prompt_halfwin << "\n";
  std::cout << "[stepwise] sideband_scale=" << sf_side
            << " E_valley=" << e_valley
            << " -> E_min=" << e_min << "\n";
  std::cout << "[stepwise] final_prompt_entries=" << hEPrompt->Integral()
            << " final_side_entries=" << hESide->Integral()
            << " excess_integral=" << hEExcess->Integral() << "\n";
  std::cout << "[stepwise] michel_fit_all_ok=" << (fit_all.fit_ok ? 1 : 0)
            << " chi2_ndf=" << fit_all.chi2ndf
            << " R2=" << fit_all.r2
            << " E0_fit=" << fit_all.e0
            << " W_fit=" << fit_all.w
            << " C_fit=" << fit_all.c << "\n";
  std::cout << "[stepwise] michel_fit_A_ok=" << (fit_a.fit_ok ? 1 : 0)
            << " chi2_ndf_A=" << fit_a.chi2ndf
            << " R2_A=" << fit_a.r2
            << " E0_A=" << fit_a.e0
            << " W_A=" << fit_a.w
            << " C_A=" << fit_a.c << "\n";
  std::cout << "[stepwise] michel_fit_B_ok=" << (fit_b.fit_ok ? 1 : 0)
            << " chi2_ndf_B=" << fit_b.chi2ndf
            << " R2_B=" << fit_b.r2
            << " E0_B=" << fit_b.e0
            << " W_B=" << fit_b.w
            << " C_B=" << fit_b.c << "\n\n";

  // --------------------------
  // PDF描画
  // --------------------------
  TCanvas* c = new TCanvas("c_stepwise", "Michel stepwise", 1400, 1000);
  c->Print(out_pdf + "[");

  c->Clear();
  TPaveText* p = new TPaveText(0.05, 0.05, 0.95, 0.95, "NDC");
  p->SetFillStyle(0);
  p->SetBorderSize(0);
  p->SetTextAlign(12);
  p->SetTextSize(0.028);
  p->AddText("=== Stepwise histogram-based decision for Michel selection (run8) ===");
  p->AddText(Form("Input: %s", kInputDir));
  p->AddText(Form("NaI clean condition: single start candidate within +/- %d samples around fit t0", kFitHalfWindow));
  p->AddText(Form("NaI clean kept: pass1 %ld/%ld, pass2 %ld/%ld",
                  nai_hits_after_clean_pass1, nai_hits_before_clean_pass1,
                  nai_hits_after_clean_pass2, nai_hits_before_clean_pass2));
  p->AddText("Decision flow:");
  p->AddText(Form(" A) PS hist -> PS_thr_A=%.2f, PS_thr_B=%.2f", ps_thr_A, ps_thr_B));
  p->AddText(Form(" B) |dt12| hist -> pair cut |dt12|<=%d", pair_cut));
  p->AddText(Form(" C) dt hist -> prompt window |dt|<=%d", prompt_halfwin));
  p->AddText(Form(" D) sideband subtraction: H_excess(E)=H_prompt(E)-s*H_side(E), s=%.3f", sf_side));
  p->AddText(Form("    threshold from excess valley -> E_sum>=%.0f", e_min));
  p->AddText(" E) Michel fit: f(E)=N*x^{2}(3-2x)+C, x=(E-E0)/W");
  p->AddText(" ");
  p->AddText(Form("Michel fit (all): chi2/ndf=%.3f, R2=%.3f, E0=%.1f, W=%.1f",
                  fit_all.chi2ndf, fit_all.r2, fit_all.e0, fit_all.w));
  p->AddText(Form("Michel fit (A/B): chi2A=%.3f, chi2B=%.3f", fit_a.chi2ndf, fit_b.chi2ndf));
  p->Draw();
  c->Print(out_pdf);
  delete p;

  // page: PS thresholds
  c->Clear();
  c->Divide(2, 1);
  c->cd(1);
  gPad->SetGridx(); gPad->SetGridy(); gPad->SetLogy();
  hPSAmpA->SetLineColor(kBlue + 1);
  hPSAmpA->SetLineWidth(2);
  hPSAmpA->Draw("hist");
  {
    const double ymax = std::max(1.0, hPSAmpA->GetMaximum());
    if (bpkA > 0) {
      TLine* lpk = new TLine(hPSAmpA->GetBinCenter(bpkA), 0.1, hPSAmpA->GetBinCenter(bpkA), ymax);
      lpk->SetLineColor(kMagenta + 2);
      lpk->SetLineStyle(3);
      lpk->Draw();
    }
    TLine* lthr = new TLine(ps_thr_A, 0.1, ps_thr_A, ymax);
    lthr->SetLineColor(kRed + 1);
    lthr->SetLineStyle(2);
    lthr->SetLineWidth(2);
    lthr->Draw();
  }

  c->cd(2);
  gPad->SetGridx(); gPad->SetGridy(); gPad->SetLogy();
  hPSAmpB->SetLineColor(kBlue + 1);
  hPSAmpB->SetLineWidth(2);
  hPSAmpB->Draw("hist");
  {
    const double ymax = std::max(1.0, hPSAmpB->GetMaximum());
    if (bpkB > 0) {
      TLine* lpk = new TLine(hPSAmpB->GetBinCenter(bpkB), 0.1, hPSAmpB->GetBinCenter(bpkB), ymax);
      lpk->SetLineColor(kMagenta + 2);
      lpk->SetLineStyle(3);
      lpk->Draw();
    }
    TLine* lthr = new TLine(ps_thr_B, 0.1, ps_thr_B, ymax);
    lthr->SetLineColor(kRed + 1);
    lthr->SetLineStyle(2);
    lthr->SetLineWidth(2);
    lthr->Draw();
  }
  c->Print(out_pdf);

  // page: pair dt12 signed (with Gaussian fit window and cut)
  c->Clear();
  c->Divide(2, 1);
  c->cd(1);
  gPad->SetGridx(); gPad->SetGridy(); gPad->SetLogy();
  hDt12A->SetLineColor(kBlue + 1);
  hDt12A->SetLineWidth(2);
  hDt12A->Draw("hist");
  {
    const double ymax = std::max(1.0, hDt12A->GetMaximum());
    TLine* l1 = new TLine(-(double)pair_cut, 0.5, -(double)pair_cut, ymax);
    TLine* l2 = new TLine(+(double)pair_cut, 0.5, +(double)pair_cut, ymax);
    l1->SetLineColor(kRed + 1); l2->SetLineColor(kRed + 1);
    l1->SetLineStyle(2); l2->SetLineStyle(2);
    l1->SetLineWidth(2); l2->SetLineWidth(2);
    l1->Draw(); l2->Draw();
    if (gPairA.ok) {
      TF1 fga("fga_pairA_draw", "gaus", -25.0, 25.0);
      fga.SetParameters(std::max(1.0, hDt12A->GetMaximum()), gPairA.mean, gPairA.sigma);
      fga.SetLineColor(kGreen + 2);
      fga.SetLineStyle(3);
      fga.SetLineWidth(2);
      fga.Draw("same");
    }
  }

  c->cd(2);
  gPad->SetGridx(); gPad->SetGridy(); gPad->SetLogy();
  hDt12B->SetLineColor(kBlue + 1);
  hDt12B->SetLineWidth(2);
  hDt12B->Draw("hist");
  {
    const double ymax = std::max(1.0, hDt12B->GetMaximum());
    TLine* l1 = new TLine(-(double)pair_cut, 0.5, -(double)pair_cut, ymax);
    TLine* l2 = new TLine(+(double)pair_cut, 0.5, +(double)pair_cut, ymax);
    l1->SetLineColor(kRed + 1); l2->SetLineColor(kRed + 1);
    l1->SetLineStyle(2); l2->SetLineStyle(2);
    l1->SetLineWidth(2); l2->SetLineWidth(2);
    l1->Draw(); l2->Draw();
    if (gPairB.ok) {
      TF1 fgb("fgb_pairB_draw", "gaus", -25.0, 25.0);
      fgb.SetParameters(std::max(1.0, hDt12B->GetMaximum()), gPairB.mean, gPairB.sigma);
      fgb.SetLineColor(kGreen + 2);
      fgb.SetLineStyle(3);
      fgb.SetLineWidth(2);
      fgb.Draw("same");
    }
  }
  c->Print(out_pdf);

  // page: pair dt12 abs and acceptance curve
  c->Clear();
  c->Divide(2, 1);
  c->cd(1);
  gPad->SetGridx(); gPad->SetGridy(); gPad->SetLogy();
  hAbsDt12A->SetLineColor(kBlue + 1);
  hAbsDt12B->SetLineColor(kOrange + 7);
  hAbsDt12A->SetLineWidth(2);
  hAbsDt12B->SetLineWidth(2);
  hAbsDt12A->Draw("hist");
  hAbsDt12B->Draw("hist same");
  {
    const double ymax = std::max(hAbsDt12A->GetMaximum(), hAbsDt12B->GetMaximum());
    TLine* lc = new TLine((double)pair_cut, 0.5, (double)pair_cut, std::max(1.0, ymax));
    lc->SetLineColor(kRed + 1);
    lc->SetLineStyle(2);
    lc->SetLineWidth(2);
    lc->Draw();
    TLegend* lg = new TLegend(0.50, 0.74, 0.92, 0.90);
    lg->SetBorderSize(0);
    lg->SetFillStyle(0);
    lg->AddEntry(hAbsDt12A, "A: |dt12|", "l");
    lg->AddEntry(hAbsDt12B, "B: |dt12|", "l");
    lg->AddEntry(lc, Form("cut=%d", pair_cut), "l");
    lg->Draw();
  }

  c->cd(2);
  gPad->SetGridx(); gPad->SetGridy();
  hEffDt12A->SetLineColor(kBlue + 1);
  hEffDt12B->SetLineColor(kOrange + 7);
  hEffDt12A->SetLineWidth(2);
  hEffDt12B->SetLineWidth(2);
  hEffDt12A->SetMinimum(0.0);
  hEffDt12A->SetMaximum(1.05);
  hEffDt12A->Draw("hist");
  hEffDt12B->Draw("hist same");
  {
    TLine* lc = new TLine((double)pair_cut, 0.0, (double)pair_cut, 1.05);
    lc->SetLineColor(kRed + 1);
    lc->SetLineStyle(2);
    lc->SetLineWidth(2);
    lc->Draw();
    TLegend* lg = new TLegend(0.48, 0.14, 0.92, 0.34);
    lg->SetBorderSize(0);
    lg->SetFillStyle(0);
    lg->AddEntry(hEffDt12A, "A: frac(|dt12|<cut)", "l");
    lg->AddEntry(hEffDt12B, "B: frac(|dt12|<cut)", "l");
    lg->Draw();
  }
  c->Print(out_pdf);

  // page: combined |dt| for prompt-window decision
  c->Clear();
  gPad->SetGridx(); gPad->SetGridy();
  hAbsDt->SetLineColor(kBlue + 1);
  hAbsDt->SetLineWidth(2);
  hAbsDt->Draw("hist");
  {
    const double ymax = std::max(1.0, hAbsDt->GetMaximum());
    if (bpkAbsDt > 0) {
      TLine* lp = new TLine(hAbsDt->GetBinCenter(bpkAbsDt), 0.0, hAbsDt->GetBinCenter(bpkAbsDt), ymax);
      lp->SetLineColor(kMagenta + 2);
      lp->SetLineStyle(3);
      lp->Draw();
    }
    TLine* lv = new TLine((double)prompt_halfwin, 0.0, (double)prompt_halfwin, ymax);
    lv->SetLineColor(kRed + 1);
    lv->SetLineStyle(2);
    lv->SetLineWidth(2);
    lv->Draw();
  }
  c->Print(out_pdf);

  // page: dt window
  c->Clear();
  c->Divide(2, 1);
  c->cd(1);
  gPad->SetGridx(); gPad->SetGridy();
  hDtA->SetLineColor(kBlue + 1);
  hDtA->SetLineWidth(2);
  hDtA->Draw("hist");
  {
    const double ymax = std::max(1.0, hDtA->GetMaximum());
    TLine* l1 = new TLine(-prompt_halfwin, 0.0, -prompt_halfwin, ymax);
    TLine* l2 = new TLine(+prompt_halfwin, 0.0, +prompt_halfwin, ymax);
    l1->SetLineColor(kRed + 1); l2->SetLineColor(kRed + 1);
    l1->SetLineStyle(2); l2->SetLineStyle(2);
    l1->Draw(); l2->Draw();
  }

  c->cd(2);
  gPad->SetGridx(); gPad->SetGridy();
  hDtB->SetLineColor(kBlue + 1);
  hDtB->SetLineWidth(2);
  hDtB->Draw("hist");
  {
    const double ymax = std::max(1.0, hDtB->GetMaximum());
    TLine* l1 = new TLine(-prompt_halfwin, 0.0, -prompt_halfwin, ymax);
    TLine* l2 = new TLine(+prompt_halfwin, 0.0, +prompt_halfwin, ymax);
    l1->SetLineColor(kRed + 1); l2->SetLineColor(kRed + 1);
    l1->SetLineStyle(2); l2->SetLineStyle(2);
    l1->Draw(); l2->Draw();
  }
  c->Print(out_pdf);

  // page: |dt| window definition for sideband subtraction
  c->Clear();
  gPad->SetGridx(); gPad->SetGridy();
  hAbsDt->SetLineColor(kBlue + 1);
  hAbsDt->SetLineWidth(2);
  hAbsDt->Draw("hist");
  {
    const double ymax = std::max(1.0, hAbsDt->GetMaximum());
    TLine* lp = new TLine((double)prompt_halfwin, 0.0, (double)prompt_halfwin, ymax);
    lp->SetLineColor(kRed + 1);
    lp->SetLineStyle(2);
    lp->SetLineWidth(2);
    lp->Draw();
    TLine* ls1 = new TLine(sideband_lo, 0.0, sideband_lo, ymax);
    TLine* ls2 = new TLine(sideband_hi, 0.0, sideband_hi, ymax);
    ls1->SetLineColor(kGreen + 2); ls2->SetLineColor(kGreen + 2);
    ls1->SetLineStyle(3); ls2->SetLineStyle(3);
    ls1->SetLineWidth(2); ls2->SetLineWidth(2);
    ls1->Draw(); ls2->Draw();
    TLatex t;
    t.SetNDC(true);
    t.SetTextSize(0.032);
    t.DrawLatex(0.14, 0.88, Form("prompt: |dt| <= %d", prompt_halfwin));
    t.DrawLatex(0.14, 0.83, Form("sideband: %.0f <= |dt| <= %.0f", sideband_lo, sideband_hi));
    t.DrawLatex(0.14, 0.78, Form("scale s = width(prompt)/width(sideband) = %.3f", sf_side));
  }
  c->Print(out_pdf);

  // page: energy threshold decision
  c->Clear();
  c->Divide(2, 1);
  c->cd(1);
  gPad->SetGridx(); gPad->SetGridy(); gPad->SetLogy();
  hEPromptPre->SetLineColor(kBlue + 1);
  hEPromptPre->SetLineWidth(2);
  hESidePre->SetLineColor(kGreen + 2);
  hESidePre->SetLineWidth(1);
  TH1D* hESideScaledPre = (TH1D*)hESidePre->Clone("hESideScaledPre_stepwise");
  hESideScaledPre->Scale(sf_side);
  hESideScaledPre->SetLineColor(kRed + 1);
  hESideScaledPre->SetLineWidth(2);
  hEPromptPre->Draw("hist");
  hESidePre->Draw("hist same");
  hESideScaledPre->Draw("hist same");
  {
    TLine* le = new TLine(e_min, 0.1, e_min, std::max(1.0, hEPromptPre->GetMaximum()));
    le->SetLineColor(kBlack);
    le->SetLineStyle(2);
    le->SetLineWidth(2);
    le->Draw();
    TLegend* lg = new TLegend(0.48, 0.66, 0.92, 0.90);
    lg->SetBorderSize(0);
    lg->SetFillStyle(0);
    lg->AddEntry(hEPromptPre, "prompt", "l");
    lg->AddEntry(hESidePre, "sideband (raw)", "l");
    lg->AddEntry(hESideScaledPre, Form("scaled sideband (x %.3f)", sf_side), "l");
    lg->Draw();
    TLatex t;
    t.SetNDC(true);
    t.SetTextSize(0.030);
    t.DrawLatex(0.48, 0.60, "H_{excess}(E)=H_{prompt}(E)-s H_{side}(E)");
  }

  c->cd(2);
  gPad->SetGridx(); gPad->SetGridy();
  hEExcessPre->SetLineColor(kBlue + 1);
  hEExcessPre->SetLineWidth(2);
  hEExcessPre->Draw("hist");
  {
    TLine* l0 = new TLine(0.0, 0.0, 2.0e5, 0.0);
    l0->SetLineColor(kGray + 2);
    l0->SetLineStyle(2);
    l0->Draw();
    TLine* le = new TLine(e_min, hEExcessPre->GetMinimum(), e_min, hEExcessPre->GetMaximum());
    le->SetLineColor(kBlack);
    le->SetLineStyle(2);
    le->SetLineWidth(2);
    le->Draw();
    TLatex t;
    t.SetNDC(true);
    t.SetTextSize(0.030);
    t.DrawLatex(0.12, 0.90, "threshold at valley of excess histogram");
  }
  c->Print(out_pdf);

  // page: final |dt| distribution after all cuts (prompt/sideband windows)
  c->Clear();
  gPad->SetGridx(); gPad->SetGridy();
  hAbsDtFinal->SetLineColor(kBlue + 1);
  hAbsDtFinal->SetLineWidth(2);
  hAbsDtFinal->Draw("hist");
  {
    const double ymax = std::max(1.0, hAbsDtFinal->GetMaximum());
    TLine* lp = new TLine((double)prompt_halfwin, 0.0, (double)prompt_halfwin, ymax);
    lp->SetLineColor(kRed + 1);
    lp->SetLineStyle(2);
    lp->SetLineWidth(2);
    lp->Draw();
    TLine* ls1 = new TLine(sideband_lo, 0.0, sideband_lo, ymax);
    TLine* ls2 = new TLine(sideband_hi, 0.0, sideband_hi, ymax);
    ls1->SetLineColor(kGreen + 2); ls2->SetLineColor(kGreen + 2);
    ls1->SetLineStyle(3); ls2->SetLineStyle(3);
    ls1->SetLineWidth(2); ls2->SetLineWidth(2);
    ls1->Draw(); ls2->Draw();
    TLegend* lg = new TLegend(0.50, 0.74, 0.92, 0.90);
    lg->SetBorderSize(0);
    lg->SetFillStyle(0);
    lg->AddEntry(hAbsDtFinal, "final |dt|", "l");
    lg->AddEntry(lp, Form("prompt: <= %d", prompt_halfwin), "l");
    lg->AddEntry(ls1, Form("sideband: %.0f-%.0f", sideband_lo, sideband_hi), "l");
    lg->Draw();
  }
  c->Print(out_pdf);

  // page: final Michel overlay
  c->Clear();
  c->Divide(2, 1);
  c->cd(1);
  gPad->SetGridx(); gPad->SetGridy(); gPad->SetLogy();
  hEPrompt->SetLineColor(kBlue + 1);
  hESide->SetLineColor(kRed + 1);
  hEPrompt->SetLineWidth(2);
  hESide->SetLineWidth(2);
  hEPrompt->Draw("hist");
  hESide->Draw("hist same");
  {
    TLegend* lg = new TLegend(0.52, 0.74, 0.92, 0.90);
    lg->SetBorderSize(0);
    lg->SetFillStyle(0);
    lg->AddEntry(hEPrompt, "prompt (final)", "l");
    lg->AddEntry(hESide, "sideband (final)", "l");
    lg->Draw();
  }

  c->cd(2);
  gPad->SetGridx(); gPad->SetGridy();
  hEExcess->SetLineColor(kBlue + 1);
  hEExcess->SetLineWidth(2);
  hEExcess->Draw("hist");
  TF1 fMichelAll("fMichelAllDraw_stepwise", MichelAffinePdfStepwise, std::max(0.0, e_min), 2.0e5, 4);
  fMichelAll.SetParameters(fit_all.n, fit_all.e0, fit_all.w, fit_all.c);
  fMichelAll.SetLineColor(kRed + 1);
  fMichelAll.SetLineWidth(2);
  fMichelAll.Draw("same");
  {
    TLatex t;
    t.SetNDC(true);
    t.SetTextSize(0.035);
    t.DrawLatex(0.12, 0.90, Form("fit ok=%d, #chi^{2}/ndf=%.3f, R^{2}=%.3f",
                                 (fit_all.fit_ok ? 1 : 0), fit_all.chi2ndf, fit_all.r2));
    t.DrawLatex(0.12, 0.84, Form("E0=%.1f, W=%.1f", fit_all.e0, fit_all.w));
    TLine* l0 = new TLine(0.0, 0.0, 2.0e5, 0.0);
    l0->SetLineColor(kGray + 2);
    l0->SetLineStyle(2);
    l0->Draw();
  }
  c->Print(out_pdf);

  // page: side-separated Michel overlay (A/B)
  c->Clear();
  c->Divide(2, 1);
  c->cd(1);
  gPad->SetGridx(); gPad->SetGridy();
  hEExcessA->SetLineColor(kBlue + 1);
  hEExcessA->SetLineWidth(2);
  hEExcessA->Draw("hist");
  {
    TF1 fMichelA("fMichelADraw_stepwise", MichelAffinePdfStepwise, std::max(0.0, e_min), 2.0e5, 4);
    fMichelA.SetParameters(fit_a.n, fit_a.e0, fit_a.w, fit_a.c);
    fMichelA.SetLineColor(kRed + 1);
    fMichelA.SetLineWidth(2);
    fMichelA.Draw("same");
    TLatex t;
    t.SetNDC(true);
    t.SetTextSize(0.032);
    t.DrawLatex(0.12, 0.90, Form("A: fit ok=%d, #chi^{2}/ndf=%.3f, R^{2}=%.3f",
                                 (fit_a.fit_ok ? 1 : 0), fit_a.chi2ndf, fit_a.r2));
    t.DrawLatex(0.12, 0.84, Form("E0=%.1f, W=%.1f", fit_a.e0, fit_a.w));
    TLine* l0 = new TLine(0.0, 0.0, 2.0e5, 0.0);
    l0->SetLineColor(kGray + 2);
    l0->SetLineStyle(2);
    l0->Draw();
  }

  c->cd(2);
  gPad->SetGridx(); gPad->SetGridy();
  hEExcessB->SetLineColor(kBlue + 1);
  hEExcessB->SetLineWidth(2);
  hEExcessB->Draw("hist");
  {
    TF1 fMichelB("fMichelBDraw_stepwise", MichelAffinePdfStepwise, std::max(0.0, e_min), 2.0e5, 4);
    fMichelB.SetParameters(fit_b.n, fit_b.e0, fit_b.w, fit_b.c);
    fMichelB.SetLineColor(kRed + 1);
    fMichelB.SetLineWidth(2);
    fMichelB.Draw("same");
    TLatex t;
    t.SetNDC(true);
    t.SetTextSize(0.032);
    t.DrawLatex(0.12, 0.90, Form("B: fit ok=%d, #chi^{2}/ndf=%.3f, R^{2}=%.3f",
                                 (fit_b.fit_ok ? 1 : 0), fit_b.chi2ndf, fit_b.r2));
    t.DrawLatex(0.12, 0.84, Form("E0=%.1f, W=%.1f", fit_b.e0, fit_b.w));
    TLine* l0 = new TLine(0.0, 0.0, 2.0e5, 0.0);
    l0->SetLineColor(kGray + 2);
    l0->SetLineStyle(2);
    l0->Draw();
  }
  c->Print(out_pdf);

  // page: Michel fit residuals (all / A / B)
  TH1D* hPullAll = new TH1D("hPullAll_stepwise", "Pull (all);E_{sum} [ADC*samples];(Data-Fit)/#sigma",
                            500, 0, 2.0e5);
  TH1D* hPullA = new TH1D("hPullA_stepwise", "Pull (A);E_{sum} [ADC*samples];(Data-Fit)/#sigma",
                          500, 0, 2.0e5);
  TH1D* hPullB = new TH1D("hPullB_stepwise", "Pull (B);E_{sum} [ADC*samples];(Data-Fit)/#sigma",
                          500, 0, 2.0e5);
  TF1 fMichelAllPull("fMichelAllPull_stepwise", MichelAffinePdfStepwise, std::max(0.0, e_min), 2.0e5, 4);
  TF1 fMichelAPull("fMichelAPull_stepwise", MichelAffinePdfStepwise, std::max(0.0, e_min), 2.0e5, 4);
  TF1 fMichelBPull("fMichelBPull_stepwise", MichelAffinePdfStepwise, std::max(0.0, e_min), 2.0e5, 4);
  fMichelAllPull.SetParameters(fit_all.n, fit_all.e0, fit_all.w, fit_all.c);
  fMichelAPull.SetParameters(fit_a.n, fit_a.e0, fit_a.w, fit_a.c);
  fMichelBPull.SetParameters(fit_b.n, fit_b.e0, fit_b.w, fit_b.c);

  const int b0e = hEExcess->GetXaxis()->FindBin(std::max(0.0, e_min));
  const int b1e = hEExcess->GetNbinsX();
  for (int b = b0e; b <= b1e; ++b) {
    const double x = hEExcess->GetBinCenter(b);
    const double y = hEExcess->GetBinContent(b);
    const double ey = std::max(1.0, hEExcess->GetBinError(b));
    hPullAll->Fill(x, (y - fMichelAllPull.Eval(x)) / ey);
    const double yA = hEExcessA->GetBinContent(b);
    const double eyA = std::max(1.0, hEExcessA->GetBinError(b));
    hPullA->Fill(x, (yA - fMichelAPull.Eval(x)) / eyA);
    const double yB = hEExcessB->GetBinContent(b);
    const double eyB = std::max(1.0, hEExcessB->GetBinError(b));
    hPullB->Fill(x, (yB - fMichelBPull.Eval(x)) / eyB);
  }

  c->Clear();
  c->Divide(3, 1);
  c->cd(1);
  gPad->SetGridx(); gPad->SetGridy();
  hPullAll->SetLineColor(kBlue + 1);
  hPullAll->SetLineWidth(2);
  hPullAll->Draw("hist");
  {
    TLine* l0 = new TLine(std::max(0.0, e_min), 0.0, 2.0e5, 0.0);
    l0->SetLineStyle(2);
    l0->SetLineColor(kGray + 2);
    l0->Draw();
  }
  c->cd(2);
  gPad->SetGridx(); gPad->SetGridy();
  hPullA->SetLineColor(kBlue + 1);
  hPullA->SetLineWidth(2);
  hPullA->Draw("hist");
  {
    TLine* l0 = new TLine(std::max(0.0, e_min), 0.0, 2.0e5, 0.0);
    l0->SetLineStyle(2);
    l0->SetLineColor(kGray + 2);
    l0->Draw();
  }
  c->cd(3);
  gPad->SetGridx(); gPad->SetGridy();
  hPullB->SetLineColor(kBlue + 1);
  hPullB->SetLineWidth(2);
  hPullB->Draw("hist");
  {
    TLine* l0 = new TLine(std::max(0.0, e_min), 0.0, 2.0e5, 0.0);
    l0->SetLineStyle(2);
    l0->SetLineColor(kGray + 2);
    l0->Draw();
  }
  c->Print(out_pdf);

  c->Print(out_pdf + "]");
  delete c;

  std::cout << "[stepwise] output pdf: " << out_pdf << std::endl;
}
