// check/study_nai_pair_timediff.C
//
// 目的:
//   NaI 2チャネルの対応付けに使う |Δt|<30 samples の妥当性を、
//   NaIヒット判定閾値（seed/end）を走査しながら検証する。
//
// 出力:
//   doc/mainexp/study_nai_pair_timediff_threshold_scan_120deg_run8.pdf
//   doc/mainexp/study_nai_pair_timediff_threshold_scan_120deg_run8.txt
//
// 実行:
//   root -l -b -q 'check/study_nai_pair_timediff.C'

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
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

// 既存解析と同一の再構成手順・定数を使う
#include "study_michel_ps_tagged.C"

struct DtStats {
  long n_ref = 0;
  long n_le10 = 0;
  long n_le20 = 0;
  long n_le30 = 0;
  long n_le40 = 0;
  std::vector<double> absdt;
};

struct AllCombStats {
  long n_all = 0;
  long n_le5 = 0;
  long n_le10 = 0;
  long n_le20 = 0;
  std::vector<double> absdt;
};

struct AccStats {
  long n_acc = 0;
  long n_le10 = 0;
  long n_le20 = 0;
  long n_le25 = 0;
  long n_le30 = 0;
  std::vector<double> absdt;
};

struct ThresholdScanConfig {
  double seed_thr = 0.0;
  double end_thr = 0.0;
  std::string tag;
  int color = kBlack;
};

struct CuspMetric {
  double n0 = 0.0;
  double side_mean = 0.0;
  double center_over_side = 0.0;
  double frac_abs1 = 0.0;
  double frac_abs2 = 0.0;
  double frac_abs5 = 0.0;
};

struct GaussianFitResult {
  bool ok = false;
  double amp = 0.0;
  double mean = 0.0;
  double sigma = 0.0;
  double chi2ndf = 0.0;
};

struct ScanResult {
  ThresholdScanConfig cfg;
  DtStats stA;
  DtStats stB;
  AllCombStats stAllA;
  AllCombStats stAllB;
  AccStats stAccA;
  AccStats stAccB;
  long n_events = 0;
  long nA_ref = 0;
  long nB_ref = 0;
  long nA_acc = 0;
  long nB_acc = 0;
  CuspMetric cuspA;
  CuspMetric cuspB;
  GaussianFitResult gausA;
  GaussianFitResult gausB;
  TH1D* hAccA_norm = nullptr; // overlay用
  TH1D* hAccB_norm = nullptr; // overlay用
};

static int FindNearestPulseIndex(const std::vector<FitPulse>& v, double tref)
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

static std::vector<double> BuildPairDt12LikeAnalysis(const std::vector<FitPulse>& c1,
                                                      const std::vector<FitPulse>& c2,
                                                      int pair_dt_max)
{
  // study_michel_ps_tagged.C の BuildModulePairs と同じ対応規則で、
  // 受理ペアの Δt12 = t2 - t1 を返す。
  std::vector<double> out;
  if (c1.empty() || c2.empty()) return out;

  std::vector<int> i1(c1.size()), i2(c2.size());
  for (size_t i = 0; i < c1.size(); ++i) i1[i] = (int)i;
  for (size_t i = 0; i < c2.size(); ++i) i2[i] = (int)i;

  std::sort(i1.begin(), i1.end(), [&](int a, int b) { return c1[a].t0_fit < c1[b].t0_fit; });
  std::sort(i2.begin(), i2.end(), [&](int a, int b) { return c2[a].t0_fit < c2[b].t0_fit; });

  std::vector<char> used2(c2.size(), 0);
  for (int ia : i1) {
    int jbest = -1;
    double dmin = 1.0e30;
    for (int ib : i2) {
      if (used2[ib]) continue;
      const double dt12 = c2[ib].t0_fit - c1[ia].t0_fit;
      const double adt = std::abs(dt12);
      if (adt <= (double)pair_dt_max && adt < dmin) {
        dmin = adt;
        jbest = ib;
      }
    }
    if (jbest < 0) continue;
    used2[jbest] = 1;
    out.push_back(c2[jbest].t0_fit - c1[ia].t0_fit);
  }
  return out;
}

static double QuantileFromSorted(const std::vector<double>& x_sorted, double q)
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

static void UpdateStats(DtStats& st, double absdt)
{
  st.n_ref++;
  if (absdt <= 10.0) st.n_le10++;
  if (absdt <= 20.0) st.n_le20++;
  if (absdt <= 30.0) st.n_le30++;
  if (absdt <= 40.0) st.n_le40++;
  st.absdt.push_back(absdt);
}

static void UpdateAccStats(AccStats& st, double absdt)
{
  st.n_acc++;
  if (absdt <= 10.0) st.n_le10++;
  if (absdt <= 20.0) st.n_le20++;
  if (absdt <= 25.0) st.n_le25++;
  if (absdt <= 30.0) st.n_le30++;
  st.absdt.push_back(absdt);
}

static void UpdateAllCombStats(AllCombStats& st, double absdt)
{
  st.n_all++;
  if (absdt <= 5.0) st.n_le5++;
  if (absdt <= 10.0) st.n_le10++;
  if (absdt <= 20.0) st.n_le20++;
  st.absdt.push_back(absdt);
}

static void PrintStats(const char* tag, const DtStats& st)
{
  std::vector<double> v = st.absdt;
  std::sort(v.begin(), v.end());
  const double q50 = QuantileFromSorted(v, 0.50);
  const double q68 = QuantileFromSorted(v, 0.68);
  const double q90 = QuantileFromSorted(v, 0.90);
  const double q95 = QuantileFromSorted(v, 0.95);
  auto frac = [&](long n) { return (st.n_ref > 0) ? (double)n / (double)st.n_ref : 0.0; };

  std::cout << "[pairdt] " << tag << " n_ref=" << st.n_ref
            << " <=10:" << st.n_le10 << "(" << frac(st.n_le10) << ")"
            << " <=20:" << st.n_le20 << "(" << frac(st.n_le20) << ")"
            << " <=30:" << st.n_le30 << "(" << frac(st.n_le30) << ")"
            << " <=40:" << st.n_le40 << "(" << frac(st.n_le40) << ")"
            << " q50=" << q50
            << " q68=" << q68
            << " q90=" << q90
            << " q95=" << q95
            << std::endl;
}

static void PrintAllCombStats(const char* tag, const AllCombStats& st)
{
  std::vector<double> v = st.absdt;
  std::sort(v.begin(), v.end());
  const double q50 = QuantileFromSorted(v, 0.50);
  const double q68 = QuantileFromSorted(v, 0.68);
  const double q90 = QuantileFromSorted(v, 0.90);
  auto frac = [&](long n) { return (st.n_all > 0) ? (double)n / (double)st.n_all : 0.0; };

  std::cout << "[pairdt] " << tag << " all-comb n=" << st.n_all
            << " <=5:" << st.n_le5 << "(" << frac(st.n_le5) << ")"
            << " <=10:" << st.n_le10 << "(" << frac(st.n_le10) << ")"
            << " <=20:" << st.n_le20 << "(" << frac(st.n_le20) << ")"
            << " q50=" << q50
            << " q68=" << q68
            << " q90=" << q90
            << std::endl;
}

static void PrintAccStats(const char* tag, const AccStats& st)
{
  std::vector<double> v = st.absdt;
  std::sort(v.begin(), v.end());
  const double q50 = QuantileFromSorted(v, 0.50);
  const double q90 = QuantileFromSorted(v, 0.90);
  const double q95 = QuantileFromSorted(v, 0.95);
  auto frac = [&](long n) { return (st.n_acc > 0) ? (double)n / (double)st.n_acc : 0.0; };

  std::cout << "[pairdt] " << tag << " accepted n=" << st.n_acc
            << " <=10:" << st.n_le10 << "(" << frac(st.n_le10) << ")"
            << " <=20:" << st.n_le20 << "(" << frac(st.n_le20) << ")"
            << " <=25:" << st.n_le25 << "(" << frac(st.n_le25) << ")"
            << " <=30:" << st.n_le30 << "(" << frac(st.n_le30) << ")"
            << " q50=" << q50
            << " q90=" << q90
            << " q95=" << q95
            << std::endl;
}

static void AddPeakSupplementCandidatesWithThreshold(const std::vector<double>& s,
                                                     std::vector<int>& starts,
                                                     double seed_thr)
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
    if (s[i] < seed_thr + 10.0) continue;

    int t0 = i;
    const int jmin = std::max(1, i - 220);
    bool found = false;
    for (int j = i; j >= jmin; --j) {
      if (s[j - 1] < seed_thr && s[j] >= seed_thr) {
        t0 = j;
        found = true;
        break;
      }
    }
    if (!found) t0 = std::max(1, i - 24);
    if (!has_near(t0, 10)) starts.push_back(t0);
  }
}

static std::vector<FitPulse> FitPulsesWithDecompositionWithThreshold(const std::vector<double>& s_input,
                                                                     double tr_fix,
                                                                     double d_fix,
                                                                     double seed_thr,
                                                                     double end_thr)
{
  // study_michel_ps_tagged.C と同じ分解手順を使い、seed/endのみ可変化する。
  std::vector<FitPulse> accepted;
  if (s_input.empty()) return accepted;

  std::vector<double> s_work = s_input;
  std::vector<int> tried;

  for (int iter = 0; iter < kMaxIterPerCh; ++iter) {
    if ((int)accepted.size() >= kMaxAcceptedPerCh) break;

    std::vector<int> starts = FindStartCandidates(s_work, seed_thr, end_thr, kNaIMinSep);
    AddPeakSupplementCandidatesWithThreshold(s_work, starts, seed_thr);
    if (starts.empty()) break;

    std::sort(starts.begin(), starts.end());
    starts.erase(std::unique(starts.begin(), starts.end()), starts.end());

    int best_t0 = -1;
    double best_sc = -1.0;
    for (int t0 : starts) {
      if (t0 - kFitPreWindow < 0) continue;
      if (t0 + kFitHalfWindow >= (int)s_work.size()) continue;
      if (HasNearbyT0(tried, t0, 10)) continue;
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

static double IntegralAbsRange(const TH1D* h, double abs_min, double abs_max)
{
  if (!h) return 0.0;
  if (abs_max < abs_min) return 0.0;

  const int nbin = h->GetNbinsX();
  double s = 0.0;
  for (int b = 1; b <= nbin; ++b) {
    const double x = std::abs(h->GetBinCenter(b));
    if (x >= abs_min && x <= abs_max) s += h->GetBinContent(b);
  }
  return s;
}

static CuspMetric EvaluateCuspMetric(const TH1D* h)
{
  CuspMetric m;
  if (!h) return m;

  const int nbin = h->GetNbinsX();
  if (nbin <= 0) return m;

  const int b0 = h->GetXaxis()->FindBin(0.0);
  m.n0 = h->GetBinContent(b0);

  double side_sum = 0.0;
  int side_cnt = 0;
  // 中心ビン近傍(±1)は除外して、裾の平均高さを作る。
  for (int k = 2; k <= 8; ++k) {
    const int bl = b0 - k;
    const int br = b0 + k;
    if (bl >= 1) {
      side_sum += h->GetBinContent(bl);
      side_cnt++;
    }
    if (br <= nbin) {
      side_sum += h->GetBinContent(br);
      side_cnt++;
    }
  }
  m.side_mean = (side_cnt > 0) ? (side_sum / (double)side_cnt) : 0.0;
  m.center_over_side = (m.side_mean > 0.0) ? (m.n0 / m.side_mean) : 0.0;

  const double nall = h->Integral(1, nbin);
  if (nall > 0.0) {
    m.frac_abs1 = IntegralAbsRange(h, 0.0, 1.0) / nall;
    m.frac_abs2 = IntegralAbsRange(h, 0.0, 2.0) / nall;
    m.frac_abs5 = IntegralAbsRange(h, 0.0, 5.0) / nall;
  }
  return m;
}

static GaussianFitResult FitGaussianAroundZero(TH1D* h, double xmin, double xmax)
{
  GaussianFitResult r;
  if (!h) return r;
  if (h->GetEntries() < 100.0) return r;

  static int fit_counter = 0;
  const TString f_name = Form("f_pairdt_gaus_%d", fit_counter++);
  TF1 f(f_name, "gaus", xmin, xmax);
  f.SetParameters(std::max(1.0, h->GetMaximum()), 0.0, 4.0);

  const int fit_status = h->Fit(&f, "Q0R");
  if (fit_status != 0) return r;
  if (!std::isfinite(f.GetParameter(2))) return r;

  r.ok = true;
  r.amp = f.GetParameter(0);
  r.mean = f.GetParameter(1);
  r.sigma = std::abs(f.GetParameter(2));
  r.chi2ndf = (f.GetNDF() > 0) ? (f.GetChisquare() / (double)f.GetNDF()) : 0.0;
  return r;
}

static void PrintScanHeadline(const ScanResult& r)
{
  std::cout << "[pairdt-scan] " << r.cfg.tag
            << " seed=" << r.cfg.seed_thr
            << " end=" << r.cfg.end_thr
            << " events=" << r.n_events
            << " A_acc=" << r.nA_acc
            << " B_acc=" << r.nB_acc
            << " A(cusp n0/side)=" << r.cuspA.center_over_side
            << " B(cusp n0/side)=" << r.cuspB.center_over_side
            << std::endl;
}

static ScanResult AnalyzeOneThreshold(const ThresholdScanConfig& cfg, size_t idx, TCanvas* c, const TString& out_pdf)
{
  ScanResult out;
  out.cfg = cfg;

  TH1D* hDtNearA = new TH1D(Form("hDtNearA_%zu", idx),
                            Form("A nearest (seed=%.0f,end=%.0f): #Delta t_{12}=t_{A2}-t_{A1};#Delta t_{12} [samples];Counts",
                                 cfg.seed_thr, cfg.end_thr),
                            400, -200, 200);
  TH1D* hDtNearB = new TH1D(Form("hDtNearB_%zu", idx),
                            Form("B nearest (seed=%.0f,end=%.0f): #Delta t_{12}=t_{B2}-t_{B1};#Delta t_{12} [samples];Counts",
                                 cfg.seed_thr, cfg.end_thr),
                            400, -200, 200);
  TH1D* hDtAccA = new TH1D(Form("hDtAccA_%zu", idx),
                           Form("A accepted pair (|#Delta t|#le%d, seed=%.0f,end=%.0f);#Delta t_{12} [samples];Counts",
                                kPairMatchMaxDt, cfg.seed_thr, cfg.end_thr),
                           400, -200, 200);
  TH1D* hDtAccB = new TH1D(Form("hDtAccB_%zu", idx),
                           Form("B accepted pair (|#Delta t|#le%d, seed=%.0f,end=%.0f);#Delta t_{12} [samples];Counts",
                                kPairMatchMaxDt, cfg.seed_thr, cfg.end_thr),
                           400, -200, 200);

  TH1D* hAbsNearA = new TH1D(Form("hAbsNearA_%zu", idx),
                             "A nearest #left|#Delta t_{12}#right|;#left|#Delta t_{12}#right| [samples];Counts",
                             200, 0, 100);
  TH1D* hAbsNearB = new TH1D(Form("hAbsNearB_%zu", idx),
                             "B nearest #left|#Delta t_{12}#right|;#left|#Delta t_{12}#right| [samples];Counts",
                             200, 0, 100);
  TH1D* hAbsAllA = new TH1D(Form("hAbsAllA_%zu", idx),
                            "A all combinations #left|#Delta t_{12}#right|;#left|#Delta t_{12}#right| [samples];Counts",
                            200, 0, 100);
  TH1D* hAbsAllB = new TH1D(Form("hAbsAllB_%zu", idx),
                            "B all combinations #left|#Delta t_{12}#right|;#left|#Delta t_{12}#right| [samples];Counts",
                            200, 0, 100);

  TString fna[4];
  for (int i = 0; i < 4; ++i) fna[i] = Form("%s/%s", kInputDir, kNaIFiles[i]);

  std::ifstream in_na[4];
  for (int i = 0; i < 4; ++i) {
    in_na[i].open(fna[i].Data());
    if (!in_na[i]) {
      std::cerr << "[pairdt] ERROR: cannot open " << fna[i] << std::endl;
      delete hDtNearA; delete hDtNearB; delete hDtAccA; delete hDtAccB;
      delete hAbsNearA; delete hAbsNearB; delete hAbsAllA; delete hAbsAllB;
      return out;
    }
  }

  std::vector<double> buf_na[4];
  for (int i = 0; i < 4; ++i) buf_na[i].reserve(kSamplesPerEvent);

  while (true) {
    double xna[4];
    bool ok = true;
    for (int i = 0; i < 4; ++i) {
      if (!(in_na[i] >> xna[i])) {
        ok = false;
        break;
      }
    }
    if (!ok) break;

    for (int i = 0; i < 4; ++i) buf_na[i].push_back(xna[i]);
    if ((int)buf_na[0].size() < kSamplesPerEvent) continue;

    out.n_events++;

    std::vector<FitPulse> fp[4];
    for (int ch = 0; ch < 4; ++ch) {
      const double b = ModeBaselineInEvent(buf_na[ch]);
      std::vector<double> s;
      BuildS(buf_na[ch], b, s);
      fp[ch] = FitPulsesWithDecompositionWithThreshold(s, kTrFix[ch], kDFix[ch], cfg.seed_thr, cfg.end_thr);
    }

    if (!fp[0].empty() && !fp[1].empty()) {
      for (const auto& p1 : fp[0]) {
        const int j = FindNearestPulseIndex(fp[1], p1.t0_fit);
        if (j < 0) continue;
        const double dt12 = fp[1][j].t0_fit - p1.t0_fit;
        const double adt = std::abs(dt12);
        hDtNearA->Fill(dt12);
        hAbsNearA->Fill(adt);
        UpdateStats(out.stA, adt);
        out.nA_ref++;
      }
      for (const auto& p1 : fp[0]) {
        for (const auto& p2 : fp[1]) {
          const double adt = std::abs(p2.t0_fit - p1.t0_fit);
          hAbsAllA->Fill(adt);
          UpdateAllCombStats(out.stAllA, adt);
        }
      }
    }

    if (!fp[2].empty() && !fp[3].empty()) {
      for (const auto& p1 : fp[2]) {
        const int j = FindNearestPulseIndex(fp[3], p1.t0_fit);
        if (j < 0) continue;
        const double dt12 = fp[3][j].t0_fit - p1.t0_fit;
        const double adt = std::abs(dt12);
        hDtNearB->Fill(dt12);
        hAbsNearB->Fill(adt);
        UpdateStats(out.stB, adt);
        out.nB_ref++;
      }
      for (const auto& p1 : fp[2]) {
        for (const auto& p2 : fp[3]) {
          const double adt = std::abs(p2.t0_fit - p1.t0_fit);
          hAbsAllB->Fill(adt);
          UpdateAllCombStats(out.stAllB, adt);
        }
      }
    }

    const std::vector<double> dt_acc_A = BuildPairDt12LikeAnalysis(fp[0], fp[1], kPairMatchMaxDt);
    const std::vector<double> dt_acc_B = BuildPairDt12LikeAnalysis(fp[2], fp[3], kPairMatchMaxDt);

    for (double dt : dt_acc_A) {
      hDtAccA->Fill(dt);
      out.nA_acc++;
      UpdateAccStats(out.stAccA, std::abs(dt));
    }
    for (double dt : dt_acc_B) {
      hDtAccB->Fill(dt);
      out.nB_acc++;
      UpdateAccStats(out.stAccB, std::abs(dt));
    }

    for (int i = 0; i < 4; ++i) buf_na[i].clear();
  }

  PrintStats(Form("A nearest (%s)", cfg.tag.c_str()), out.stA);
  PrintStats(Form("B nearest (%s)", cfg.tag.c_str()), out.stB);
  PrintAllCombStats(Form("A (%s)", cfg.tag.c_str()), out.stAllA);
  PrintAllCombStats(Form("B (%s)", cfg.tag.c_str()), out.stAllB);
  PrintAccStats(Form("A (%s)", cfg.tag.c_str()), out.stAccA);
  PrintAccStats(Form("B (%s)", cfg.tag.c_str()), out.stAccB);

  out.cuspA = EvaluateCuspMetric(hDtAccA);
  out.cuspB = EvaluateCuspMetric(hDtAccB);
  out.gausA = FitGaussianAroundZero(hDtAccA, -20.0, 20.0);
  out.gausB = FitGaussianAroundZero(hDtAccB, -20.0, 20.0);
  PrintScanHeadline(out);

  c->Clear();
  c->Divide(2, 2);

  c->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  hDtNearA->SetLineColor(kGray + 2);
  hDtNearA->SetLineWidth(2);
  hDtAccA->SetLineColor(cfg.color);
  hDtAccA->SetLineWidth(2);
  hDtNearA->Draw("hist");
  hDtAccA->Draw("hist same");
  {
    TLine* l1 = new TLine(-kPairMatchMaxDt, 0.0, -kPairMatchMaxDt, std::max(1.0, hDtNearA->GetMaximum()));
    TLine* l2 = new TLine(+kPairMatchMaxDt, 0.0, +kPairMatchMaxDt, std::max(1.0, hDtNearA->GetMaximum()));
    l1->SetLineStyle(2); l2->SetLineStyle(2);
    l1->SetLineColor(kRed + 1); l2->SetLineColor(kRed + 1);
    l1->Draw(); l2->Draw();
    TLegend* lg = new TLegend(0.44, 0.72, 0.92, 0.90);
    lg->SetBorderSize(0);
    lg->SetFillStyle(0);
    lg->AddEntry(hDtNearA, "nearest pair (no |#Delta t| cut)", "l");
    lg->AddEntry(hDtAccA, "accepted pair (|#Delta t| cut)", "l");
    lg->Draw();
    TLatex tx;
    tx.SetNDC(true);
    tx.SetTextSize(0.038);
    tx.DrawLatex(0.12, 0.92, Form("%s: seed=%.0f end=%.0f", cfg.tag.c_str(), cfg.seed_thr, cfg.end_thr));
  }

  c->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  hDtNearB->SetLineColor(kGray + 2);
  hDtNearB->SetLineWidth(2);
  hDtAccB->SetLineColor(cfg.color);
  hDtAccB->SetLineWidth(2);
  hDtNearB->Draw("hist");
  hDtAccB->Draw("hist same");
  {
    TLine* l1 = new TLine(-kPairMatchMaxDt, 0.0, -kPairMatchMaxDt, std::max(1.0, hDtNearB->GetMaximum()));
    TLine* l2 = new TLine(+kPairMatchMaxDt, 0.0, +kPairMatchMaxDt, std::max(1.0, hDtNearB->GetMaximum()));
    l1->SetLineStyle(2); l2->SetLineStyle(2);
    l1->SetLineColor(kRed + 1); l2->SetLineColor(kRed + 1);
    l1->Draw(); l2->Draw();
    TLegend* lg = new TLegend(0.44, 0.72, 0.92, 0.90);
    lg->SetBorderSize(0);
    lg->SetFillStyle(0);
    lg->AddEntry(hDtNearB, "nearest pair (no |#Delta t| cut)", "l");
    lg->AddEntry(hDtAccB, "accepted pair (|#Delta t| cut)", "l");
    lg->Draw();
  }

  c->cd(3);
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogy();
  hAbsAllA->SetLineColor(kGray + 2);
  hAbsAllA->SetLineWidth(2);
  hAbsNearA->SetLineColor(cfg.color);
  hAbsNearA->SetLineWidth(2);
  hAbsAllA->Draw("hist");
  hAbsNearA->Draw("hist same");
  {
    TLine* l = new TLine(kPairMatchMaxDt, 0.5, kPairMatchMaxDt, std::max(1.0, hAbsNearA->GetMaximum()));
    l->SetLineColor(kRed + 1);
    l->SetLineStyle(2);
    l->SetLineWidth(2);
    l->Draw();
    TLegend* lg = new TLegend(0.50, 0.74, 0.92, 0.90);
    lg->SetBorderSize(0);
    lg->SetFillStyle(0);
    lg->AddEntry(hAbsAllA, "all combinations", "l");
    lg->AddEntry(hAbsNearA, "nearest pair", "l");
    lg->Draw();
  }

  c->cd(4);
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogy();
  hAbsAllB->SetLineColor(kGray + 2);
  hAbsAllB->SetLineWidth(2);
  hAbsNearB->SetLineColor(cfg.color);
  hAbsNearB->SetLineWidth(2);
  hAbsAllB->Draw("hist");
  hAbsNearB->Draw("hist same");
  {
    TLine* l = new TLine(kPairMatchMaxDt, 0.5, kPairMatchMaxDt, std::max(1.0, hAbsNearB->GetMaximum()));
    l->SetLineColor(kRed + 1);
    l->SetLineStyle(2);
    l->SetLineWidth(2);
    l->Draw();
    TLegend* lg = new TLegend(0.50, 0.74, 0.92, 0.90);
    lg->SetBorderSize(0);
    lg->SetFillStyle(0);
    lg->AddEntry(hAbsAllB, "all combinations", "l");
    lg->AddEntry(hAbsNearB, "nearest pair", "l");
    lg->Draw();
  }
  c->Print(out_pdf);

  c->Clear();
  c->Divide(2, 1);
  c->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  hDtAccA->SetLineColor(cfg.color);
  hDtAccA->SetLineWidth(2);
  hDtAccA->GetXaxis()->SetRangeUser(-40.0, 40.0);
  hDtAccA->Draw("hist");
  {
    TLine* l0 = new TLine(0.0, 0.0, 0.0, std::max(1.0, hDtAccA->GetMaximum()));
    l0->SetLineColor(kBlack);
    l0->SetLineStyle(2);
    l0->Draw();
    if (out.gausA.ok) {
      TF1* f = new TF1(Form("f_draw_A_%zu", idx), "gaus", -20.0, 20.0);
      f->SetParameters(out.gausA.amp, out.gausA.mean, out.gausA.sigma);
      f->SetLineColor(kRed + 1);
      f->SetLineWidth(2);
      f->Draw("same");
    }
    TLatex tx;
    tx.SetNDC(true);
    tx.SetTextSize(0.036);
    tx.DrawLatex(0.12, 0.90, Form("A: n0/side=%.2f, |dt|<=1 frac=%.3f", out.cuspA.center_over_side, out.cuspA.frac_abs1));
    if (out.gausA.ok) {
      tx.DrawLatex(0.12, 0.84, Form("Gaussian: #mu=%.2f, #sigma=%.2f, #chi^{2}/ndf=%.2f",
                                    out.gausA.mean, out.gausA.sigma, out.gausA.chi2ndf));
    } else {
      tx.DrawLatex(0.12, 0.84, "Gaussian fit: failed / insufficient entries");
    }
  }

  c->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  hDtAccB->SetLineColor(cfg.color);
  hDtAccB->SetLineWidth(2);
  hDtAccB->GetXaxis()->SetRangeUser(-40.0, 40.0);
  hDtAccB->Draw("hist");
  {
    TLine* l0 = new TLine(0.0, 0.0, 0.0, std::max(1.0, hDtAccB->GetMaximum()));
    l0->SetLineColor(kBlack);
    l0->SetLineStyle(2);
    l0->Draw();
    if (out.gausB.ok) {
      TF1* f = new TF1(Form("f_draw_B_%zu", idx), "gaus", -20.0, 20.0);
      f->SetParameters(out.gausB.amp, out.gausB.mean, out.gausB.sigma);
      f->SetLineColor(kRed + 1);
      f->SetLineWidth(2);
      f->Draw("same");
    }
    TLatex tx;
    tx.SetNDC(true);
    tx.SetTextSize(0.036);
    tx.DrawLatex(0.12, 0.90, Form("B: n0/side=%.2f, |dt|<=1 frac=%.3f", out.cuspB.center_over_side, out.cuspB.frac_abs1));
    if (out.gausB.ok) {
      tx.DrawLatex(0.12, 0.84, Form("Gaussian: #mu=%.2f, #sigma=%.2f, #chi^{2}/ndf=%.2f",
                                    out.gausB.mean, out.gausB.sigma, out.gausB.chi2ndf));
    } else {
      tx.DrawLatex(0.12, 0.84, "Gaussian fit: failed / insufficient entries");
    }
  }
  c->Print(out_pdf);

  out.hAccA_norm = (TH1D*)hDtAccA->Clone(Form("hAccA_norm_%zu", idx));
  out.hAccB_norm = (TH1D*)hDtAccB->Clone(Form("hAccB_norm_%zu", idx));
  out.hAccA_norm->SetDirectory(nullptr);
  out.hAccB_norm->SetDirectory(nullptr);

  const double intA = out.hAccA_norm->Integral(1, out.hAccA_norm->GetNbinsX());
  const double intB = out.hAccB_norm->Integral(1, out.hAccB_norm->GetNbinsX());
  if (intA > 0.0) out.hAccA_norm->Scale(1.0 / intA);
  if (intB > 0.0) out.hAccB_norm->Scale(1.0 / intB);
  out.hAccA_norm->SetLineColor(cfg.color);
  out.hAccA_norm->SetLineWidth(2);
  out.hAccB_norm->SetLineColor(cfg.color);
  out.hAccB_norm->SetLineWidth(2);

  delete hDtNearA;
  delete hDtNearB;
  delete hDtAccA;
  delete hDtAccB;
  delete hAbsNearA;
  delete hAbsNearB;
  delete hAbsAllA;
  delete hAbsAllB;
  return out;
}

void study_nai_pair_timediff()
{
  gStyle->SetOptStat(1110);
  gSystem->mkdir("doc/mainexp", true);

  const TString out_pdf = "doc/mainexp/study_nai_pair_timediff_threshold_scan_120deg_run8.pdf";
  const TString out_txt = "doc/mainexp/study_nai_pair_timediff_threshold_scan_120deg_run8.txt";

  // 既定値(60/30)を含めて、NaIヒット閾値を厳しくして走査する。
  const std::vector<ThresholdScanConfig> cfgs = {
    {50.0, 25.0, "thr50_25",  kGray + 2},
    {60.0, 30.0, "thr60_30",  kBlack},
    {70.0, 35.0, "thr70_35",  kBlue + 1},
    {80.0, 40.0, "thr80_40",  kRed + 1},
    {90.0, 45.0, "thr90_45",  kGreen + 2},
    {100.0, 50.0, "thr100_50", kMagenta + 2},
    {120.0, 60.0, "thr120_60", kOrange + 7}
  };

  std::vector<ScanResult> results;
  results.reserve(cfgs.size());

  TCanvas* c = new TCanvas("c_pairdt_scan", "nai pair timediff threshold scan", 1400, 1000);
  c->Print(out_pdf + "[");

  for (size_t i = 0; i < cfgs.size(); ++i) {
    std::cout << "[pairdt-scan] start " << cfgs[i].tag
              << " (seed=" << cfgs[i].seed_thr
              << ", end=" << cfgs[i].end_thr << ")" << std::endl;
    ScanResult r = AnalyzeOneThreshold(cfgs[i], i, c, out_pdf);
    results.push_back(r);
  }

  // A/B accepted分布の閾値依存を重ね描き。
  c->Clear();
  c->Divide(2, 1);

  c->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  TLegend* lgA = new TLegend(0.58, 0.60, 0.92, 0.90);
  lgA->SetBorderSize(0);
  lgA->SetFillStyle(0);
  bool drawnA = false;
  for (const auto& r : results) {
    if (!r.hAccA_norm) continue;
    r.hAccA_norm->SetTitle("A accepted #Delta t overlay (normalized);#Delta t_{12} [samples];Normalized counts");
    r.hAccA_norm->GetXaxis()->SetRangeUser(-40.0, 40.0);
    if (!drawnA) {
      r.hAccA_norm->Draw("hist");
      drawnA = true;
    } else {
      r.hAccA_norm->Draw("hist same");
    }
    lgA->AddEntry(r.hAccA_norm, Form("seed=%.0f end=%.0f", r.cfg.seed_thr, r.cfg.end_thr), "l");
  }
  if (drawnA) lgA->Draw();

  c->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  TLegend* lgB = new TLegend(0.58, 0.60, 0.92, 0.90);
  lgB->SetBorderSize(0);
  lgB->SetFillStyle(0);
  bool drawnB = false;
  for (const auto& r : results) {
    if (!r.hAccB_norm) continue;
    r.hAccB_norm->SetTitle("B accepted #Delta t overlay (normalized);#Delta t_{12} [samples];Normalized counts");
    r.hAccB_norm->GetXaxis()->SetRangeUser(-40.0, 40.0);
    if (!drawnB) {
      r.hAccB_norm->Draw("hist");
      drawnB = true;
    } else {
      r.hAccB_norm->Draw("hist same");
    }
    lgB->AddEntry(r.hAccB_norm, Form("seed=%.0f end=%.0f", r.cfg.seed_thr, r.cfg.end_thr), "l");
  }
  if (drawnB) lgB->Draw();
  c->Print(out_pdf);

  // まとめページ
  c->Clear();
  TPaveText* pt = new TPaveText(0.04, 0.04, 0.96, 0.96, "NDC");
  pt->SetFillStyle(0);
  pt->SetBorderSize(0);
  pt->SetTextAlign(12);
  pt->SetTextSize(0.025);
  pt->AddText("=== NaI pair #Delta t threshold scan summary ===");
  pt->AddText(Form("Input: %s", kInputDir));
  pt->AddText(Form("Pair criterion: |#Delta t_{12}| <= %d samples", kPairMatchMaxDt));
  pt->AddText("columns: seed/end, accepted(A/B), cusp ratio n0/<side>, frac(|dt|<=1), Gaussian fit");
  pt->AddText(" ");
  for (const auto& r : results) {
    const TString ga = r.gausA.ok ? Form("A(#mu=%.2f,#sigma=%.2f,#chi2/ndf=%.2f)", r.gausA.mean, r.gausA.sigma, r.gausA.chi2ndf)
                                   : "A(fit failed)";
    const TString gb = r.gausB.ok ? Form("B(#mu=%.2f,#sigma=%.2f,#chi2/ndf=%.2f)", r.gausB.mean, r.gausB.sigma, r.gausB.chi2ndf)
                                   : "B(fit failed)";
    pt->AddText(Form("seed=%.0f end=%.0f | acc A/B=%ld/%ld | cusp A/B=%.2f/%.2f | frac|dt|<=1 A/B=%.3f/%.3f | %s %s",
                     r.cfg.seed_thr, r.cfg.end_thr, r.nA_acc, r.nB_acc,
                     r.cuspA.center_over_side, r.cuspB.center_over_side,
                     r.cuspA.frac_abs1, r.cuspB.frac_abs1,
                     ga.Data(), gb.Data()));
  }
  pt->AddText(" ");
  pt->AddText("Interpretation note:");
  pt->AddText(" - n0/<side> ~ 1 かつ #chi2/ndf が過大でなければ、中心が尖り過ぎない(ガウシアン寄り)とみなせる。");
  pt->AddText(" - n0/<side> >> 1 は t=0 カスプ優勢を示唆する。");
  pt->Draw();
  c->Print(out_pdf);
  delete pt;

  c->Print(out_pdf + "]");
  delete c;
  delete lgA;
  delete lgB;

  std::ofstream fout(out_txt.Data());
  if (!fout) {
    std::cerr << "[pairdt-scan] ERROR: cannot open " << out_txt << std::endl;
  } else {
    fout << "# NaI pair timediff threshold scan\n";
    fout << "# input_dir: " << kInputDir << "\n";
    fout << "# pair_cut : |dt| <= " << kPairMatchMaxDt << " samples\n";
    fout << "# columns:\n";
    fout << "# seed end n_events nA_acc nB_acc cuspA cuspB fracA_abs1 fracB_abs1 "
         << "gausA_ok gausA_mu gausA_sigma gausA_chi2ndf "
         << "gausB_ok gausB_mu gausB_sigma gausB_chi2ndf\n";
    fout << std::fixed << std::setprecision(6);
    for (const auto& r : results) {
      fout << r.cfg.seed_thr << " "
           << r.cfg.end_thr << " "
           << r.n_events << " "
           << r.nA_acc << " "
           << r.nB_acc << " "
           << r.cuspA.center_over_side << " "
           << r.cuspB.center_over_side << " "
           << r.cuspA.frac_abs1 << " "
           << r.cuspB.frac_abs1 << " "
           << (r.gausA.ok ? 1 : 0) << " "
           << r.gausA.mean << " "
           << r.gausA.sigma << " "
           << r.gausA.chi2ndf << " "
           << (r.gausB.ok ? 1 : 0) << " "
           << r.gausB.mean << " "
           << r.gausB.sigma << " "
           << r.gausB.chi2ndf
           << "\n";
    }
  }

  for (auto& r : results) {
    delete r.hAccA_norm;
    delete r.hAccB_norm;
    r.hAccA_norm = nullptr;
    r.hAccB_norm = nullptr;
  }

  std::cout << "[pairdt-scan] output pdf: " << out_pdf << std::endl;
  std::cout << "[pairdt-scan] output txt: " << out_txt << std::endl;
}
