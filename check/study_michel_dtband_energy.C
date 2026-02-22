// check/study_michel_dtband_energy.C
//
// 目的:
//   ModuleA(A1+A2), ModuleB(B1+B2) の nearest-pair dt 分布で見える
//   「主ピーク」と「遅延側ピーク(例: dt~170 samp)」を分けて、
//   それぞれの NaI 和積分スペクトラムを比較する。
//
// 実行:
//   root -l -q 'check/study_michel_dtband_energy.C'

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"

// ============================================================
// 手調整パラメータ
// ============================================================

static const char* kInputDir = "data/rawdata/mainexp/120°/run8";
static const char* kInputTag = "120°_run8";
static const char* kOutputDir = "doc/mainexp";

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

static const int    kSamplesPerEvent = 1000;
static const double kDtNs = 4.0;

// PS / NaI の簡易パルス抽出条件
static const double kPSSeedThr = 22.0;
static const double kPSEndThr  = 10.0;
static const int    kPSMinSep  = 12;
static const int    kPSMinWidth = 2;

static const double kNaISeedThr = 35.0;
static const double kNaIEndThr  = 12.0;
static const int    kNaIMinSep  = 12;
static const int    kNaIMinWidth = 4;

// NaI 積分
static const int kIntPreStart = 80;
static const int kIntPreEnd   = 20;
static const int kIntWinPre   = 8;
static const int kIntWinPost  = 240;

// A1-A2, B1-B2 ペアの一致幅
static const int kPairMatchMaxDt = 30;

// nearest pair を作るための最小 PS 条件
static const double kPSAmpMinForDt = 25.0;

// 「主ピーク」探索範囲（samples）
static const int kPromptSearchMin = -60;
static const int kPromptSearchMax = 60;

// 「遅延側ピーク」探索範囲（samples）
static const int kDelayedSearchMin = 80;
static const int kDelayedSearchMax = 260;

// 切り出し窓半幅（samples）
static const int kBandHalfWindow = 20;
static const double kManualDelayedCenter = 170.0;

static const int    kDtBins = 600;
static const double kDtMin = -600.0;
static const double kDtMax = 600.0;

static const int    kEBins = 2500;
static const double kEMin = 0.0;
static const double kEMax = 2.0e5;

static const bool kLogY = true;
static const int  kCanvasW = 1400;
static const int  kCanvasH = 1000;

// ============================================================

struct Pulse {
  int start = -1;
  int peak = -1;
  int end = -1;
  double amp = 0.0;
  double area = 0.0;
  double pre_offset = 0.0;
};

struct ModulePulse {
  int t = -1;
  int t1 = -1;
  int t2 = -1;
  double dt12 = 0.0;
  double area_sum = 0.0;
};

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

static double ComputeLocalArea(const std::vector<double>& s, int peak, double& pre_offset)
{
  const int n = (int)s.size();
  if (n <= 0 || peak < 0 || peak >= n) {
    pre_offset = 0.0;
    return 0.0;
  }
  pre_offset = MeanInRange(s, peak - kIntPreStart, peak - kIntPreEnd);
  const int i0 = std::max(0, peak - kIntWinPre);
  const int i1 = std::min(n - 1, peak + kIntWinPost);
  double area = 0.0;
  for (int i = i0; i <= i1; ++i) {
    const double y = s[i] - pre_offset;
    if (y > 0.0) area += y;
  }
  return area;
}

static std::vector<Pulse> FindPulses(const std::vector<double>& s,
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
      if (width >= min_width && cur.amp > 0.0) {
        cur.area = ComputeLocalArea(s, cur.peak, cur.pre_offset);
        out.push_back(cur);
      }
      in = false;
      dead = min_sep;
      amax = -1.0;
    }
  }

  if (in) {
    cur.end = n - 1;
    cur.amp = std::max(0.0, amax);
    const int width = cur.end - cur.start + 1;
    if (width >= min_width && cur.amp > 0.0) {
      cur.area = ComputeLocalArea(s, cur.peak, cur.pre_offset);
      out.push_back(cur);
    }
  }
  return out;
}

static std::vector<ModulePulse> BuildModulePairs(const std::vector<Pulse>& p1,
                                                 const std::vector<Pulse>& p2,
                                                 int pair_dt_max)
{
  std::vector<ModulePulse> out;
  if (p1.empty() || p2.empty()) return out;

  std::vector<int> i1(p1.size()), i2(p2.size());
  for (size_t i = 0; i < p1.size(); ++i) i1[i] = (int)i;
  for (size_t i = 0; i < p2.size(); ++i) i2[i] = (int)i;
  std::sort(i1.begin(), i1.end(), [&](int a, int b) { return p1[a].start < p1[b].start; });
  std::sort(i2.begin(), i2.end(), [&](int a, int b) { return p2[a].start < p2[b].start; });

  std::vector<char> used2(p2.size(), 0);
  for (int a : i1) {
    int jbest = -1;
    int dmin = std::numeric_limits<int>::max();
    for (int b : i2) {
      if (used2[b]) continue;
      const int d = std::abs(p1[a].start - p2[b].start);
      if (d <= pair_dt_max && d < dmin) {
        dmin = d;
        jbest = b;
      }
    }
    if (jbest < 0) continue;
    used2[jbest] = 1;
    ModulePulse m;
    m.t1 = p1[a].start;
    m.t2 = p2[jbest].start;
    m.t = (int)std::lround(0.5 * (double)(m.t1 + m.t2));
    m.dt12 = (double)(m.t1 - m.t2);
    m.area_sum = p1[a].area + p2[jbest].area;
    out.push_back(m);
  }

  std::sort(out.begin(), out.end(), [](const ModulePulse& a, const ModulePulse& b) {
    return a.t < b.t;
  });
  return out;
}

static int FindNearestModuleIndex(const std::vector<ModulePulse>& v, int tps)
{
  if (v.empty()) return -1;
  int ibest = -1;
  int dmin = std::numeric_limits<int>::max();
  for (size_t i = 0; i < v.size(); ++i) {
    const int d = std::abs(v[i].t - tps);
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

void study_michel_dtband_energy()
{
  gStyle->SetOptStat(1110);
  gStyle->SetTitleFontSize(0.040);

  gSystem->mkdir(kOutputDir, true);
  TString out_pdf = Form("%s/study_michel_dtband_energy_%s.pdf", kOutputDir, kInputTag);
  TString out_root = Form("%s/study_michel_dtband_energy_%s.root", kOutputDir, kInputTag);

  TString f_ps[2], f_na[4];
  for (int i = 0; i < 2; ++i) f_ps[i] = Form("%s/%s", kInputDir, kPSFiles[i]);
  for (int i = 0; i < 4; ++i) f_na[i] = Form("%s/%s", kInputDir, kNaIFiles[i]);

  std::ifstream in_ps[2], in_na[4];
  for (int i = 0; i < 2; ++i) {
    in_ps[i].open(f_ps[i].Data());
    if (!in_ps[i]) {
      std::cerr << "[dtband] ERROR: cannot open " << f_ps[i] << std::endl;
      return;
    }
  }
  for (int i = 0; i < 4; ++i) {
    in_na[i].open(f_na[i].Data());
    if (!in_na[i]) {
      std::cerr << "[dtband] ERROR: cannot open " << f_na[i] << std::endl;
      return;
    }
  }

  TH1D* hDtA = new TH1D("hDtA", "ModuleA nearest pair dt;dt=t_{mod}-t_{PS} [samples];Events",
                        kDtBins, kDtMin, kDtMax);
  TH1D* hDtB = new TH1D("hDtB", "ModuleB nearest pair dt;dt=t_{mod}-t_{PS} [samples];Events",
                        kDtBins, kDtMin, kDtMax);

  std::vector<double> dtA_vals;
  std::vector<double> dtB_vals;
  std::vector<double> EA_vals;
  std::vector<double> EB_vals;
  dtA_vals.reserve(20000);
  dtB_vals.reserve(6000);
  EA_vals.reserve(20000);
  EB_vals.reserve(6000);

  std::vector<double> buf_ps[2], buf_na[4];
  for (int i = 0; i < 2; ++i) buf_ps[i].reserve(kSamplesPerEvent);
  for (int i = 0; i < 4; ++i) buf_na[i].reserve(kSamplesPerEvent);

  long n_events = 0;
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

    std::vector<Pulse> ps_p[2];
    for (int i = 0; i < 2; ++i) {
      const double b = ModeBaselineInEvent(buf_ps[i]);
      std::vector<double> s;
      BuildS(buf_ps[i], b, s);
      ps_p[i] = FindPulses(s, kPSSeedThr, kPSEndThr, kPSMinSep, kPSMinWidth);
    }

    // PS earliest
    bool has_ps = false;
    int tps = -1;
    double aps = 0.0;
    for (int i = 0; i < 2; ++i) {
      for (const auto& p : ps_p[i]) {
        if (!has_ps || p.start < tps) {
          has_ps = true;
          tps = p.start;
          aps = p.amp;
        }
      }
    }

    if (has_ps && aps >= kPSAmpMinForDt) {
      std::vector<Pulse> na_p[4];
      for (int ch = 0; ch < 4; ++ch) {
        const double b = ModeBaselineInEvent(buf_na[ch]);
        std::vector<double> s;
        BuildS(buf_na[ch], b, s);
        na_p[ch] = FindPulses(s, kNaISeedThr, kNaIEndThr, kNaIMinSep, kNaIMinWidth);
      }

      std::vector<ModulePulse> modA = BuildModulePairs(na_p[0], na_p[1], kPairMatchMaxDt);
      std::vector<ModulePulse> modB = BuildModulePairs(na_p[2], na_p[3], kPairMatchMaxDt);

      const int ia = FindNearestModuleIndex(modA, tps);
      if (ia >= 0) {
        const double dt = (double)(modA[ia].t - tps);
        hDtA->Fill(dt);
        dtA_vals.push_back(dt);
        EA_vals.push_back(modA[ia].area_sum);
      }
      const int ib = FindNearestModuleIndex(modB, tps);
      if (ib >= 0) {
        const double dt = (double)(modB[ib].t - tps);
        hDtB->Fill(dt);
        dtB_vals.push_back(dt);
        EB_vals.push_back(modB[ib].area_sum);
      }
    }

    for (int i = 0; i < 2; ++i) buf_ps[i].clear();
    for (int i = 0; i < 4; ++i) buf_na[i].clear();
  }

  const double dtA_prompt = FindPeakXInRange(hDtA, kPromptSearchMin, kPromptSearchMax);
  const double dtB_prompt = FindPeakXInRange(hDtB, kPromptSearchMin, kPromptSearchMax);
  const double dtA_delay  = FindPeakXInRange(hDtA, kDelayedSearchMin, kDelayedSearchMax);
  const double dtB_delay  = FindPeakXInRange(hDtB, kDelayedSearchMin, kDelayedSearchMax);

  TH1D* hEA_all = new TH1D("hEA_all", "ModuleA energy (all nearest pair);I_{A1}+I_{A2} [ADC*samples];Pairs",
                           kEBins, kEMin, kEMax);
  TH1D* hEB_all = new TH1D("hEB_all", "ModuleB energy (all nearest pair);I_{B1}+I_{B2} [ADC*samples];Pairs",
                           kEBins, kEMin, kEMax);
  TH1D* hEA_prompt = new TH1D("hEA_prompt",
                              Form("ModuleA energy (prompt dt=%.1f#pm%d);I_{A1}+I_{A2} [ADC*samples];Pairs",
                                   dtA_prompt, kBandHalfWindow),
                              kEBins, kEMin, kEMax);
  TH1D* hEA_delay = new TH1D("hEA_delay",
                             Form("ModuleA energy (delayed dt=%.1f#pm%d);I_{A1}+I_{A2} [ADC*samples];Pairs",
                                  dtA_delay, kBandHalfWindow),
                             kEBins, kEMin, kEMax);
  TH1D* hEA_d170 = new TH1D("hEA_d170",
                            Form("ModuleA energy (manual delayed dt=%.0f#pm%d);I_{A1}+I_{A2} [ADC*samples];Pairs",
                                 kManualDelayedCenter, kBandHalfWindow),
                            kEBins, kEMin, kEMax);
  TH1D* hEB_prompt = new TH1D("hEB_prompt",
                              Form("ModuleB energy (prompt dt=%.1f#pm%d);I_{B1}+I_{B2} [ADC*samples];Pairs",
                                   dtB_prompt, kBandHalfWindow),
                              kEBins, kEMin, kEMax);
  TH1D* hEB_delay = new TH1D("hEB_delay",
                             Form("ModuleB energy (delayed dt=%.1f#pm%d);I_{B1}+I_{B2} [ADC*samples];Pairs",
                                  dtB_delay, kBandHalfWindow),
                             kEBins, kEMin, kEMax);
  TH1D* hEB_d170 = new TH1D("hEB_d170",
                            Form("ModuleB energy (manual delayed dt=%.0f#pm%d);I_{B1}+I_{B2} [ADC*samples];Pairs",
                                 kManualDelayedCenter, kBandHalfWindow),
                            kEBins, kEMin, kEMax);

  for (size_t i = 0; i < EA_vals.size(); ++i) {
    const double dt = dtA_vals[i];
    const double e = EA_vals[i];
    hEA_all->Fill(e);
    if (std::abs(dt - dtA_prompt) <= (double)kBandHalfWindow) hEA_prompt->Fill(e);
    if (std::abs(dt - dtA_delay) <= (double)kBandHalfWindow) hEA_delay->Fill(e);
    if (std::abs(dt - kManualDelayedCenter) <= (double)kBandHalfWindow) hEA_d170->Fill(e);
  }
  for (size_t i = 0; i < EB_vals.size(); ++i) {
    const double dt = dtB_vals[i];
    const double e = EB_vals[i];
    hEB_all->Fill(e);
    if (std::abs(dt - dtB_prompt) <= (double)kBandHalfWindow) hEB_prompt->Fill(e);
    if (std::abs(dt - dtB_delay) <= (double)kBandHalfWindow) hEB_delay->Fill(e);
    if (std::abs(dt - kManualDelayedCenter) <= (double)kBandHalfWindow) hEB_d170->Fill(e);
  }

  std::cout << "[dtband] events = " << n_events << std::endl;
  std::cout << "[dtband] ModuleA prompt peak dt = " << dtA_prompt
            << " samp (" << dtA_prompt * kDtNs << " ns)" << std::endl;
  std::cout << "[dtband] ModuleA delayed peak dt = " << dtA_delay
            << " samp (" << dtA_delay * kDtNs << " ns)" << std::endl;
  std::cout << "[dtband] ModuleB prompt peak dt = " << dtB_prompt
            << " samp (" << dtB_prompt * kDtNs << " ns)" << std::endl;
  std::cout << "[dtband] ModuleB delayed peak dt = " << dtB_delay
            << " samp (" << dtB_delay * kDtNs << " ns)" << std::endl;
  std::cout << "[dtband] Entries A(all/prompt/delay) = "
            << (long)hEA_all->GetEntries() << " / "
            << (long)hEA_prompt->GetEntries() << " / "
            << (long)hEA_delay->GetEntries() << std::endl;
  std::cout << "[dtband] Entries A(manual dt~170) = "
            << (long)hEA_d170->GetEntries() << std::endl;
  std::cout << "[dtband] Entries B(all/prompt/delay) = "
            << (long)hEB_all->GetEntries() << " / "
            << (long)hEB_prompt->GetEntries() << " / "
            << (long)hEB_delay->GetEntries() << std::endl;
  std::cout << "[dtband] Entries B(manual dt~170) = "
            << (long)hEB_d170->GetEntries() << std::endl;

  TCanvas* c = new TCanvas("c_dtband", "dt band energy", kCanvasW, kCanvasH);
  c->Print(out_pdf + "[");

  c->Clear();
  TPaveText* pt = new TPaveText(0.05, 0.05, 0.95, 0.95, "NDC");
  pt->SetFillStyle(0);
  pt->SetBorderSize(0);
  pt->SetTextAlign(12);
  pt->SetTextSize(0.028);
  pt->AddText("=== Michel dt-band energy study ===");
  pt->AddText(Form("InputDir: %s", kInputDir));
  pt->AddText(Form("Tag: %s", kInputTag));
  pt->AddText(Form("PS ref = earliest, PS amp >= %.1f", kPSAmpMinForDt));
  pt->AddText(Form("Prompt band: dt peak in [%d,%d], delayed band: peak in [%d,%d], half-width=+-%d",
                   kPromptSearchMin, kPromptSearchMax,
                   kDelayedSearchMin, kDelayedSearchMax,
                   kBandHalfWindow));
  pt->AddText(Form("ModuleA prompt/delayed dt = %.1f / %.1f samples", dtA_prompt, dtA_delay));
  pt->AddText(Form("ModuleB prompt/delayed dt = %.1f / %.1f samples", dtB_prompt, dtB_delay));
  pt->AddText(Form("Entries A(all/prompt/delayed) = %ld / %ld / %ld",
                   (long)hEA_all->GetEntries(), (long)hEA_prompt->GetEntries(), (long)hEA_delay->GetEntries()));
  pt->AddText(Form("Entries A(manual delayed dt=%.0f#pm%d) = %ld",
                   kManualDelayedCenter, kBandHalfWindow, (long)hEA_d170->GetEntries()));
  pt->AddText(Form("Entries B(all/prompt/delayed) = %ld / %ld / %ld",
                   (long)hEB_all->GetEntries(), (long)hEB_prompt->GetEntries(), (long)hEB_delay->GetEntries()));
  pt->AddText(Form("Entries B(manual delayed dt=%.0f#pm%d) = %ld",
                   kManualDelayedCenter, kBandHalfWindow, (long)hEB_d170->GetEntries()));
  pt->Draw();
  c->Print(out_pdf);
  delete pt;

  c->Clear();
  c->Divide(2, 1);
  c->cd(1);
  hDtA->SetLineColor(kBlue + 1);
  hDtA->SetLineWidth(2);
  hDtA->Draw("hist");
  {
    const double ymax = hDtA->GetMaximum();
    TLine* lp = new TLine(dtA_prompt, 0.0, dtA_prompt, ymax);
    TLine* ld = new TLine(dtA_delay, 0.0, dtA_delay, ymax);
    lp->SetLineColor(kRed + 1); lp->SetLineStyle(2); lp->Draw();
    ld->SetLineColor(kMagenta + 2); ld->SetLineStyle(2); ld->Draw();
  }
  c->cd(2);
  hDtB->SetLineColor(kBlue + 1);
  hDtB->SetLineWidth(2);
  hDtB->Draw("hist");
  {
    const double ymax = hDtB->GetMaximum();
    TLine* lp = new TLine(dtB_prompt, 0.0, dtB_prompt, ymax);
    TLine* ld = new TLine(dtB_delay, 0.0, dtB_delay, ymax);
    lp->SetLineColor(kRed + 1); lp->SetLineStyle(2); lp->Draw();
    ld->SetLineColor(kMagenta + 2); ld->SetLineStyle(2); ld->Draw();
  }
  c->Print(out_pdf);

  c->Clear();
  c->Divide(2, 1);
  c->cd(1);
  if (kLogY) gPad->SetLogy();
  hEA_all->SetLineColor(kGray + 2); hEA_all->SetLineWidth(2); hEA_all->Draw("hist");
  hEA_prompt->SetLineColor(kRed + 1); hEA_prompt->SetLineWidth(2); hEA_prompt->Draw("hist same");
  hEA_delay->SetLineColor(kMagenta + 2); hEA_delay->SetLineWidth(2); hEA_delay->Draw("hist same");
  hEA_d170->SetLineColor(kGreen + 2); hEA_d170->SetLineWidth(2); hEA_d170->Draw("hist same");
  {
    TLegend* lg = new TLegend(0.58, 0.70, 0.92, 0.90);
    lg->SetBorderSize(0);
    lg->SetFillStyle(0);
    lg->AddEntry(hEA_all, "A all nearest", "l");
    lg->AddEntry(hEA_prompt, "A prompt band", "l");
    lg->AddEntry(hEA_delay, "A delayed band", "l");
    lg->AddEntry(hEA_d170, "A dt=170#pm20", "l");
    lg->Draw();
  }
  c->cd(2);
  if (kLogY) gPad->SetLogy();
  hEB_all->SetLineColor(kGray + 2); hEB_all->SetLineWidth(2); hEB_all->Draw("hist");
  hEB_prompt->SetLineColor(kRed + 1); hEB_prompt->SetLineWidth(2); hEB_prompt->Draw("hist same");
  hEB_delay->SetLineColor(kMagenta + 2); hEB_delay->SetLineWidth(2); hEB_delay->Draw("hist same");
  hEB_d170->SetLineColor(kGreen + 2); hEB_d170->SetLineWidth(2); hEB_d170->Draw("hist same");
  {
    TLegend* lg = new TLegend(0.58, 0.70, 0.92, 0.90);
    lg->SetBorderSize(0);
    lg->SetFillStyle(0);
    lg->AddEntry(hEB_all, "B all nearest", "l");
    lg->AddEntry(hEB_prompt, "B prompt band", "l");
    lg->AddEntry(hEB_delay, "B delayed band", "l");
    lg->AddEntry(hEB_d170, "B dt=170#pm20", "l");
    lg->Draw();
  }
  c->Print(out_pdf);

  // 形比較（正規化）
  c->Clear();
  c->Divide(2, 1);
  c->cd(1);
  if (kLogY) gPad->SetLogy();
  TH1D* hEA_prompt_n = (TH1D*)hEA_prompt->Clone("hEA_prompt_n");
  TH1D* hEA_delay_n = (TH1D*)hEA_delay->Clone("hEA_delay_n");
  TH1D* hEA_d170_n = (TH1D*)hEA_d170->Clone("hEA_d170_n");
  if (hEA_prompt_n->Integral() > 0) hEA_prompt_n->Scale(1.0 / hEA_prompt_n->Integral());
  if (hEA_delay_n->Integral() > 0) hEA_delay_n->Scale(1.0 / hEA_delay_n->Integral());
  if (hEA_d170_n->Integral() > 0) hEA_d170_n->Scale(1.0 / hEA_d170_n->Integral());
  hEA_prompt_n->SetLineColor(kRed + 1); hEA_prompt_n->SetLineWidth(2); hEA_prompt_n->Draw("hist");
  hEA_delay_n->SetLineColor(kMagenta + 2); hEA_delay_n->SetLineWidth(2); hEA_delay_n->Draw("hist same");
  hEA_d170_n->SetLineColor(kGreen + 2); hEA_d170_n->SetLineWidth(2); hEA_d170_n->Draw("hist same");
  hEA_prompt_n->SetTitle("ModuleA normalized shape;I_{A1}+I_{A2} [ADC*samples];a.u.");
  {
    TLegend* lg = new TLegend(0.58, 0.74, 0.92, 0.90);
    lg->SetBorderSize(0); lg->SetFillStyle(0);
    lg->AddEntry(hEA_prompt_n, "A prompt norm", "l");
    lg->AddEntry(hEA_delay_n, "A delayed norm", "l");
    lg->AddEntry(hEA_d170_n, "A dt170 norm", "l");
    lg->Draw();
  }
  c->cd(2);
  if (kLogY) gPad->SetLogy();
  TH1D* hEB_prompt_n = (TH1D*)hEB_prompt->Clone("hEB_prompt_n");
  TH1D* hEB_delay_n = (TH1D*)hEB_delay->Clone("hEB_delay_n");
  TH1D* hEB_d170_n = (TH1D*)hEB_d170->Clone("hEB_d170_n");
  if (hEB_prompt_n->Integral() > 0) hEB_prompt_n->Scale(1.0 / hEB_prompt_n->Integral());
  if (hEB_delay_n->Integral() > 0) hEB_delay_n->Scale(1.0 / hEB_delay_n->Integral());
  if (hEB_d170_n->Integral() > 0) hEB_d170_n->Scale(1.0 / hEB_d170_n->Integral());
  hEB_prompt_n->SetLineColor(kRed + 1); hEB_prompt_n->SetLineWidth(2); hEB_prompt_n->Draw("hist");
  hEB_delay_n->SetLineColor(kMagenta + 2); hEB_delay_n->SetLineWidth(2); hEB_delay_n->Draw("hist same");
  hEB_d170_n->SetLineColor(kGreen + 2); hEB_d170_n->SetLineWidth(2); hEB_d170_n->Draw("hist same");
  hEB_prompt_n->SetTitle("ModuleB normalized shape;I_{B1}+I_{B2} [ADC*samples];a.u.");
  {
    TLegend* lg = new TLegend(0.58, 0.74, 0.92, 0.90);
    lg->SetBorderSize(0); lg->SetFillStyle(0);
    lg->AddEntry(hEB_prompt_n, "B prompt norm", "l");
    lg->AddEntry(hEB_delay_n, "B delayed norm", "l");
    lg->AddEntry(hEB_d170_n, "B dt170 norm", "l");
    lg->Draw();
  }
  c->Print(out_pdf);

  c->Print(out_pdf + "]");
  delete c;

  TFile fout(out_root, "RECREATE");
  hDtA->Write();
  hDtB->Write();
  hEA_all->Write();
  hEA_prompt->Write();
  hEA_delay->Write();
  hEA_d170->Write();
  hEB_all->Write();
  hEB_prompt->Write();
  hEB_delay->Write();
  hEB_d170->Write();
  fout.Close();

  std::cout << "[dtband] output pdf  : " << out_pdf << std::endl;
  std::cout << "[dtband] output root : " << out_root << std::endl;
}
