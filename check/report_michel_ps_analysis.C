// check/report_michel_ps_analysis.C
//
// 目的:
//   study_michel_ps_tagged.C の結果を、再現可能な形で定量化し、
//   レポート用の図・数値サマリを出力する。
//
// 出力:
//   doc/mainexp/report_michel_ps_analysis_run8.pdf
//   doc/mainexp/report_michel_ps_analysis_run8.root
//   doc/mainexp/report_michel_ps_analysis_run8.txt
//
// 実行:
//   root -l -b -q 'check/report_michel_ps_analysis.C'

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TStyle.h"
#include "TSystem.h"

// 既存の解析ロジック（PS閾値、NaI分解、moduleペア化）をそのまま使う。
#include "study_michel_ps_tagged.C"

struct ShapeStats {
  double entries = 0.0;
  double n_low = 0.0;   // [0, 20k]
  double n_mid = 0.0;   // [40k, 120k]
  double n_high = 0.0;  // [120k, 200k]
  double f_low = 0.0;
  double f_mid = 0.0;
  double f_high = 0.0;
  double first_nonzero_x = -1.0;
};

static double IntegralInRange(const TH1D* h, double xmin, double xmax)
{
  if (!h) return 0.0;
  int b0 = h->GetXaxis()->FindBin(xmin);
  int b1 = h->GetXaxis()->FindBin(xmax);
  if (b1 < b0) std::swap(b0, b1);
  return h->Integral(b0, b1);
}

static double FindFirstNonzeroX(const TH1D* h)
{
  if (!h) return -1.0;
  for (int b = 1; b <= h->GetNbinsX(); ++b) {
    if (h->GetBinContent(b) > 0.0) return h->GetBinLowEdge(b);
  }
  return -1.0;
}

static ShapeStats ComputeShapeStats(const TH1D* h)
{
  ShapeStats s;
  if (!h) return s;

  s.entries = h->GetEntries();
  s.n_low = IntegralInRange(h, 0.0, 2.0e4);
  s.n_mid = IntegralInRange(h, 4.0e4, 1.2e5);
  s.n_high = IntegralInRange(h, 1.2e5, 2.0e5);
  s.first_nonzero_x = FindFirstNonzeroX(h);

  const double n_all = IntegralInRange(h, 0.0, 2.0e5);
  if (n_all > 0.0) {
    s.f_low = s.n_low / n_all;
    s.f_mid = s.n_mid / n_all;
    s.f_high = s.n_high / n_all;
  }
  return s;
}

static double PeakXInRange(const TH1D* h, double xmin, double xmax)
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

static bool FindBestDtToAnyPS(const std::vector<ModulePulse>& mod,
                              const std::vector<double>& tps_all,
                              double& dt_best)
{
  dt_best = 0.0;
  if (mod.empty() || tps_all.empty()) return false;

  bool found = false;
  double dmin = 1.0e30;
  for (const auto& m : mod) {
    for (double tps : tps_all) {
      const double d = std::abs(m.t - tps);
      if (!found || d < dmin) {
        found = true;
        dmin = d;
        dt_best = m.t - tps;
      }
    }
  }
  return found;
}

static double MinAbsDtToPS(double t_mod, const std::vector<double>& tps_all)
{
  if (tps_all.empty()) return 1.0e30;
  double best = 1.0e30;
  for (double tps : tps_all) {
    const double d = t_mod - tps;
    if (std::abs(d) < std::abs(best)) best = d;
  }
  return best;
}

static void ScaleToUnity(TH1D* h)
{
  if (!h) return;
  const double s = h->Integral();
  if (s > 0.0) h->Scale(1.0 / s);
}

static void PrintShapeLine(std::ofstream& ofs,
                           const char* label,
                           const ShapeStats& s)
{
  ofs << std::left << std::setw(34) << label
      << " entries=" << std::setw(8) << (long)s.entries
      << " first_x=" << std::setw(8) << std::fixed << std::setprecision(0) << s.first_nonzero_x
      << " low=" << std::setw(8) << (long)s.n_low
      << " mid=" << std::setw(8) << (long)s.n_mid
      << " high=" << std::setw(8) << (long)s.n_high
      << " f_low=" << std::setprecision(4) << s.f_low
      << " f_mid=" << s.f_mid
      << " f_high=" << s.f_high
      << "\n";
}

void report_michel_ps_analysis()
{
  gStyle->SetOptStat(1110);
  gSystem->mkdir("doc/mainexp", true);

  const TString out_pdf = "doc/mainexp/report_michel_ps_analysis_run8.pdf";
  const TString out_root = "doc/mainexp/report_michel_ps_analysis_run8.root";
  const TString out_txt = "doc/mainexp/report_michel_ps_analysis_run8.txt";

  // ============================================================
  // Pass 1: PS閾値を study_michel_ps_tagged.C と同じ手順で自動決定
  // ============================================================
  TH1D* hPSAmp[2] = {
    new TH1D("hPSAmp_A_report", "PS_A amp;counts;Pulses", kPSAmpBins, kPSAmpMin, kPSAmpMax),
    new TH1D("hPSAmp_B_report", "PS_B amp;counts;Pulses", kPSAmpBins, kPSAmpMin, kPSAmpMax)
  };

  {
    TString fps[2];
    for (int i = 0; i < 2; ++i) fps[i] = Form("%s/%s", kInputDir, kPSFiles[i]);

    std::ifstream in_ps[2];
    for (int i = 0; i < 2; ++i) {
      in_ps[i].open(fps[i].Data());
      if (!in_ps[i]) {
        std::cerr << "[report] ERROR: cannot open " << fps[i] << std::endl;
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

      for (int ip = 0; ip < 2; ++ip) {
        const double b = ModeBaselineInEvent(buf_ps[ip]);
        std::vector<double> s;
        BuildS(buf_ps[ip], b, s);
        const std::vector<Pulse> ps = FindPulsesSimple(s, kPSSeedThr, kPSEndThr, kPSMinSep, kPSMinWidth);
        for (const auto& q : ps) {
          if (q.amp > 0.0) hPSAmp[ip]->Fill(q.amp);
        }
      }

      for (int i = 0; i < 2; ++i) buf_ps[i].clear();
    }
  }

  int ps_peak_bin[2] = {-1, -1};
  int ps_valley_bin[2] = {-1, -1};
  double ps_thr_auto[2] = {
    FindPSValleyThreshold(hPSAmp[0], ps_peak_bin[0], ps_valley_bin[0]),
    FindPSValleyThreshold(hPSAmp[1], ps_peak_bin[1], ps_valley_bin[1])
  };

  // ============================================================
  // Pass 2: NaI再構成 + 条件別のヒスト作成
  // ============================================================

  TH1D* hDtA_earliest = new TH1D("hDtA_earliest", "A dt (earliest PS);dt [samples];Events", kDtBins, kDtMin, kDtMax);
  TH1D* hDtB_earliest = new TH1D("hDtB_earliest", "B dt (earliest PS);dt [samples];Events", kDtBins, kDtMin, kDtMax);
  TH1D* hDtA_anyps = new TH1D("hDtA_anyps", "A dt (best among any PS);dt [samples];Events", kDtBins, kDtMin, kDtMax);
  TH1D* hDtB_anyps = new TH1D("hDtB_anyps", "B dt (best among any PS);dt [samples];Events", kDtBins, kDtMin, kDtMax);
  TH1D* hPSGap = new TH1D("hPSGap", "PS selected pulses gap (2nd-1st);#Delta t [samples];Events", 600, -600, 600);

  TH1D* hA1_all = new TH1D("hA1_all", "A1 fit area (no PS);I [ADC*samples];Pulses", kEnergyBins, kEnergyMin, kEnergyMax);
  TH1D* hA2_all = new TH1D("hA2_all", "A2 fit area (no PS);I [ADC*samples];Pulses", kEnergyBins, kEnergyMin, kEnergyMax);
  TH1D* hAsum_allpairs_nops = new TH1D("hAsum_allpairs_nops", "A1+A2 all pairs (no PS);I_{A1}+I_{A2} [ADC*samples];Pairs", kEnergyBins, kEnergyMin, kEnergyMax);
  TH1D* hAsum_allpairs_ps = new TH1D("hAsum_allpairs_ps", "A1+A2 all pairs (PS-selected events);I_{A1}+I_{A2} [ADC*samples];Pairs", kEnergyBins, kEnergyMin, kEnergyMax);
  TH1D* hAsum_nearest_ps = new TH1D("hAsum_nearest_ps", "A1+A2 nearest pair to earliest PS;I_{A1}+I_{A2} [ADC*samples];Pairs", kEnergyBins, kEnergyMin, kEnergyMax);
  TH1D* hAsum_prompt_ps = new TH1D("hAsum_prompt_ps", "A1+A2 nearest + prompt dt;I_{A1}+I_{A2} [ADC*samples];Pairs", kEnergyBins, kEnergyMin, kEnergyMax);
  TH1D* hAsum_delay_ps = new TH1D("hAsum_delay_ps", "A1+A2 nearest + delayed dt;I_{A1}+I_{A2} [ADC*samples];Pairs", kEnergyBins, kEnergyMin, kEnergyMax);

  std::vector<double> dtA_vals;
  std::vector<double> dtB_vals;
  std::vector<double> EA_vals;
  dtA_vals.reserve(30000);
  dtB_vals.reserve(12000);
  EA_vals.reserve(30000);

  long n_events = 0;
  long n_ps_selected = 0;
  long nA_neg_band = 0;
  long nA_delay_band = 0;
  long nA_neg_reprompt_sel = 0;
  long nA_delay_reprompt_sel = 0;
  long nA_neg_reprompt_seed = 0;
  long nA_delay_reprompt_seed = 0;

  {
    TString fps[2], fna[4];
    for (int i = 0; i < 2; ++i) fps[i] = Form("%s/%s", kInputDir, kPSFiles[i]);
    for (int i = 0; i < 4; ++i) fna[i] = Form("%s/%s", kInputDir, kNaIFiles[i]);

    std::ifstream in_ps[2];
    std::ifstream in_na[4];

    for (int i = 0; i < 2; ++i) {
      in_ps[i].open(fps[i].Data());
      if (!in_ps[i]) {
        std::cerr << "[report] ERROR: cannot open " << fps[i] << std::endl;
        return;
      }
    }
    for (int i = 0; i < 4; ++i) {
      in_na[i].open(fna[i].Data());
      if (!in_na[i]) {
        std::cerr << "[report] ERROR: cannot open " << fna[i] << std::endl;
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

      // NaI分解（studyマクロと同一）
      std::vector<FitPulse> fp[4];
      for (int ch = 0; ch < 4; ++ch) {
        const double b = ModeBaselineInEvent(buf_na[ch]);
        std::vector<double> s;
        BuildS(buf_na[ch], b, s);
        fp[ch] = FitPulsesWithDecomposition(s, kTrFix[ch], kDFix[ch]);
      }

      for (const auto& q : fp[0]) hA1_all->Fill(q.I);
      for (const auto& q : fp[1]) hA2_all->Fill(q.I);

      const std::vector<ModulePulse> modA = BuildModulePairs(fp[0], fp[1], kPairMatchMaxDt);
      const std::vector<ModulePulse> modB = BuildModulePairs(fp[2], fp[3], kPairMatchMaxDt);
      for (const auto& m : modA) hAsum_allpairs_nops->Fill(m.area_sum);

      // PS選別パルス（最早だけでなく、全候補時刻も保持）
      bool has_ps = false;
      double tps_earliest = -1.0;
      std::vector<double> tps_all_sel;
      std::vector<double> tps_all_seed;

      for (int ip = 0; ip < 2; ++ip) {
        const double b = ModeBaselineInEvent(buf_ps[ip]);
        std::vector<double> s;
        BuildS(buf_ps[ip], b, s);
        const std::vector<Pulse> ps = FindPulsesSimple(s, kPSSeedThr, kPSEndThr, kPSMinSep, kPSMinWidth);

        for (const auto& q : ps) {
          const double t = (double)q.start;
          tps_all_seed.push_back(t);
          if (q.amp < ps_thr_auto[ip]) continue;
          tps_all_sel.push_back(t);
          if (!has_ps || t < tps_earliest) {
            has_ps = true;
            tps_earliest = t;
          }
        }
      }

      if (has_ps) {
        n_ps_selected++;

        std::sort(tps_all_sel.begin(), tps_all_sel.end());
        if (tps_all_sel.size() >= 2) hPSGap->Fill(tps_all_sel[1] - tps_all_sel[0]);

        for (const auto& m : modA) hAsum_allpairs_ps->Fill(m.area_sum);

        // studyマクロと同じ earliest-PS基準
        const int ia = FindNearestModuleIndex(modA, tps_earliest);
        if (ia >= 0) {
          const double t_mod = modA[ia].t;
          const double dt = t_mod - tps_earliest;
          hDtA_earliest->Fill(dt);
          dtA_vals.push_back(dt);
          EA_vals.push_back(modA[ia].area_sum);
          hAsum_nearest_ps->Fill(modA[ia].area_sum);

          if (-220.0 <= dt && dt <= -100.0) {
            nA_neg_band++;
            const double dt_best_sel = MinAbsDtToPS(t_mod, tps_all_sel);
            const double dt_best_seed = MinAbsDtToPS(t_mod, tps_all_seed);
            if (std::abs(dt_best_sel) <= 20.0) nA_neg_reprompt_sel++;
            if (std::abs(dt_best_seed) <= 20.0) nA_neg_reprompt_seed++;
          }
          if (110.0 <= dt && dt <= 170.0) {
            nA_delay_band++;
            const double dt_best_sel = MinAbsDtToPS(t_mod, tps_all_sel);
            const double dt_best_seed = MinAbsDtToPS(t_mod, tps_all_seed);
            if (std::abs(dt_best_sel) <= 20.0) nA_delay_reprompt_sel++;
            if (std::abs(dt_best_seed) <= 20.0) nA_delay_reprompt_seed++;
          }
        }

        const int ib = FindNearestModuleIndex(modB, tps_earliest);
        if (ib >= 0) {
          const double dt = modB[ib].t - tps_earliest;
          hDtB_earliest->Fill(dt);
          dtB_vals.push_back(dt);
        }

        // 参考: event内の任意PS時刻との最短差
        double dt_any = 0.0;
        if (FindBestDtToAnyPS(modA, tps_all_sel, dt_any)) hDtA_anyps->Fill(dt_any);
        if (FindBestDtToAnyPS(modB, tps_all_sel, dt_any)) hDtB_anyps->Fill(dt_any);
      }

      for (int i = 0; i < 2; ++i) buf_ps[i].clear();
      for (int i = 0; i < 4; ++i) buf_na[i].clear();
    }
  }

  // dtピーク（studyと同じ探索窓）
  const double dtA_prompt = PeakXInRange(hDtA_earliest, kPromptSearchMin, kPromptSearchMax);
  const double dtA_delay = PeakXInRange(hDtA_earliest, kDelayedSearchMin, kDelayedSearchMax);
  const double dtB_prompt = PeakXInRange(hDtB_earliest, kPromptSearchMin, kPromptSearchMax);
  const double dtB_delay = PeakXInRange(hDtB_earliest, kDelayedSearchMin, kDelayedSearchMax);

  // PS間隔ピーク（参考）
  const double ps_gap_peak = PeakXInRange(hPSGap, 80.0, 260.0);

  for (size_t i = 0; i < EA_vals.size(); ++i) {
    const double dt = dtA_vals[i];
    if (std::abs(dt - dtA_prompt) <= (double)kDtBandHalfWindow) hAsum_prompt_ps->Fill(EA_vals[i]);
    if (std::abs(dt - dtA_delay) <= (double)kDtBandHalfWindow) hAsum_delay_ps->Fill(EA_vals[i]);
  }

  const ShapeStats st_a1 = ComputeShapeStats(hA1_all);
  const ShapeStats st_a2 = ComputeShapeStats(hA2_all);
  const ShapeStats st_nops = ComputeShapeStats(hAsum_allpairs_nops);
  const ShapeStats st_psall = ComputeShapeStats(hAsum_allpairs_ps);
  const ShapeStats st_near = ComputeShapeStats(hAsum_nearest_ps);
  const ShapeStats st_prompt = ComputeShapeStats(hAsum_prompt_ps);
  const ShapeStats st_delay = ComputeShapeStats(hAsum_delay_ps);

  // ============================================================
  // 数値サマリを書き出し
  // ============================================================
  {
    std::ofstream ofs(out_txt.Data());
    ofs << "# report_michel_ps_analysis_run8 summary\n";
    ofs << "input_dir = " << kInputDir << "\n";
    ofs << "events_total = " << n_events << "\n";
    ofs << "events_ps_selected = " << n_ps_selected << "\n";
    ofs << "ps_selection_fraction = "
        << (n_events > 0 ? (double)n_ps_selected / (double)n_events : 0.0) << "\n";
    ofs << "ps_threshold_A = " << ps_thr_auto[0] << "\n";
    ofs << "ps_threshold_B = " << ps_thr_auto[1] << "\n";
    ofs << "\n";

    ofs << "[dt_peaks_earliestPS]\n";
    ofs << "dtA_prompt = " << dtA_prompt << " samples\n";
    ofs << "dtA_delay  = " << dtA_delay << " samples\n";
    ofs << "dtB_prompt = " << dtB_prompt << " samples\n";
    ofs << "dtB_delay  = " << dtB_delay << " samples\n";
    ofs << "ps_gap_peak(80-260) = " << ps_gap_peak << " samples\n";
    ofs << "\n";

    ofs << "[dt_yields_earliestPS]\n";
    ofs << "A_neg(-220,-100)   = " << IntegralInRange(hDtA_earliest, -220.0, -100.0) << "\n";
    ofs << "A_prompt(-20,20)   = " << IntegralInRange(hDtA_earliest, -20.0, 20.0) << "\n";
    ofs << "A_delay(110,170)   = " << IntegralInRange(hDtA_earliest, 110.0, 170.0) << "\n";
    ofs << "B_neg(-220,-100)   = " << IntegralInRange(hDtB_earliest, -220.0, -100.0) << "\n";
    ofs << "B_prompt(-20,20)   = " << IntegralInRange(hDtB_earliest, -20.0, 20.0) << "\n";
    ofs << "B_delay(110,170)   = " << IntegralInRange(hDtB_earliest, 110.0, 170.0) << "\n";
    ofs << "\n";

    ofs << "[dt_yields_anyPS_reference]\n";
    ofs << "A_delay(110,170)   = " << IntegralInRange(hDtA_anyps, 110.0, 170.0) << "\n";
    ofs << "A_prompt(-20,20)   = " << IntegralInRange(hDtA_anyps, -20.0, 20.0) << "\n";
    ofs << "B_delay(110,170)   = " << IntegralInRange(hDtB_anyps, 110.0, 170.0) << "\n";
    ofs << "B_prompt(-20,20)   = " << IntegralInRange(hDtB_anyps, -20.0, 20.0) << "\n";
    ofs << "\n";

    ofs << "[dt_reassociation_A_earliest]\n";
    ofs << "A_neg_band_count = " << nA_neg_band << "\n";
    ofs << "A_neg_reprompt_with_selectedPS(|dt|<=20) = " << nA_neg_reprompt_sel
        << " (frac=" << (nA_neg_band > 0 ? (double)nA_neg_reprompt_sel / (double)nA_neg_band : 0.0) << ")\n";
    ofs << "A_neg_reprompt_with_allSeedPS(|dt|<=20) = " << nA_neg_reprompt_seed
        << " (frac=" << (nA_neg_band > 0 ? (double)nA_neg_reprompt_seed / (double)nA_neg_band : 0.0) << ")\n";
    ofs << "A_delay_band_count = " << nA_delay_band << "\n";
    ofs << "A_delay_reprompt_with_selectedPS(|dt|<=20) = " << nA_delay_reprompt_sel
        << " (frac=" << (nA_delay_band > 0 ? (double)nA_delay_reprompt_sel / (double)nA_delay_band : 0.0) << ")\n";
    ofs << "A_delay_reprompt_with_allSeedPS(|dt|<=20) = " << nA_delay_reprompt_seed
        << " (frac=" << (nA_delay_band > 0 ? (double)nA_delay_reprompt_seed / (double)nA_delay_band : 0.0) << ")\n";
    ofs << "\n";

    ofs << "[shape_stats]\n";
    PrintShapeLine(ofs, "A1 all fits (no PS)", st_a1);
    PrintShapeLine(ofs, "A2 all fits (no PS)", st_a2);
    PrintShapeLine(ofs, "A1+A2 all pairs (no PS)", st_nops);
    PrintShapeLine(ofs, "A1+A2 all pairs (PS event)", st_psall);
    PrintShapeLine(ofs, "A1+A2 nearest to PS", st_near);
    PrintShapeLine(ofs, "A1+A2 nearest + prompt dt", st_prompt);
    PrintShapeLine(ofs, "A1+A2 nearest + delayed dt", st_delay);
  }

  // ============================================================
  // 図を作成
  // ============================================================
  TCanvas* c = new TCanvas("c_report_michel_ps", "report michel ps", 1400, 1000);
  c->Print(out_pdf + "[");

  // page 1: dt（earliest PS）
  c->Clear();
  c->Divide(2, 1);
  c->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  hDtA_earliest->SetLineColor(kBlue + 1);
  hDtA_earliest->SetLineWidth(2);
  hDtA_earliest->Draw("hist");
  {
    const double ymax = hDtA_earliest->GetMaximum();
    TLine* l0 = new TLine(0.0, 0.0, 0.0, ymax);
    TLine* lp = new TLine(dtA_prompt, 0.0, dtA_prompt, ymax);
    TLine* ld = new TLine(dtA_delay, 0.0, dtA_delay, ymax);
    l0->SetLineColor(kGray + 2); l0->SetLineStyle(3); l0->Draw();
    lp->SetLineColor(kRed + 1); lp->SetLineStyle(2); lp->SetLineWidth(2); lp->Draw();
    ld->SetLineColor(kMagenta + 2); ld->SetLineStyle(2); ld->SetLineWidth(2); ld->Draw();
  }

  c->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  hDtB_earliest->SetLineColor(kBlue + 1);
  hDtB_earliest->SetLineWidth(2);
  hDtB_earliest->Draw("hist");
  {
    const double ymax = hDtB_earliest->GetMaximum();
    TLine* l0 = new TLine(0.0, 0.0, 0.0, ymax);
    TLine* lp = new TLine(dtB_prompt, 0.0, dtB_prompt, ymax);
    TLine* ld = new TLine(dtB_delay, 0.0, dtB_delay, ymax);
    l0->SetLineColor(kGray + 2); l0->SetLineStyle(3); l0->Draw();
    lp->SetLineColor(kRed + 1); lp->SetLineStyle(2); lp->SetLineWidth(2); lp->Draw();
    ld->SetLineColor(kMagenta + 2); ld->SetLineStyle(2); ld->SetLineWidth(2); ld->Draw();
  }
  c->Print(out_pdf);

  // page 2: dt 参照方式の比較
  c->Clear();
  c->Divide(2, 1);
  c->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  hDtA_earliest->SetLineColor(kBlue + 1);
  hDtA_anyps->SetLineColor(kRed + 1);
  hDtA_earliest->SetLineWidth(2);
  hDtA_anyps->SetLineWidth(2);
  hDtA_earliest->Draw("hist");
  hDtA_anyps->Draw("hist same");
  {
    TLegend* lg = new TLegend(0.53, 0.74, 0.92, 0.90);
    lg->SetBorderSize(0);
    lg->SetFillStyle(0);
    lg->AddEntry(hDtA_earliest, "A: earliest PS reference", "l");
    lg->AddEntry(hDtA_anyps, "A: best among any PS", "l");
    lg->Draw();
  }

  c->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  hDtB_earliest->SetLineColor(kBlue + 1);
  hDtB_anyps->SetLineColor(kRed + 1);
  hDtB_earliest->SetLineWidth(2);
  hDtB_anyps->SetLineWidth(2);
  hDtB_earliest->Draw("hist");
  hDtB_anyps->Draw("hist same");
  {
    TLegend* lg = new TLegend(0.53, 0.74, 0.92, 0.90);
    lg->SetBorderSize(0);
    lg->SetFillStyle(0);
    lg->AddEntry(hDtB_earliest, "B: earliest PS reference", "l");
    lg->AddEntry(hDtB_anyps, "B: best among any PS", "l");
    lg->Draw();
  }
  c->Print(out_pdf);

  // page 3: PS間隔
  c->Clear();
  gPad->SetGridx();
  gPad->SetGridy();
  hPSGap->SetLineColor(kBlue + 1);
  hPSGap->SetLineWidth(2);
  hPSGap->Draw("hist");
  {
    const double ymax = hPSGap->GetMaximum();
    TLine* lpk = new TLine(ps_gap_peak, 0.0, ps_gap_peak, ymax);
    lpk->SetLineColor(kRed + 1);
    lpk->SetLineStyle(2);
    lpk->SetLineWidth(2);
    lpk->Draw();
  }
  c->Print(out_pdf);

  // page 4: エネルギー分布（絶対）
  c->Clear();
  c->Divide(2, 1);
  c->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogy();
  hA1_all->SetLineColor(kGray + 2);
  hAsum_allpairs_nops->SetLineColor(kBlue + 1);
  hAsum_allpairs_ps->SetLineColor(kGreen + 2);
  hAsum_nearest_ps->SetLineColor(kRed + 1);
  hA1_all->SetLineWidth(2);
  hAsum_allpairs_nops->SetLineWidth(2);
  hAsum_allpairs_ps->SetLineWidth(2);
  hAsum_nearest_ps->SetLineWidth(2);
  hA1_all->SetTitle("A energy distributions (absolute);I [ADC*samples];Counts");
  hA1_all->Draw("hist");
  hAsum_allpairs_nops->Draw("hist same");
  hAsum_allpairs_ps->Draw("hist same");
  hAsum_nearest_ps->Draw("hist same");
  {
    TLegend* lg = new TLegend(0.45, 0.65, 0.92, 0.90);
    lg->SetBorderSize(0);
    lg->SetFillStyle(0);
    lg->AddEntry(hA1_all, "A1 all fits", "l");
    lg->AddEntry(hAsum_allpairs_nops, "A1+A2 all pairs (no PS)", "l");
    lg->AddEntry(hAsum_allpairs_ps, "A1+A2 all pairs (PS)", "l");
    lg->AddEntry(hAsum_nearest_ps, "A1+A2 nearest to PS", "l");
    lg->Draw();
  }

  c->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogy();
  TH1D* h_nops_norm = (TH1D*)hAsum_allpairs_nops->Clone("h_nops_norm");
  TH1D* h_psall_norm = (TH1D*)hAsum_allpairs_ps->Clone("h_psall_norm");
  TH1D* h_near_norm = (TH1D*)hAsum_nearest_ps->Clone("h_near_norm");
  TH1D* h_prompt_norm = (TH1D*)hAsum_prompt_ps->Clone("h_prompt_norm");
  ScaleToUnity(h_nops_norm);
  ScaleToUnity(h_psall_norm);
  ScaleToUnity(h_near_norm);
  ScaleToUnity(h_prompt_norm);
  h_nops_norm->SetLineColor(kBlue + 1);
  h_psall_norm->SetLineColor(kGreen + 2);
  h_near_norm->SetLineColor(kRed + 1);
  h_prompt_norm->SetLineColor(kMagenta + 2);
  h_nops_norm->SetLineWidth(2);
  h_psall_norm->SetLineWidth(2);
  h_near_norm->SetLineWidth(2);
  h_prompt_norm->SetLineWidth(2);
  h_nops_norm->SetTitle("A energy shape (normalized);I [ADC*samples];Arbitrary unit");
  h_nops_norm->Draw("hist");
  h_psall_norm->Draw("hist same");
  h_near_norm->Draw("hist same");
  h_prompt_norm->Draw("hist same");
  {
    TLegend* lg = new TLegend(0.42, 0.65, 0.92, 0.90);
    lg->SetBorderSize(0);
    lg->SetFillStyle(0);
    lg->AddEntry(h_nops_norm, "all pairs (no PS)", "l");
    lg->AddEntry(h_psall_norm, "all pairs (PS)", "l");
    lg->AddEntry(h_near_norm, "nearest to PS", "l");
    lg->AddEntry(h_prompt_norm, "nearest + prompt dt", "l");
    lg->Draw();
  }
  c->Print(out_pdf);

  c->Print(out_pdf + "]");
  delete c;

  // ROOTファイル保存
  {
    TFile fout(out_root.Data(), "RECREATE");
    hPSAmp[0]->Write();
    hPSAmp[1]->Write();
    hDtA_earliest->Write();
    hDtB_earliest->Write();
    hDtA_anyps->Write();
    hDtB_anyps->Write();
    hPSGap->Write();

    hA1_all->Write();
    hA2_all->Write();
    hAsum_allpairs_nops->Write();
    hAsum_allpairs_ps->Write();
    hAsum_nearest_ps->Write();
    hAsum_prompt_ps->Write();
    hAsum_delay_ps->Write();
    fout.Close();
  }

  std::cout << "[report] output pdf  : " << out_pdf << "\n";
  std::cout << "[report] output root : " << out_root << "\n";
  std::cout << "[report] output txt  : " << out_txt << "\n";
}
