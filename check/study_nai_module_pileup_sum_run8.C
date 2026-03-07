// check/study_nai_module_pileup_sum_run8.C
//
// 目的:
//   pileup をできるだけ拾う緩い NaI 判定でパルスを抽出し、
//   同一モジュール内（A1-A2, B1-B2）で
//   |t1 - t2| <= 2 samples の同時判定があれば
//   エネルギー和 (I1+I2) をヒストグラムへ追加する。
//
// 出力:
//   doc/mainexp/study_nai_module_pileup_sum_run8.pdf
//
// 実行:
//   root -l -b -q 'check/study_nai_module_pileup_sum_run8.C'

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"

// 既存の波形処理・フィット関数群を再利用する
#include "study_michel_ps_tagged.C"

// ============================================================
// pileup を拾いやすくするための緩い判定パラメータ
// ============================================================
static const double kNaISeedThrLoose = 35.0;
static const double kNaIEndThrLoose = 12.0;
static const int    kNaIMinSepLoose = 8;
static const int    kTriedSepLoose = 4;
static const int    kNearDupSepLoose = 4;
static const int    kMaxAcceptedPerChLoose = 14;
static const int    kMaxIterPerChLoose = 48;
static const double kChi2LooseMax = 250.0;
static const int    kCoincMaxDtSamples = 2; // 要求仕様

// ヒスト設定
static const int    kEBins = 2500;
static const double kEMin = 0.0;
static const double kEMax = 2.0e5;
static const int    kPileDtBins = 41;
static const double kPileDtMin = -20.5;
static const double kPileDtMax = 20.5;

static const int    kPileMulBins = 30;
static const double kPileMulMin = -0.5;
static const double kPileMulMax = 29.5;

// 局所ピーク由来候補の追加（seed 可変版）
static void AddPeakSupplementCandidatesLoose(const std::vector<double>& s,
                                             std::vector<int>& starts,
                                             double seed_thr)
{
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
    if (s[i] < seed_thr + 8.0) continue;

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
    if (!has_near(t0, 6)) starts.push_back(t0);
  }
}

// 緩い NaI 分解:
//   - seed/end/min_sep を緩める
//   - 近接重複の抑制を弱める
//   - chi2/ndf は緩い上限だけ掛ける
static std::vector<FitPulse> FitPulsesWithDecompositionLoose(const std::vector<double>& s_input,
                                                             double tr_fix,
                                                             double d_fix)
{
  std::vector<FitPulse> accepted;
  if (s_input.empty()) return accepted;

  std::vector<double> s_work = s_input;
  std::vector<int> tried;

  for (int iter = 0; iter < kMaxIterPerChLoose; ++iter) {
    if ((int)accepted.size() >= kMaxAcceptedPerChLoose) break;

    std::vector<int> starts = FindStartCandidates(s_work,
                                                  kNaISeedThrLoose,
                                                  kNaIEndThrLoose,
                                                  kNaIMinSepLoose);
    AddPeakSupplementCandidatesLoose(s_work, starts, kNaISeedThrLoose);
    if (starts.empty()) break;

    std::sort(starts.begin(), starts.end());
    starts.erase(std::unique(starts.begin(), starts.end()), starts.end());

    int best_t0 = -1;
    double best_sc = -1.0;
    for (int t0 : starts) {
      if (t0 - kFitPreWindow < 0) continue;
      if (t0 + kFitHalfWindow >= (int)s_work.size()) continue;
      if (HasNearbyT0(tried, t0, kTriedSepLoose)) continue;

      const double sc = CandidatePeakScore(s_work, t0);
      if (sc > best_sc) {
        best_sc = sc;
        best_t0 = t0;
      }
    }
    if (best_t0 < 0) break;

    tried.push_back(best_t0);

    // 既存 FitOnePulse は内部で標準 chi2 判定もするが、
    // ここでは A>0 と有限性を主条件にして緩く受理する。
    FitPulse fp = FitOnePulse(s_work, best_t0, tr_fix, d_fix);
    if (!(fp.A > 0.0)) continue;
    if (!std::isfinite(fp.t0_fit) || !std::isfinite(fp.I) || !std::isfinite(fp.chi2ndf)) continue;
    if (fp.chi2ndf > kChi2LooseMax) continue;

    bool near_dup = false;
    for (const auto& ap : accepted) {
      if (std::abs((int)std::lround(ap.t0_fit) - (int)std::lround(fp.t0_fit)) < kNearDupSepLoose) {
        near_dup = true;
        break;
      }
    }
    if (near_dup) continue;

    fp.ok = true;
    accepted.push_back(fp);
    SubtractAcceptedPulse(s_work, fp, tr_fix, d_fix);
  }

  std::sort(accepted.begin(), accepted.end(),
            [](const FitPulse& a, const FitPulse& b) { return a.t0_fit < b.t0_fit; });
  return accepted;
}

static void FillCoincidentPairSums(const std::vector<FitPulse>& c1,
                                   const std::vector<FitPulse>& c2,
                                   TH1D* hE,
                                   TH1D* hDt,
                                   long& n_pairs)
{
  if (!hE || !hDt) return;
  if (c1.empty() || c2.empty()) return;

  for (const auto& p1 : c1) {
    for (const auto& p2 : c2) {
      const double dt = p2.t0_fit - p1.t0_fit;
      if (std::abs(dt) > (double)kCoincMaxDtSamples) continue;
      hE->Fill(p1.I + p2.I);
      hDt->Fill(dt);
      n_pairs++;
    }
  }
}

void study_nai_module_pileup_sum_run8()
{
  gStyle->SetOptStat(1110);
  gSystem->mkdir(kOutputDir, true);

  const TString out_pdf = Form("%s/study_nai_module_pileup_sum_run8.pdf", kOutputDir);

  TH1D* hEA = new TH1D("hEA_moduleA",
                       "Module A pileup-coincident E_{sum} (|#Deltat|#le2);I_{A1}+I_{A2} [ADC*samples];Pairs",
                       kEBins, kEMin, kEMax);
  TH1D* hEB = new TH1D("hEB_moduleB",
                       "Module B pileup-coincident E_{sum} (|#Deltat|#le2);I_{B1}+I_{B2} [ADC*samples];Pairs",
                       kEBins, kEMin, kEMax);
  TH1D* hDtA = new TH1D("hDtA_moduleA",
                        "Module A #Deltat=t_{A2}-t_{A1} (coincident);#Deltat [samples];Pairs",
                        kPileDtBins, kPileDtMin, kPileDtMax);
  TH1D* hDtB = new TH1D("hDtB_moduleB",
                        "Module B #Deltat=t_{B2}-t_{B1} (coincident);#Deltat [samples];Pairs",
                        kPileDtBins, kPileDtMin, kPileDtMax);
  TH1D* hMulA1 = new TH1D("hMulA1", "A1 pulses/event (loose);N_{pulse};Events", kPileMulBins, kPileMulMin, kPileMulMax);
  TH1D* hMulA2 = new TH1D("hMulA2", "A2 pulses/event (loose);N_{pulse};Events", kPileMulBins, kPileMulMin, kPileMulMax);
  TH1D* hMulB1 = new TH1D("hMulB1", "B1 pulses/event (loose);N_{pulse};Events", kPileMulBins, kPileMulMin, kPileMulMax);
  TH1D* hMulB2 = new TH1D("hMulB2", "B2 pulses/event (loose);N_{pulse};Events", kPileMulBins, kPileMulMin, kPileMulMax);
  TH1D* hPairPerEvA = new TH1D("hPairPerEvA", "Module A coincident pairs/event;N_{pair};Events", kPileMulBins, kPileMulMin, kPileMulMax);
  TH1D* hPairPerEvB = new TH1D("hPairPerEvB", "Module B coincident pairs/event;N_{pair};Events", kPileMulBins, kPileMulMin, kPileMulMax);

  long n_events = 0;
  long n_pulses[4] = {0, 0, 0, 0};
  long n_pairs_A = 0;
  long n_pairs_B = 0;

  TString fna[4];
  for (int i = 0; i < 4; ++i) fna[i] = Form("%s/%s", kInputDir, kNaIFiles[i]);

  std::ifstream in_na[4];
  for (int i = 0; i < 4; ++i) {
    in_na[i].open(fna[i].Data());
    if (!in_na[i]) {
      std::cerr << "[pileup-sum] ERROR: cannot open " << fna[i] << std::endl;
      return;
    }
  }

  std::vector<double> buf_na[4];
  for (int i = 0; i < 4; ++i) buf_na[i].reserve(kSamplesPerEvent);

  while (true) {
    double xna[4];
    bool ok = true;
    for (int i = 0; i < 4; ++i) {
      if (!(in_na[i] >> xna[i])) { ok = false; break; }
    }
    if (!ok) break;

    for (int i = 0; i < 4; ++i) buf_na[i].push_back(xna[i]);
    if ((int)buf_na[0].size() < kSamplesPerEvent) continue;

    n_events++;

    std::vector<FitPulse> fp[4];
    for (int ch = 0; ch < 4; ++ch) {
      const double b = ModeBaselineInEvent(buf_na[ch]);
      std::vector<double> s;
      BuildS(buf_na[ch], b, s);
      fp[ch] = FitPulsesWithDecompositionLoose(s, kTrFix[ch], kDFix[ch]);
      n_pulses[ch] += (long)fp[ch].size();
    }

    hMulA1->Fill((double)fp[0].size());
    hMulA2->Fill((double)fp[1].size());
    hMulB1->Fill((double)fp[2].size());
    hMulB2->Fill((double)fp[3].size());

    const long pair_before_A = n_pairs_A;
    const long pair_before_B = n_pairs_B;
    FillCoincidentPairSums(fp[0], fp[1], hEA, hDtA, n_pairs_A);
    FillCoincidentPairSums(fp[2], fp[3], hEB, hDtB, n_pairs_B);
    hPairPerEvA->Fill((double)(n_pairs_A - pair_before_A));
    hPairPerEvB->Fill((double)(n_pairs_B - pair_before_B));

    for (int i = 0; i < 4; ++i) buf_na[i].clear();
  }

  for (int i = 0; i < 4; ++i) in_na[i].close();

  std::cout << "[pileup-sum] events = " << n_events << "\n";
  std::cout << "[pileup-sum] pulses A1/A2/B1/B2 = "
            << n_pulses[0] << " / "
            << n_pulses[1] << " / "
            << n_pulses[2] << " / "
            << n_pulses[3] << "\n";
  std::cout << "[pileup-sum] coincident pairs A/B (|dt|<=2) = "
            << n_pairs_A << " / " << n_pairs_B << "\n";

  // PDF 描画
  TCanvas* c = new TCanvas("c_pileup_sum", "NaI module pileup sum", 1400, 1000);
  c->Print(out_pdf + "[");

  c->Clear();
  c->Divide(2, 1);
  c->cd(1); gPad->SetLogy(1); hEA->SetLineColor(kBlue + 1); hEA->SetLineWidth(2); hEA->Draw("hist");
  c->cd(2); gPad->SetLogy(1); hEB->SetLineColor(kOrange + 7); hEB->SetLineWidth(2); hEB->Draw("hist");
  c->Print(out_pdf);

  c->Clear();
  c->Divide(2, 1);
  c->cd(1);
  hDtA->SetLineColor(kBlue + 1);
  hDtA->SetLineWidth(2);
  hDtA->Draw("hist");
  c->cd(2);
  hDtB->SetLineColor(kOrange + 7);
  hDtB->SetLineWidth(2);
  hDtB->Draw("hist");
  c->Print(out_pdf);

  c->Clear();
  c->Divide(2, 2);
  c->cd(1); hMulA1->SetLineColor(kBlue + 1); hMulA1->Draw("hist");
  c->cd(2); hMulA2->SetLineColor(kBlue + 2); hMulA2->Draw("hist");
  c->cd(3); hMulB1->SetLineColor(kOrange + 7); hMulB1->Draw("hist");
  c->cd(4); hMulB2->SetLineColor(kOrange + 9); hMulB2->Draw("hist");
  c->Print(out_pdf);

  c->Clear();
  c->Divide(2, 1);
  c->cd(1); hPairPerEvA->SetLineColor(kBlue + 1); hPairPerEvA->Draw("hist");
  c->cd(2); hPairPerEvB->SetLineColor(kOrange + 7); hPairPerEvB->Draw("hist");
  c->Print(out_pdf);

  c->Clear();
  TPaveText* pt = new TPaveText(0.06, 0.06, 0.94, 0.94, "NDC");
  pt->SetFillStyle(0);
  pt->SetBorderSize(0);
  pt->SetTextAlign(12);
  pt->SetTextSize(0.030);
  pt->AddText("=== NaI module pileup-coincident sum (run8) ===");
  pt->AddText(Form("InputDir: %s", kInputDir));
  pt->AddText("Loose NaI decision (pileup pickup priority):");
  pt->AddText(Form("  seed/end/min_sep = %.1f / %.1f / %d samples",
                   kNaISeedThrLoose, kNaIEndThrLoose, kNaIMinSepLoose));
  pt->AddText(Form("  near-dup sep = %d samples, chi2/ndf <= %.1f",
                   kNearDupSepLoose, kChi2LooseMax));
  pt->AddText(Form("Coincidence condition in one module: |t1 - t2| <= %d samples",
                   kCoincMaxDtSamples));
  pt->AddText(Form("Processed events = %ld", n_events));
  pt->AddText(Form("Pulses A1/A2/B1/B2 = %ld / %ld / %ld / %ld",
                   n_pulses[0], n_pulses[1], n_pulses[2], n_pulses[3]));
  pt->AddText(Form("Coincident pairs A/B = %ld / %ld",
                   n_pairs_A, n_pairs_B));
  pt->AddText(Form("Output PDF : %s", out_pdf.Data()));
  pt->Draw();
  c->Print(out_pdf);

  c->Print(out_pdf + "]");

  std::cout << "[pileup-sum] output pdf  : " << out_pdf << "\n";
}
