// check/reassess_michel_selection_run8.C
//
// 目的:
//   これまでのMichel選別考察の問題点（過学習、偶発見積の粗さ、A/B非対称）を
//   追加検証し、より頑健な結論を得る。
//
// 実施内容:
//   1) even/odd splitで最適化と再現性評価（疑似cross validation）
//   2) same-event と event-mixing の dt 比較で偶発成分を見積
//   3) 最終選別での prompt excess を定量化
//
// 出力:
//   doc/mainexp/reassess_michel_selection_run8.pdf
//   （数値は標準出力にのみ表示。doc 配下に txt/root は保存しない）
//
// 実行:
//   root -l -b -q 'check/reassess_michel_selection_run8.C'

#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

#include "TCanvas.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TSystem.h"

#include "study_michel_ps_tagged.C"

struct SideEventData {
  std::vector<double> ps_t;     // PS候補時刻（seed条件を満たす全候補）
  std::vector<double> ps_amp;   // 振幅
  std::vector<ModulePulse> mod; // NaI module パルス
};

struct SideCandidate {
  bool ok = false;
  double t_ps = 0.0;
  double t_mod = 0.0;
  double dt = 0.0;
  double ps_amp = 0.0;
  double e_sum = 0.0;
  int n_ps_sel = 0;
  int n_mod = 0;
};

struct CutDef {
  double scale_A = 1.0;
  double scale_B = 1.0;
  double dtw = 20.0;
  double e_min = 0.0;
  bool require_single_ps = false;
  bool require_single_mod = false;
};

struct EvalMetrics {
  long n_total = 0;
  long n_prompt = 0;
  long n_prompt_A = 0;
  long n_prompt_B = 0;
  long n_delay_pos = 0;
  long n_delay_neg = 0;
  long n_sb1 = 0; // 80<=|dt|<=160
  long n_sb2 = 0; // 200<=|dt|<=280
  double purity_sb1 = 0.0;
  double purity_sb2 = 0.0;
  double prompt_over_delay = 0.0;
  double side_balance = 0.0;
  double score = -1.0e30;
};

static bool BuildCandidate(const SideEventData& ev,
                           double ps_thr,
                           bool require_single_ps,
                           bool require_single_mod,
                           SideCandidate& out)
{
  out = SideCandidate{};
  if (ev.ps_t.size() != ev.ps_amp.size()) return false;
  if (ev.ps_t.empty()) return false;
  if (ev.mod.empty()) return false;

  std::vector<int> isel;
  isel.reserve(ev.ps_t.size());
  for (size_t i = 0; i < ev.ps_t.size(); ++i) {
    if (ev.ps_amp[i] >= ps_thr) isel.push_back((int)i);
  }
  if (isel.empty()) return false;
  if (require_single_ps && (int)isel.size() != 1) return false;
  if (require_single_mod && (int)ev.mod.size() != 1) return false;

  int ibest = isel[0];
  double tmin = ev.ps_t[ibest];
  for (int idx : isel) {
    if (ev.ps_t[idx] < tmin) {
      tmin = ev.ps_t[idx];
      ibest = idx;
    }
  }

  const int imod = FindNearestModuleIndex(ev.mod, tmin);
  if (imod < 0) return false;

  out.ok = true;
  out.t_ps = tmin;
  out.t_mod = ev.mod[imod].t;
  out.dt = out.t_mod - out.t_ps;
  out.ps_amp = ev.ps_amp[ibest];
  out.e_sum = ev.mod[imod].area_sum;
  out.n_ps_sel = (int)isel.size();
  out.n_mod = (int)ev.mod.size();
  return true;
}

static double ScoreFromMetrics(const EvalMetrics& m)
{
  if (m.n_prompt < 200) return -1.0e30;
  const double purity = std::max(0.0, m.purity_sb1);
  const double stat = std::sqrt((double)m.n_prompt);
  const double pod = std::max(0.1, m.prompt_over_delay);
  const double bal = std::max(0.0, m.side_balance);
  return purity * stat * std::pow(pod, 0.25) * (0.5 + 0.5 * bal);
}

static EvalMetrics EvaluateCut(const std::vector<SideEventData>& evA,
                               const std::vector<SideEventData>& evB,
                               double ps_thr_A,
                               double ps_thr_B,
                               const CutDef& cut,
                               int parity_filter)
{
  // parity_filter:
  //   -1: 全イベント
  //    0: even event index
  //    1: odd event index
  EvalMetrics m;
  const size_t n = std::min(evA.size(), evB.size());
  for (size_t i = 0; i < n; ++i) {
    const int p = (int)(i % 2);
    if (parity_filter >= 0 && p != parity_filter) continue;

    SideCandidate cA;
    if (BuildCandidate(evA[i], ps_thr_A, cut.require_single_ps, cut.require_single_mod, cA)) {
      if (cA.e_sum >= cut.e_min) {
        m.n_total++;
        const double adt = std::abs(cA.dt);
        if (adt <= cut.dtw) {
          m.n_prompt++;
          m.n_prompt_A++;
        }
        if (110.0 <= cA.dt && cA.dt <= 170.0) m.n_delay_pos++;
        if (-170.0 <= cA.dt && cA.dt <= -110.0) m.n_delay_neg++;
        if (80.0 <= adt && adt <= 160.0) m.n_sb1++;
        if (200.0 <= adt && adt <= 280.0) m.n_sb2++;
      }
    }

    SideCandidate cB;
    if (BuildCandidate(evB[i], ps_thr_B, cut.require_single_ps, cut.require_single_mod, cB)) {
      if (cB.e_sum >= cut.e_min) {
        m.n_total++;
        const double adt = std::abs(cB.dt);
        if (adt <= cut.dtw) {
          m.n_prompt++;
          m.n_prompt_B++;
        }
        if (110.0 <= cB.dt && cB.dt <= 170.0) m.n_delay_pos++;
        if (-170.0 <= cB.dt && cB.dt <= -110.0) m.n_delay_neg++;
        if (80.0 <= adt && adt <= 160.0) m.n_sb1++;
        if (200.0 <= adt && adt <= 280.0) m.n_sb2++;
      }
    }
  }

  const double w_prompt = cut.dtw;
  const double w_sb = 80.0;
  const double b1 = (double)m.n_sb1 * (w_prompt / w_sb);
  const double b2 = (double)m.n_sb2 * (w_prompt / w_sb);
  if (m.n_prompt > 0) {
    m.purity_sb1 = ((double)m.n_prompt - b1) / (double)m.n_prompt;
    m.purity_sb2 = ((double)m.n_prompt - b2) / (double)m.n_prompt;
  }

  const long n_delay = m.n_delay_pos + m.n_delay_neg;
  m.prompt_over_delay = (double)m.n_prompt / (double)std::max(1L, n_delay);
  const long nmax = std::max(m.n_prompt_A, m.n_prompt_B);
  const long nmin = std::min(m.n_prompt_A, m.n_prompt_B);
  m.side_balance = (nmax > 0) ? (double)nmin / (double)nmax : 0.0;
  m.score = ScoreFromMetrics(m);
  return m;
}

static CutDef OptimizeCutOnSplit(const std::vector<SideEventData>& evA,
                                 const std::vector<SideEventData>& evB,
                                 double ps_auto_A,
                                 double ps_auto_B,
                                 int parity_train,
                                 EvalMetrics& met_train)
{
  // 過学習を抑えるため探索空間を物理的に妥当な範囲へ限定する。
  const std::vector<double> scales_A = {0.35, 0.40, 0.45, 0.50};
  const std::vector<double> scales_B = {0.20, 0.25, 0.30, 0.35};
  const std::vector<double> dtw_set = {16.0, 20.0, 24.0};
  const std::vector<double> e_set = {15000.0, 30000.0, 45000.0};

  CutDef best;
  EvalMetrics bestm;
  bool found = false;

  for (double sA : scales_A) {
    for (double sB : scales_B) {
      for (double w : dtw_set) {
        for (double e : e_set) {
          CutDef c;
          c.scale_A = sA;
          c.scale_B = sB;
          c.dtw = w;
          c.e_min = e;
          c.require_single_ps = true;
          c.require_single_mod = true;
          const double thrA = ps_auto_A * c.scale_A;
          const double thrB = ps_auto_B * c.scale_B;
          EvalMetrics m = EvaluateCut(evA, evB, thrA, thrB, c, parity_train);

          // ロバスト性を強制
          if (m.n_prompt < 700) continue;
          if (m.purity_sb1 < 0.80) continue;
          if (m.prompt_over_delay < 1.80) continue;
          if (m.side_balance < 0.30) continue;

          if (!found || m.score > bestm.score) {
            found = true;
            best = c;
            bestm = m;
          }
        }
      }
    }
  }

  if (!found) {
    // 条件を満たすものがない場合は妥協せず、最もscoreの高い設定を返す。
    for (double sA : scales_A) {
      for (double sB : scales_B) {
        for (double w : dtw_set) {
          for (double e : e_set) {
            CutDef c;
            c.scale_A = sA;
            c.scale_B = sB;
            c.dtw = w;
            c.e_min = e;
            c.require_single_ps = true;
            c.require_single_mod = true;
            const double thrA = ps_auto_A * c.scale_A;
            const double thrB = ps_auto_B * c.scale_B;
            EvalMetrics m = EvaluateCut(evA, evB, thrA, thrB, c, parity_train);
            if (!found || m.score > bestm.score) {
              found = true;
              best = c;
              bestm = m;
            }
          }
        }
      }
    }
  }

  met_train = bestm;
  return best;
}

static void BuildCandidateArrays(const std::vector<SideEventData>& ev,
                                 double ps_thr,
                                 const CutDef& cut,
                                 std::vector<SideCandidate>& cand)
{
  cand.resize(ev.size());
  for (size_t i = 0; i < ev.size(); ++i) {
    SideCandidate c;
    if (BuildCandidate(ev[i], ps_thr, cut.require_single_ps, cut.require_single_mod, c)) {
      if (c.e_sum >= cut.e_min) cand[i] = c;
      else cand[i] = SideCandidate{};
    } else {
      cand[i] = SideCandidate{};
    }
  }
}

static double ComputeMixScale(const TH1D* h_real, const TH1D* h_mix)
{
  // 偶発優位とみなす遠側で正規化
  const double sb_real = h_real->Integral(h_real->GetXaxis()->FindBin(-280.0), h_real->GetXaxis()->FindBin(-200.0))
                       + h_real->Integral(h_real->GetXaxis()->FindBin(-160.0), h_real->GetXaxis()->FindBin(-80.0))
                       + h_real->Integral(h_real->GetXaxis()->FindBin(80.0), h_real->GetXaxis()->FindBin(160.0))
                       + h_real->Integral(h_real->GetXaxis()->FindBin(200.0), h_real->GetXaxis()->FindBin(280.0));
  const double sb_mix = h_mix->Integral(h_mix->GetXaxis()->FindBin(-280.0), h_mix->GetXaxis()->FindBin(-200.0))
                      + h_mix->Integral(h_mix->GetXaxis()->FindBin(-160.0), h_mix->GetXaxis()->FindBin(-80.0))
                      + h_mix->Integral(h_mix->GetXaxis()->FindBin(80.0), h_mix->GetXaxis()->FindBin(160.0))
                      + h_mix->Integral(h_mix->GetXaxis()->FindBin(200.0), h_mix->GetXaxis()->FindBin(280.0));
  if (sb_mix <= 0.0) return 0.0;
  return sb_real / sb_mix;
}

void reassess_michel_selection_run8()
{
  gStyle->SetOptStat(1110);
  gSystem->mkdir("doc/mainexp", true);

  const TString out_pdf = "doc/mainexp/reassess_michel_selection_run8.pdf";

  std::vector<SideEventData> evA;
  std::vector<SideEventData> evB;
  evA.reserve(40000);
  evB.reserve(40000);

  TH1D* hPSAmpA = new TH1D("hPSAmpA_reassess", "PS_A amplitude;counts;Pulses",
                           kPSAmpBins, kPSAmpMin, kPSAmpMax);
  TH1D* hPSAmpB = new TH1D("hPSAmpB_reassess", "PS_B amplitude;counts;Pulses",
                           kPSAmpBins, kPSAmpMin, kPSAmpMax);

  long n_events = 0;
  {
    TString fps[2], fna[4];
    for (int i = 0; i < 2; ++i) fps[i] = Form("%s/%s", kInputDir, kPSFiles[i]);
    for (int i = 0; i < 4; ++i) fna[i] = Form("%s/%s", kInputDir, kNaIFiles[i]);

    std::ifstream in_ps[2];
    std::ifstream in_na[4];
    for (int i = 0; i < 2; ++i) {
      in_ps[i].open(fps[i].Data());
      if (!in_ps[i]) {
        std::cerr << "[reassess] ERROR: cannot open " << fps[i] << std::endl;
        return;
      }
    }
    for (int i = 0; i < 4; ++i) {
      in_na[i].open(fna[i].Data());
      if (!in_na[i]) {
        std::cerr << "[reassess] ERROR: cannot open " << fna[i] << std::endl;
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
      SideEventData eA;
      SideEventData eB;

      for (int ip = 0; ip < 2; ++ip) {
        const double b = ModeBaselineInEvent(buf_ps[ip]);
        std::vector<double> s;
        BuildS(buf_ps[ip], b, s);
        const std::vector<Pulse> ps = FindPulsesSimple(s, kPSSeedThr, kPSEndThr, kPSMinSep, kPSMinWidth);
        for (const auto& q : ps) {
          if (q.amp <= 0.0) continue;
          if (ip == 0) {
            eA.ps_t.push_back((double)q.start);
            eA.ps_amp.push_back(q.amp);
            hPSAmpA->Fill(q.amp);
          } else {
            eB.ps_t.push_back((double)q.start);
            eB.ps_amp.push_back(q.amp);
            hPSAmpB->Fill(q.amp);
          }
        }
      }

      std::vector<FitPulse> fp[4];
      for (int ch = 0; ch < 4; ++ch) {
        const double b = ModeBaselineInEvent(buf_na[ch]);
        std::vector<double> s;
        BuildS(buf_na[ch], b, s);
        fp[ch] = FitPulsesWithDecomposition(s, kTrFix[ch], kDFix[ch]);
      }
      eA.mod = BuildModulePairs(fp[0], fp[1], kPairMatchMaxDt);
      eB.mod = BuildModulePairs(fp[2], fp[3], kPairMatchMaxDt);

      evA.push_back(eA);
      evB.push_back(eB);

      for (int i = 0; i < 2; ++i) buf_ps[i].clear();
      for (int i = 0; i < 4; ++i) buf_na[i].clear();
    }
  }

  int peakA = -1;
  int valleyA = -1;
  int peakB = -1;
  int valleyB = -1;
  const double ps_auto_A = FindPSValleyThreshold(hPSAmpA, peakA, valleyA);
  const double ps_auto_B = FindPSValleyThreshold(hPSAmpB, peakB, valleyB);

  EvalMetrics m_train_even;
  EvalMetrics m_train_odd;
  const CutDef cut_even = OptimizeCutOnSplit(evA, evB, ps_auto_A, ps_auto_B, 0, m_train_even);
  const CutDef cut_odd = OptimizeCutOnSplit(evA, evB, ps_auto_A, ps_auto_B, 1, m_train_odd);

  const double thr_even_A = ps_auto_A * cut_even.scale_A;
  const double thr_even_B = ps_auto_B * cut_even.scale_B;
  const double thr_odd_A = ps_auto_A * cut_odd.scale_A;
  const double thr_odd_B = ps_auto_B * cut_odd.scale_B;

  // split再現性評価
  const EvalMetrics m_even_to_odd = EvaluateCut(evA, evB, thr_even_A, thr_even_B, cut_even, 1);
  const EvalMetrics m_odd_to_even = EvaluateCut(evA, evB, thr_odd_A, thr_odd_B, cut_odd, 0);

  // 最終採用: even学習cutを全体に適用（第三者が再現しやすいよう固定）
  const CutDef cut_final = cut_even;
  const double thr_final_A = thr_even_A;
  const double thr_final_B = thr_even_B;
  const EvalMetrics m_final_all = EvaluateCut(evA, evB, thr_final_A, thr_final_B, cut_final, -1);

  // event-mixing偶発見積
  std::vector<SideCandidate> cA;
  std::vector<SideCandidate> cB;
  BuildCandidateArrays(evA, thr_final_A, cut_final, cA);
  BuildCandidateArrays(evB, thr_final_B, cut_final, cB);

  TH1D* hDtReal = new TH1D("hDtReal_reassess", "dt real (A+B);dt [samples];Candidates", 600, -300, 300);
  TH1D* hDtMix = new TH1D("hDtMix_reassess", "dt mixed (A+B);dt [samples];Candidates", 600, -300, 300);
  TH1D* hEPrompt = new TH1D("hEPrompt_reassess", "E prompt real;E_{sum} [ADC*samples];Counts", 2500, 0, 2.0e5);
  TH1D* hEMixPrompt = new TH1D("hEMixPrompt_reassess", "E prompt mixed;E_{sum} [ADC*samples];Counts", 2500, 0, 2.0e5);

  const std::vector<int> mix_offsets = {137, 503, 997};
  const int n = (int)std::min(cA.size(), cB.size());

  for (int i = 0; i < n; ++i) {
    if (cA[i].ok) {
      hDtReal->Fill(cA[i].dt);
      if (std::abs(cA[i].dt) <= cut_final.dtw) hEPrompt->Fill(cA[i].e_sum);
    }
    if (cB[i].ok) {
      hDtReal->Fill(cB[i].dt);
      if (std::abs(cB[i].dt) <= cut_final.dtw) hEPrompt->Fill(cB[i].e_sum);
    }
  }

  for (int off : mix_offsets) {
    for (int i = 0; i < n; ++i) {
      const int j = (i + off) % n;
      if (cA[i].ok && cA[j].ok) {
        const double dtm = cA[j].t_mod - cA[i].t_ps;
        hDtMix->Fill(dtm);
        if (std::abs(dtm) <= cut_final.dtw) hEMixPrompt->Fill(cA[j].e_sum);
      }
      if (cB[i].ok && cB[j].ok) {
        const double dtm = cB[j].t_mod - cB[i].t_ps;
        hDtMix->Fill(dtm);
        if (std::abs(dtm) <= cut_final.dtw) hEMixPrompt->Fill(cB[j].e_sum);
      }
    }
  }

  const double mix_scale = ComputeMixScale(hDtReal, hDtMix);
  TH1D* hDtMixScaled = (TH1D*)hDtMix->Clone("hDtMixScaled_reassess");
  hDtMixScaled->Scale(mix_scale);

  const double n_real_prompt = hDtReal->Integral(hDtReal->FindBin(-cut_final.dtw), hDtReal->FindBin(+cut_final.dtw));
  const double n_mix_prompt = hDtMixScaled->Integral(hDtMixScaled->FindBin(-cut_final.dtw), hDtMixScaled->FindBin(+cut_final.dtw));
  const double n_real_delay = hDtReal->Integral(hDtReal->FindBin(110.0), hDtReal->FindBin(170.0))
                            + hDtReal->Integral(hDtReal->FindBin(-170.0), hDtReal->FindBin(-110.0));
  const double n_mix_delay = hDtMixScaled->Integral(hDtMixScaled->FindBin(110.0), hDtMixScaled->FindBin(170.0))
                           + hDtMixScaled->Integral(hDtMixScaled->FindBin(-170.0), hDtMixScaled->FindBin(-110.0));
  const double prompt_excess = n_real_prompt - n_mix_prompt;
  const double delay_excess = n_real_delay - n_mix_delay;

  TH1D* hEPromptSub = (TH1D*)hEPrompt->Clone("hEPromptSub_reassess");
  hEPromptSub->SetTitle("E prompt excess (real - scaled mixed);E_{sum} [ADC*samples];Excess");
  hEPromptSub->Add(hEMixPrompt, -mix_scale);

  // 数値を標準出力へ（texに引用するため）
  std::cout << "\n[reassess] === dataset ===\n";
  std::cout << "[reassess] events=" << n_events << "\n";
  std::cout << "[reassess] ps_auto_A=" << ps_auto_A << " ps_auto_B=" << ps_auto_B << "\n";

  auto dump_cut = [&](const char* tag, const CutDef& c, const EvalMetrics& m) {
    std::cout << "[reassess] " << tag
              << " scaleA=" << c.scale_A
              << " scaleB=" << c.scale_B
              << " dtw=" << c.dtw
              << " emin=" << c.e_min
              << " singlePS=" << (c.require_single_ps ? 1 : 0)
              << " singleMod=" << (c.require_single_mod ? 1 : 0)
              << " n_prompt=" << m.n_prompt
              << " n_prompt_A=" << m.n_prompt_A
              << " n_prompt_B=" << m.n_prompt_B
              << " n_delay=" << (m.n_delay_pos + m.n_delay_neg)
              << " puritySB1=" << m.purity_sb1
              << " puritySB2=" << m.purity_sb2
              << " p/d=" << m.prompt_over_delay
              << " balance=" << m.side_balance
              << " score=" << m.score
              << "\n";
  };

  std::cout << "\n[reassess] === split optimization ===\n";
  dump_cut("train-even(best)", cut_even, m_train_even);
  dump_cut("train-even -> test-odd", cut_even, m_even_to_odd);
  dump_cut("train-odd(best)", cut_odd, m_train_odd);
  dump_cut("train-odd -> test-even", cut_odd, m_odd_to_even);
  dump_cut("final(all, using even-cut)", cut_final, m_final_all);

  std::cout << "\n[reassess] === event mixing ===\n";
  std::cout << "[reassess] mix_scale=" << mix_scale << "\n";
  std::cout << "[reassess] real_prompt=" << n_real_prompt
            << " mix_prompt_scaled=" << n_mix_prompt
            << " prompt_excess=" << prompt_excess << "\n";
  std::cout << "[reassess] real_delay=" << n_real_delay
            << " mix_delay_scaled=" << n_mix_delay
            << " delay_excess=" << delay_excess << "\n";
  std::cout << "[reassess] excess_prompt_over_delay="
            << (delay_excess > 0.0 ? prompt_excess / delay_excess : 0.0) << "\n";
  std::cout << "[reassess] E_prompt_mean_real=" << hEPrompt->GetMean()
            << " E_prompt_mean_mix=" << hEMixPrompt->GetMean()
            << " E_prompt_mean_sub=" << hEPromptSub->GetMean()
            << "\n\n";

  // 図作成
  TCanvas* c = new TCanvas("c_reassess", "reassess michel selection", 1400, 1000);
  c->Print(out_pdf + "[");

  c->Clear();
  TPaveText* p = new TPaveText(0.05, 0.05, 0.95, 0.95, "NDC");
  p->SetFillStyle(0);
  p->SetBorderSize(0);
  p->SetTextAlign(12);
  p->SetTextSize(0.027);
  p->AddText("=== Reassessment of Michel-dominant selection (run8) ===");
  p->AddText(Form("Input: %s", kInputDir));
  p->AddText(Form("Events: %ld", n_events));
  p->AddText(Form("PS auto thresholds: A=%.1f, B=%.1f", ps_auto_A, ps_auto_B));
  p->AddText(" ");
  p->AddText("Final cut (chosen from even-train optimization):");
  p->AddText(Form(" - PS_A >= %.2f, PS_B >= %.2f", thr_final_A, thr_final_B));
  p->AddText(Form(" - |dt| <= %.0f samples, E_sum >= %.0f", cut_final.dtw, cut_final.e_min));
  p->AddText(" - single PS per side = ON, single module per side = ON");
  p->AddText(Form(" - prompt=%ld (A=%ld,B=%ld), delay=%ld",
                  m_final_all.n_prompt, m_final_all.n_prompt_A, m_final_all.n_prompt_B,
                  m_final_all.n_delay_pos + m_final_all.n_delay_neg));
  p->AddText(Form(" - purity(SB1)=%.3f, prompt/delay=%.3f, side balance=%.3f",
                  m_final_all.purity_sb1, m_final_all.prompt_over_delay, m_final_all.side_balance));
  p->AddText(" ");
  p->AddText("Event-mixing check:");
  p->AddText(Form(" - prompt excess = %.1f, delay excess = %.1f, ratio = %.2f",
                  prompt_excess, delay_excess,
                  (delay_excess > 0.0 ? prompt_excess / delay_excess : 0.0)));
  p->Draw();
  c->Print(out_pdf);
  delete p;

  c->Clear();
  c->Divide(2, 1);
  c->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  hDtReal->SetLineColor(kBlue + 1);
  hDtReal->SetLineWidth(2);
  hDtMixScaled->SetLineColor(kRed + 1);
  hDtMixScaled->SetLineWidth(2);
  hDtReal->Draw("hist");
  hDtMixScaled->Draw("hist same");
  {
    const double ymax = std::max(1.0, hDtReal->GetMaximum());
    TLine* l1 = new TLine(-cut_final.dtw, 0.0, -cut_final.dtw, ymax);
    TLine* l2 = new TLine(+cut_final.dtw, 0.0, +cut_final.dtw, ymax);
    l1->SetLineStyle(2); l2->SetLineStyle(2);
    l1->SetLineColor(kBlack); l2->SetLineColor(kBlack);
    l1->Draw(); l2->Draw();
    TLegend* lg = new TLegend(0.50, 0.74, 0.92, 0.90);
    lg->SetBorderSize(0);
    lg->SetFillStyle(0);
    lg->AddEntry(hDtReal, "real dt", "l");
    lg->AddEntry(hDtMixScaled, "mixed dt (scaled by far-sideband)", "l");
    lg->Draw();
  }

  c->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  TH1D* hSplit = new TH1D("hSplit_reassess", "Split reproducibility;metric;value", 6, 0.5, 6.5);
  hSplit->GetXaxis()->SetBinLabel(1, "even->odd purity");
  hSplit->GetXaxis()->SetBinLabel(2, "odd->even purity");
  hSplit->GetXaxis()->SetBinLabel(3, "even->odd p/d");
  hSplit->GetXaxis()->SetBinLabel(4, "odd->even p/d");
  hSplit->GetXaxis()->SetBinLabel(5, "even->odd bal");
  hSplit->GetXaxis()->SetBinLabel(6, "odd->even bal");
  hSplit->SetBinContent(1, m_even_to_odd.purity_sb1);
  hSplit->SetBinContent(2, m_odd_to_even.purity_sb1);
  hSplit->SetBinContent(3, m_even_to_odd.prompt_over_delay);
  hSplit->SetBinContent(4, m_odd_to_even.prompt_over_delay);
  hSplit->SetBinContent(5, m_even_to_odd.side_balance);
  hSplit->SetBinContent(6, m_odd_to_even.side_balance);
  hSplit->SetLineColor(kBlue + 1);
  hSplit->SetLineWidth(2);
  hSplit->Draw("hist");
  c->Print(out_pdf);

  c->Clear();
  c->Divide(2, 2);
  c->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogy();
  hEPrompt->SetLineColor(kBlue + 1);
  hEMixPrompt->SetLineColor(kRed + 1);
  hEPrompt->SetLineWidth(2);
  hEMixPrompt->SetLineWidth(2);
  hEPrompt->Draw("hist");
  hEMixPrompt->Draw("hist same");
  {
    TLegend* lg = new TLegend(0.52, 0.74, 0.92, 0.90);
    lg->SetBorderSize(0);
    lg->SetFillStyle(0);
    lg->AddEntry(hEPrompt, "prompt real", "l");
    lg->AddEntry(hEMixPrompt, "prompt mixed", "l");
    lg->Draw();
  }

  c->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  hEPromptSub->SetLineColor(kBlue + 1);
  hEPromptSub->SetLineWidth(2);
  hEPromptSub->Draw("hist");
  {
    TLine* l0 = new TLine(0.0, 0.0, 2.0e5, 0.0);
    l0->SetLineStyle(2);
    l0->SetLineColor(kGray + 2);
    l0->Draw();
  }

  c->cd(3);
  gPad->SetGridx();
  gPad->SetGridy();
  hPSAmpA->SetLineColor(kBlue + 1);
  hPSAmpA->SetLineWidth(2);
  hPSAmpA->Draw("hist");
  {
    TLine* l = new TLine(thr_final_A, 0.0, thr_final_A, std::max(1.0, hPSAmpA->GetMaximum()));
    l->SetLineColor(kRed + 1);
    l->SetLineStyle(2);
    l->Draw();
  }

  c->cd(4);
  gPad->SetGridx();
  gPad->SetGridy();
  hPSAmpB->SetLineColor(kBlue + 1);
  hPSAmpB->SetLineWidth(2);
  hPSAmpB->Draw("hist");
  {
    TLine* l = new TLine(thr_final_B, 0.0, thr_final_B, std::max(1.0, hPSAmpB->GetMaximum()));
    l->SetLineColor(kRed + 1);
    l->SetLineStyle(2);
    l->Draw();
  }
  c->Print(out_pdf);

  c->Print(out_pdf + "]");
  delete c;

  std::cout << "[reassess] output pdf: " << out_pdf << std::endl;
}
