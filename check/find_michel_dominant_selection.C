// check/find_michel_dominant_selection.C
//
// 目的:
//   run8データから「Michel崩壊由来が支配的」と言える選別条件を
//   データ駆動で探索する。
//
// 方針:
//   - PS_A と NaI_A, PS_B と NaI_B を同側対応で扱う
//   - 変える条件:
//       1) PS閾値スケール（自動閾値に対する倍率）
//       2) prompt窓 |dt|<=w
//       3) NaI和積分しきい値 E_sum >= E_min
//       4) side内で単一PS要求 / 単一module要求
//   - 指標:
//       * prompt収量 N_prompt
//       * 遅延帯収量 N_delay
//       * sideband背景見積からの純度 proxy
//       * PS振幅とNaI和積分の相関（prompt / delay）
//
// 出力:
//   doc/mainexp/find_michel_dominant_selection_run8.pdf
//   doc/mainexp/find_michel_dominant_selection_run8.txt
//
// 実行:
//   root -l -b -q 'check/find_michel_dominant_selection.C'

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
#include "TH2D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TSystem.h"

// 既存再構成ロジックを再利用
#include "study_michel_ps_tagged.C"

struct SideEventData {
  std::vector<double> ps_t;    // PS候補時刻（seed条件を満たす全候補）
  std::vector<double> ps_amp;  // 同じ順序で振幅
  std::vector<ModulePulse> mod; // NaI module対応パルス
};

struct Candidate {
  double tps = 0.0;
  double ps_amp = 0.0;
  double dt = 0.0;
  double e_sum = 0.0;
  int n_ps_sel = 0;
  int n_mod = 0;
  bool ok = false;
};

struct EvalResult {
  double scale_A = 1.0;
  double scale_B = 1.0;
  double dtw = 20.0;
  double e_min = 0.0;
  bool require_single_ps = false;
  bool require_single_mod = false;

  long n_candidates = 0;
  long n_prompt = 0;
  long n_prompt_A = 0;
  long n_prompt_B = 0;
  long n_delay_pos = 0;
  long n_delay_neg = 0;
  long n_sb1 = 0; // |dt| in [80,160]
  long n_sb2 = 0; // |dt| in [200,280]

  double b_est_sb1 = 0.0;
  double b_est_sb2 = 0.0;
  double s_est_sb1 = 0.0;
  double s_est_sb2 = 0.0;
  double purity_sb1 = 0.0;
  double purity_sb2 = 0.0;
  double prompt_over_delay = 0.0;
  double side_balance = 0.0; // min(A,B)/max(A,B)
  double corr_prompt = 0.0;  // corr(PS amp, E_sum) in prompt
  double corr_delay = 0.0;   // corr(PS amp, E_sum) in delay
  double score = -1.0e30;
};

static double PearsonCorr(const std::vector<double>& x, const std::vector<double>& y)
{
  if (x.size() != y.size()) return 0.0;
  if (x.size() < 3) return 0.0;

  const size_t n = x.size();
  double mx = 0.0;
  double my = 0.0;
  for (size_t i = 0; i < n; ++i) {
    mx += x[i];
    my += y[i];
  }
  mx /= (double)n;
  my /= (double)n;

  double sxx = 0.0;
  double syy = 0.0;
  double sxy = 0.0;
  for (size_t i = 0; i < n; ++i) {
    const double dx = x[i] - mx;
    const double dy = y[i] - my;
    sxx += dx * dx;
    syy += dy * dy;
    sxy += dx * dy;
  }
  if (sxx <= 0.0 || syy <= 0.0) return 0.0;
  return sxy / std::sqrt(sxx * syy);
}

static bool BuildSideCandidate(const SideEventData& ev,
                               double ps_thr,
                               bool require_single_ps,
                               bool require_single_mod,
                               Candidate& out)
{
  out = Candidate{};

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

  int ibest_ps = isel[0];
  double tmin = ev.ps_t[ibest_ps];
  for (int idx : isel) {
    if (ev.ps_t[idx] < tmin) {
      tmin = ev.ps_t[idx];
      ibest_ps = idx;
    }
  }
  const double tps = ev.ps_t[ibest_ps];

  const int imod = FindNearestModuleIndex(ev.mod, tps);
  if (imod < 0) return false;

  out.tps = tps;
  out.ps_amp = ev.ps_amp[ibest_ps];
  out.dt = ev.mod[imod].t - tps;
  out.e_sum = ev.mod[imod].area_sum;
  out.n_ps_sel = (int)isel.size();
  out.n_mod = (int)ev.mod.size();
  out.ok = true;
  return true;
}

static double ScoreResult(const EvalResult& r)
{
  // 物理的に望ましい条件:
  // 1) prompt純度が高い
  // 2) prompt統計が十分ある
  // 3) 遅延ピークが小さい
  // 4) A/Bの偏りが極端でない
  // 5) promptでPS-E相関が見える（delayより強い）
  if (r.n_prompt < 200) return -1.0e30;
  const double purity = std::max(0.0, r.purity_sb1);
  const double stat = std::sqrt((double)r.n_prompt);
  const double pd = std::max(0.1, r.prompt_over_delay);
  const double bal = std::max(0.0, r.side_balance);
  const double corr_term = std::max(0.0, r.corr_prompt - std::max(0.0, r.corr_delay));
  return purity * stat * std::pow(pd, 0.25) * (0.6 + 0.4 * bal) * (1.0 + 0.6 * corr_term);
}

static EvalResult EvaluateCut(const std::vector<SideEventData>& evA,
                              const std::vector<SideEventData>& evB,
                              double ps_thr_A,
                              double ps_thr_B,
                              double scale_A,
                              double scale_B,
                              double dtw,
                              double e_min,
                              bool require_single_ps,
                              bool require_single_mod)
{
  EvalResult r;
  r.scale_A = scale_A;
  r.scale_B = scale_B;
  r.dtw = dtw;
  r.e_min = e_min;
  r.require_single_ps = require_single_ps;
  r.require_single_mod = require_single_mod;

  std::vector<double> x_prompt;
  std::vector<double> y_prompt;
  std::vector<double> x_delay;
  std::vector<double> y_delay;
  x_prompt.reserve(50000);
  y_prompt.reserve(50000);
  x_delay.reserve(30000);
  y_delay.reserve(30000);

  const size_t n = std::min(evA.size(), evB.size());
  for (size_t i = 0; i < n; ++i) {
    Candidate cA;
    if (BuildSideCandidate(evA[i], ps_thr_A, require_single_ps, require_single_mod, cA)) {
      if (cA.e_sum >= e_min) {
        r.n_candidates++;
        const double adt = std::abs(cA.dt);
        if (adt <= dtw) {
          r.n_prompt++;
          r.n_prompt_A++;
          x_prompt.push_back(cA.ps_amp);
          y_prompt.push_back(cA.e_sum);
        }
        if (110.0 <= cA.dt && cA.dt <= 170.0) {
          r.n_delay_pos++;
          x_delay.push_back(cA.ps_amp);
          y_delay.push_back(cA.e_sum);
        }
        if (-170.0 <= cA.dt && cA.dt <= -110.0) {
          r.n_delay_neg++;
          x_delay.push_back(cA.ps_amp);
          y_delay.push_back(cA.e_sum);
        }
        if (80.0 <= adt && adt <= 160.0) r.n_sb1++;
        if (200.0 <= adt && adt <= 280.0) r.n_sb2++;
      }
    }

    Candidate cB;
    if (BuildSideCandidate(evB[i], ps_thr_B, require_single_ps, require_single_mod, cB)) {
      if (cB.e_sum >= e_min) {
        r.n_candidates++;
        const double adt = std::abs(cB.dt);
        if (adt <= dtw) {
          r.n_prompt++;
          r.n_prompt_B++;
          x_prompt.push_back(cB.ps_amp);
          y_prompt.push_back(cB.e_sum);
        }
        if (110.0 <= cB.dt && cB.dt <= 170.0) {
          r.n_delay_pos++;
          x_delay.push_back(cB.ps_amp);
          y_delay.push_back(cB.e_sum);
        }
        if (-170.0 <= cB.dt && cB.dt <= -110.0) {
          r.n_delay_neg++;
          x_delay.push_back(cB.ps_amp);
          y_delay.push_back(cB.e_sum);
        }
        if (80.0 <= adt && adt <= 160.0) r.n_sb1++;
        if (200.0 <= adt && adt <= 280.0) r.n_sb2++;
      }
    }
  }

  const double w_prompt = dtw;
  const double w_sb = 80.0; // [80,160] は幅80
  r.b_est_sb1 = (double)r.n_sb1 * (w_prompt / w_sb);
  r.b_est_sb2 = (double)r.n_sb2 * (w_prompt / w_sb);
  r.s_est_sb1 = (double)r.n_prompt - r.b_est_sb1;
  r.s_est_sb2 = (double)r.n_prompt - r.b_est_sb2;

  if (r.n_prompt > 0) {
    r.purity_sb1 = r.s_est_sb1 / (double)r.n_prompt;
    r.purity_sb2 = r.s_est_sb2 / (double)r.n_prompt;
  }

  const long n_delay = r.n_delay_pos + r.n_delay_neg;
  r.prompt_over_delay = (double)r.n_prompt / (double)std::max(1L, n_delay);

  const long nmax_side = std::max(r.n_prompt_A, r.n_prompt_B);
  const long nmin_side = std::min(r.n_prompt_A, r.n_prompt_B);
  r.side_balance = (nmax_side > 0) ? (double)nmin_side / (double)nmax_side : 0.0;

  r.corr_prompt = PearsonCorr(x_prompt, y_prompt);
  r.corr_delay = PearsonCorr(x_delay, y_delay);
  r.score = ScoreResult(r);
  return r;
}

static void FillBestCutHists(const std::vector<SideEventData>& evA,
                             const std::vector<SideEventData>& evB,
                             double ps_thr_A,
                             double ps_thr_B,
                             const EvalResult& best,
                             TH1D* hDtAll,
                             TH1D* hDtA,
                             TH1D* hDtB,
                             TH1D* hEPrompt,
                             TH1D* hESb1,
                             TH1D* hESb2,
                             TH2D* h2Prompt,
                             TH2D* h2Delay)
{
  const size_t n = std::min(evA.size(), evB.size());
  for (size_t i = 0; i < n; ++i) {
    Candidate cA;
    if (BuildSideCandidate(evA[i], ps_thr_A, best.require_single_ps, best.require_single_mod, cA) &&
        cA.e_sum >= best.e_min) {
      const double adt = std::abs(cA.dt);
      hDtAll->Fill(cA.dt);
      hDtA->Fill(cA.dt);
      if (adt <= best.dtw) {
        hEPrompt->Fill(cA.e_sum);
        h2Prompt->Fill(cA.ps_amp, cA.e_sum);
      }
      if (80.0 <= adt && adt <= 160.0) hESb1->Fill(cA.e_sum);
      if (200.0 <= adt && adt <= 280.0) hESb2->Fill(cA.e_sum);
      if ((110.0 <= cA.dt && cA.dt <= 170.0) || (-170.0 <= cA.dt && cA.dt <= -110.0)) {
        h2Delay->Fill(cA.ps_amp, cA.e_sum);
      }
    }

    Candidate cB;
    if (BuildSideCandidate(evB[i], ps_thr_B, best.require_single_ps, best.require_single_mod, cB) &&
        cB.e_sum >= best.e_min) {
      const double adt = std::abs(cB.dt);
      hDtAll->Fill(cB.dt);
      hDtB->Fill(cB.dt);
      if (adt <= best.dtw) {
        hEPrompt->Fill(cB.e_sum);
        h2Prompt->Fill(cB.ps_amp, cB.e_sum);
      }
      if (80.0 <= adt && adt <= 160.0) hESb1->Fill(cB.e_sum);
      if (200.0 <= adt && adt <= 280.0) hESb2->Fill(cB.e_sum);
      if ((110.0 <= cB.dt && cB.dt <= 170.0) || (-170.0 <= cB.dt && cB.dt <= -110.0)) {
        h2Delay->Fill(cB.ps_amp, cB.e_sum);
      }
    }
  }
}

void find_michel_dominant_selection()
{
  gStyle->SetOptStat(1110);
  gSystem->mkdir("doc/mainexp", true);

  const TString out_pdf = "doc/mainexp/find_michel_dominant_selection_run8.pdf";
  const TString out_txt = "doc/mainexp/find_michel_dominant_selection_run8.txt";

  std::vector<SideEventData> evA;
  std::vector<SideEventData> evB;
  evA.reserve(40000);
  evB.reserve(40000);

  TH1D* hPSAmpA = new TH1D("hPSAmpA_scan", "PS_A amp;counts;Pulses", kPSAmpBins, kPSAmpMin, kPSAmpMax);
  TH1D* hPSAmpB = new TH1D("hPSAmpB_scan", "PS_B amp;counts;Pulses", kPSAmpBins, kPSAmpMin, kPSAmpMax);

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
        std::cerr << "[michel-find] ERROR: cannot open " << fps[i] << std::endl;
        return;
      }
    }
    for (int i = 0; i < 4; ++i) {
      in_na[i].open(fna[i].Data());
      if (!in_na[i]) {
        std::cerr << "[michel-find] ERROR: cannot open " << fna[i] << std::endl;
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

      // PS候補（seed条件を満たす全候補）を sideごとに保持
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

      // NaI分解
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
  const double ps_thr_auto_A = FindPSValleyThreshold(hPSAmpA, peakA, valleyA);
  const double ps_thr_auto_B = FindPSValleyThreshold(hPSAmpB, peakB, valleyB);

  std::cout << "[michel-find] events loaded: " << n_events << std::endl;
  std::cout << "[michel-find] auto PS thresholds A/B = "
            << ps_thr_auto_A << " / " << ps_thr_auto_B << std::endl;

  const std::vector<double> scales_A = {0.4, 0.5, 0.6, 0.7, 0.85, 1.0, 1.15, 1.3};
  const std::vector<double> scales_B = {0.25, 0.35, 0.45, 0.55, 0.7, 0.85, 1.0, 1.2};
  const std::vector<double> dtws = {8.0, 10.0, 12.0, 16.0, 20.0, 24.0};
  const std::vector<double> e_mins = {0.0, 15000.0, 30000.0, 45000.0, 60000.0, 80000.0, 100000.0};
  const std::vector<bool> single_ps_opts = {false, true};
  const std::vector<bool> single_mod_opts = {false, true};

  std::vector<EvalResult> results;
  results.reserve(scales_A.size() * scales_B.size() * dtws.size() * e_mins.size() * single_ps_opts.size() * single_mod_opts.size());

  for (double scA : scales_A) {
    const double thrA = ps_thr_auto_A * scA;
    for (double scB : scales_B) {
      const double thrB = ps_thr_auto_B * scB;
      for (double w : dtws) {
        for (double emin : e_mins) {
          for (bool sp : single_ps_opts) {
            for (bool sm : single_mod_opts) {
              EvalResult r = EvaluateCut(evA, evB, thrA, thrB, scA, scB, w, emin, sp, sm);
              results.push_back(r);
            }
          }
        }
      }
    }
  }

  std::sort(results.begin(), results.end(),
            [](const EvalResult& a, const EvalResult& b) {
              if (a.score != b.score) return a.score > b.score;
              if (a.purity_sb1 != b.purity_sb1) return a.purity_sb1 > b.purity_sb1;
              return a.n_prompt > b.n_prompt;
            });

  EvalResult best;
  EvalResult best_default;
  bool found_best = false;
  bool found_default = false;

  for (const auto& r : results) {
    if (std::abs(r.scale_A - 1.0) < 1e-9 &&
        std::abs(r.scale_B - 1.0) < 1e-9 &&
        std::abs(r.dtw - 20.0) < 1e-9 &&
        std::abs(r.e_min - 0.0) < 1e-9 &&
        !r.require_single_ps &&
        !r.require_single_mod) {
      best_default = r;
      found_default = true;
      break;
    }
  }

  // まず「Michel支配」を主張しやすいロバスト条件を優先して採用する。
  // 条件:
  //  - sideごと単一PS/単一module
  //  - A/Bの極端な非対称を避ける
  //  - prompt優位が明確
  for (const auto& r : results) {
    if (!r.require_single_ps) continue;
    if (!r.require_single_mod) continue;
    if (r.n_prompt < 700) continue;
    if (r.purity_sb1 < 0.80) continue;
    if (r.prompt_over_delay < 2.00) continue;
    if (r.side_balance < 0.30) continue;
    best = r;
    found_best = true;
    break;
  }

  // ロバスト条件が満たせない場合は、従来の高スコアを採用。
  if (!found_best) {
    for (const auto& r : results) {
      if (r.n_prompt < 400) continue;
      if (r.purity_sb1 < 0.50) continue;
      best = r;
      found_best = true;
      break;
    }
  }
  if (!found_best && !results.empty()) {
    best = results.front();
    found_best = true;
  }

  if (!found_best) {
    std::cerr << "[michel-find] ERROR: no valid result" << std::endl;
    return;
  }

  const double best_thr_A = ps_thr_auto_A * best.scale_A;
  const double best_thr_B = ps_thr_auto_B * best.scale_B;

  TH1D* hDtAll = new TH1D("hDtAll_best", "dt all sides (best cut);dt [samples];Candidates", 600, -300, 300);
  TH1D* hDtA = new TH1D("hDtA_best", "dt A side (best cut);dt [samples];Candidates", 600, -300, 300);
  TH1D* hDtB = new TH1D("hDtB_best", "dt B side (best cut);dt [samples];Candidates", 600, -300, 300);
  TH1D* hEPrompt = new TH1D("hEPrompt_best", "E prompt (|dt|<=w);E_{sum} [ADC*samples];Counts", 2500, 0, 2.0e5);
  TH1D* hESb1 = new TH1D("hESb1_best", "E sideband1 (80<=|dt|<=160);E_{sum} [ADC*samples];Counts", 2500, 0, 2.0e5);
  TH1D* hESb2 = new TH1D("hESb2_best", "E sideband2 (200<=|dt|<=280);E_{sum} [ADC*samples];Counts", 2500, 0, 2.0e5);
  TH2D* h2Prompt = new TH2D("h2Prompt_best", "Prompt: PS amp vs E_{sum};PS amp [counts];E_{sum} [ADC*samples]",
                            160, 0, 800, 160, 0, 2.0e5);
  TH2D* h2Delay = new TH2D("h2Delay_best", "Delay: PS amp vs E_{sum};PS amp [counts];E_{sum} [ADC*samples]",
                           160, 0, 800, 160, 0, 2.0e5);

  FillBestCutHists(evA, evB, best_thr_A, best_thr_B, best,
                   hDtAll, hDtA, hDtB, hEPrompt, hESb1, hESb2, h2Prompt, h2Delay);

  TH1D* hEsub1 = (TH1D*)hEPrompt->Clone("hEsub1_best");
  TH1D* hEsub2 = (TH1D*)hEPrompt->Clone("hEsub2_best");
  hEsub1->SetTitle("Prompt - scaled SB1;E_{sum} [ADC*samples];Counts (arb.)");
  hEsub2->SetTitle("Prompt - scaled SB2;E_{sum} [ADC*samples];Counts (arb.)");
  const double sf1 = best.dtw / 80.0;
  const double sf2 = best.dtw / 80.0;
  hEsub1->Add(hESb1, -sf1);
  hEsub2->Add(hESb2, -sf2);

  // 出力テキスト
  {
    std::ofstream ofs(out_txt.Data());
    ofs << "# Michel-dominant selection search (run8)\n";
    ofs << "# input_dir = " << kInputDir << "\n";
    ofs << "# events = " << n_events << "\n";
    ofs << "# PS auto thresholds A/B = " << ps_thr_auto_A << " / " << ps_thr_auto_B << "\n";
    ofs << "#\n";
    ofs << "# Default reference cut (scaleA=1, scaleB=1, dtw=20, e_min=0, singlePS=0, singleMod=0)\n";
    if (found_default) {
      ofs << "default.n_prompt = " << best_default.n_prompt << "\n";
      ofs << "default.n_prompt_A = " << best_default.n_prompt_A << "\n";
      ofs << "default.n_prompt_B = " << best_default.n_prompt_B << "\n";
      ofs << "default.purity_sb1 = " << best_default.purity_sb1 << "\n";
      ofs << "default.purity_sb2 = " << best_default.purity_sb2 << "\n";
      ofs << "default.prompt_over_delay = " << best_default.prompt_over_delay << "\n";
      ofs << "default.side_balance = " << best_default.side_balance << "\n";
    } else {
      ofs << "default.not_found = 1\n";
    }
    ofs << "\n";

    ofs << "# Best selection\n";
    ofs << "best.scale_A = " << best.scale_A << "\n";
    ofs << "best.scale_B = " << best.scale_B << "\n";
    ofs << "best.ps_thr_A = " << best_thr_A << "\n";
    ofs << "best.ps_thr_B = " << best_thr_B << "\n";
    ofs << "best.dtw = " << best.dtw << "\n";
    ofs << "best.e_min = " << best.e_min << "\n";
    ofs << "best.require_single_ps = " << (best.require_single_ps ? 1 : 0) << "\n";
    ofs << "best.require_single_mod = " << (best.require_single_mod ? 1 : 0) << "\n";
    ofs << "best.n_candidates = " << best.n_candidates << "\n";
    ofs << "best.n_prompt = " << best.n_prompt << "\n";
    ofs << "best.n_prompt_A = " << best.n_prompt_A << "\n";
    ofs << "best.n_prompt_B = " << best.n_prompt_B << "\n";
    ofs << "best.n_delay_pos = " << best.n_delay_pos << "\n";
    ofs << "best.n_delay_neg = " << best.n_delay_neg << "\n";
    ofs << "best.n_sb1 = " << best.n_sb1 << "\n";
    ofs << "best.n_sb2 = " << best.n_sb2 << "\n";
    ofs << "best.b_est_sb1 = " << best.b_est_sb1 << "\n";
    ofs << "best.b_est_sb2 = " << best.b_est_sb2 << "\n";
    ofs << "best.purity_sb1 = " << best.purity_sb1 << "\n";
    ofs << "best.purity_sb2 = " << best.purity_sb2 << "\n";
    ofs << "best.prompt_over_delay = " << best.prompt_over_delay << "\n";
    ofs << "best.side_balance = " << best.side_balance << "\n";
    ofs << "best.corr_prompt = " << best.corr_prompt << "\n";
    ofs << "best.corr_delay = " << best.corr_delay << "\n";
    ofs << "best.score = " << best.score << "\n";
    ofs << "\n";

    ofs << "# Top 20 by score\n";
    ofs << "# rank scale_A scale_B dtw e_min single_ps single_mod n_prompt purity_sb1 purity_sb2 "
        << "prompt_over_delay side_balance corr_prompt corr_delay score\n";
    for (size_t i = 0; i < std::min<size_t>(20, results.size()); ++i) {
      const auto& r = results[i];
      ofs << (i + 1) << " "
          << r.scale_A << " "
          << r.scale_B << " "
          << r.dtw << " "
          << r.e_min << " "
          << (r.require_single_ps ? 1 : 0) << " "
          << (r.require_single_mod ? 1 : 0) << " "
          << r.n_prompt << " "
          << r.purity_sb1 << " "
          << r.purity_sb2 << " "
          << r.prompt_over_delay << " "
          << r.side_balance << " "
          << r.corr_prompt << " "
          << r.corr_delay << " "
          << r.score << "\n";
    }
    ofs << "\n";
    ofs << "# All results (sorted by score)\n";
    ofs << "# rank scale_A scale_B dtw e_min single_ps single_mod "
        << "n_candidates n_prompt n_prompt_A n_prompt_B n_delay_pos n_delay_neg "
        << "n_sb1 n_sb2 purity_sb1 purity_sb2 prompt_over_delay side_balance "
        << "corr_prompt corr_delay score\n";
    for (size_t i = 0; i < results.size(); ++i) {
      const auto& r = results[i];
      ofs << (i + 1) << " "
          << r.scale_A << " "
          << r.scale_B << " "
          << r.dtw << " "
          << r.e_min << " "
          << (r.require_single_ps ? 1 : 0) << " "
          << (r.require_single_mod ? 1 : 0) << " "
          << r.n_candidates << " "
          << r.n_prompt << " "
          << r.n_prompt_A << " "
          << r.n_prompt_B << " "
          << r.n_delay_pos << " "
          << r.n_delay_neg << " "
          << r.n_sb1 << " "
          << r.n_sb2 << " "
          << r.purity_sb1 << " "
          << r.purity_sb2 << " "
          << r.prompt_over_delay << " "
          << r.side_balance << " "
          << r.corr_prompt << " "
          << r.corr_delay << " "
          << r.score << "\n";
    }
  }

  // 図出力
  TCanvas* c = new TCanvas("c_find_michel", "find michel dominant", 1400, 1000);
  c->Print(out_pdf + "[");

  c->Clear();
  TPaveText* p = new TPaveText(0.05, 0.05, 0.95, 0.95, "NDC");
  p->SetFillStyle(0);
  p->SetBorderSize(0);
  p->SetTextAlign(12);
  p->SetTextSize(0.028);
  p->AddText("=== Michel-dominant selection search (run8, side-matched PS-NaI) ===");
  p->AddText(Form("Input: %s", kInputDir));
  p->AddText(Form("Events loaded: %ld", n_events));
  p->AddText(Form("PS auto thresholds: A=%.1f, B=%.1f", ps_thr_auto_A, ps_thr_auto_B));
  p->AddText(" ");
  p->AddText("Best cut (robust criteria):");
  p->AddText(Form(" - scale A/B = %.2f / %.2f -> PSthr A/B = %.1f / %.1f",
                  best.scale_A, best.scale_B, best_thr_A, best_thr_B));
  p->AddText(Form(" - |dt| <= %.0f, E_sum >= %.0f", best.dtw, best.e_min));
  p->AddText(Form(" - single PS = %s, single module = %s",
                  (best.require_single_ps ? "ON" : "OFF"),
                  (best.require_single_mod ? "ON" : "OFF")));
  p->AddText(Form(" - N_prompt = %ld (A=%ld, B=%ld), N_delay(pos+neg)=%ld",
                  best.n_prompt, best.n_prompt_A, best.n_prompt_B, best.n_delay_pos + best.n_delay_neg));
  p->AddText(Form(" - prompt/delay = %.2f", best.prompt_over_delay));
  p->AddText(Form(" - purity proxy (SB1[80,160]) = %.3f", best.purity_sb1));
  p->AddText(Form(" - purity proxy (SB2[200,280]) = %.3f", best.purity_sb2));
  p->AddText(Form(" - side balance min/max = %.3f", best.side_balance));
  p->AddText(Form(" - corr(PSamp,E) prompt/delay = %.3f / %.3f", best.corr_prompt, best.corr_delay));
  if (found_default) {
    p->AddText(" ");
    p->AddText("Default reference (scaleA=1, scaleB=1, dtw=20, e_min=0, no single constraints):");
    p->AddText(Form(" - N_prompt=%ld (A=%ld,B=%ld), prompt/delay=%.2f",
                    best_default.n_prompt, best_default.n_prompt_A, best_default.n_prompt_B,
                    best_default.prompt_over_delay));
    p->AddText(Form(" - purity SB1/SB2 = %.3f / %.3f, side balance=%.3f",
                    best_default.purity_sb1, best_default.purity_sb2, best_default.side_balance));
  }
  p->Draw();
  c->Print(out_pdf);
  delete p;

  c->Clear();
  c->Divide(2, 2);
  c->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  hDtAll->SetLineColor(kBlue + 1);
  hDtAll->SetLineWidth(2);
  hDtAll->Draw("hist");
  {
    const double ymax = std::max(1.0, hDtAll->GetMaximum());
    TLine* l0 = new TLine(0.0, 0.0, 0.0, ymax);
    TLine* lp1 = new TLine(-best.dtw, 0.0, -best.dtw, ymax);
    TLine* lp2 = new TLine(+best.dtw, 0.0, +best.dtw, ymax);
    TLine* ld1 = new TLine(110.0, 0.0, 110.0, ymax);
    TLine* ld2 = new TLine(170.0, 0.0, 170.0, ymax);
    l0->SetLineStyle(3);
    lp1->SetLineStyle(2); lp2->SetLineStyle(2);
    ld1->SetLineStyle(2); ld2->SetLineStyle(2);
    lp1->SetLineColor(kRed + 1); lp2->SetLineColor(kRed + 1);
    ld1->SetLineColor(kMagenta + 2); ld2->SetLineColor(kMagenta + 2);
    l0->Draw(); lp1->Draw(); lp2->Draw(); ld1->Draw(); ld2->Draw();
  }

  c->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  hDtA->SetLineColor(kBlue + 1);
  hDtA->SetLineWidth(2);
  hDtA->Draw("hist");

  c->cd(3);
  gPad->SetGridx();
  gPad->SetGridy();
  hDtB->SetLineColor(kRed + 1);
  hDtB->SetLineWidth(2);
  hDtB->Draw("hist");

  c->cd(4);
  gPad->SetGridx();
  gPad->SetGridy();
  TH1D* hTopPrompt = new TH1D("hTopPrompt", "Top20: N_{prompt};rank;N_{prompt}", 20, 0.5, 20.5);
  TH1D* hTopPurity = new TH1D("hTopPurity", "Top20: purity(SB1);rank;purity", 20, 0.5, 20.5);
  for (int i = 0; i < 20 && i < (int)results.size(); ++i) {
    hTopPrompt->SetBinContent(i + 1, results[i].n_prompt);
    hTopPurity->SetBinContent(i + 1, results[i].purity_sb1);
  }
  hTopPrompt->SetLineColor(kBlue + 1);
  hTopPrompt->SetLineWidth(2);
  hTopPrompt->Draw("hist");
  hTopPurity->SetLineColor(kRed + 1);
  hTopPurity->SetLineWidth(2);
  hTopPurity->Draw("hist same");
  {
    TLegend* lg = new TLegend(0.52, 0.74, 0.92, 0.90);
    lg->SetBorderSize(0);
    lg->SetFillStyle(0);
    lg->AddEntry(hTopPrompt, "N_prompt", "l");
    lg->AddEntry(hTopPurity, "purity(SB1)", "l");
    lg->Draw();
  }
  c->Print(out_pdf);

  c->Clear();
  c->Divide(2, 2);
  c->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogy();
  hEPrompt->SetLineColor(kBlue + 1);
  hESb1->SetLineColor(kRed + 1);
  hESb2->SetLineColor(kGreen + 2);
  hEPrompt->SetLineWidth(2);
  hESb1->SetLineWidth(2);
  hESb2->SetLineWidth(2);
  hEPrompt->Draw("hist");
  hESb1->Draw("hist same");
  hESb2->Draw("hist same");
  {
    TLegend* lg = new TLegend(0.50, 0.72, 0.92, 0.90);
    lg->SetBorderSize(0);
    lg->SetFillStyle(0);
    lg->AddEntry(hEPrompt, Form("prompt |dt|<=%.0f", best.dtw), "l");
    lg->AddEntry(hESb1, "SB1: 80<=|dt|<=160", "l");
    lg->AddEntry(hESb2, "SB2: 200<=|dt|<=280", "l");
    lg->Draw();
  }

  c->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  hEsub1->SetLineColor(kRed + 1);
  hEsub2->SetLineColor(kBlue + 1);
  hEsub1->SetLineWidth(2);
  hEsub2->SetLineWidth(2);
  hEsub1->Draw("hist");
  hEsub2->Draw("hist same");
  {
    TLine* l0 = new TLine(0.0, 0.0, 2.0e5, 0.0);
    l0->SetLineStyle(2);
    l0->SetLineColor(kGray + 2);
    l0->Draw();
    TLegend* lg = new TLegend(0.50, 0.74, 0.92, 0.90);
    lg->SetBorderSize(0);
    lg->SetFillStyle(0);
    lg->AddEntry(hEsub1, Form("prompt - (%.2f)*SB1", sf1), "l");
    lg->AddEntry(hEsub2, Form("prompt - (%.2f)*SB2", sf2), "l");
    lg->Draw();
  }

  c->cd(3);
  gPad->SetGridx();
  gPad->SetGridy();
  h2Prompt->Draw("colz");

  c->cd(4);
  gPad->SetGridx();
  gPad->SetGridy();
  h2Delay->Draw("colz");
  c->Print(out_pdf);

  c->Print(out_pdf + "]");
  delete c;

  // 参考root保存
  {
    TFile froot("doc/mainexp/find_michel_dominant_selection_run8.root", "RECREATE");
    hPSAmpA->Write();
    hPSAmpB->Write();
    hDtAll->Write();
    hDtA->Write();
    hDtB->Write();
    hEPrompt->Write();
    hESb1->Write();
    hESb2->Write();
    hEsub1->Write();
    hEsub2->Write();
    h2Prompt->Write();
    h2Delay->Write();
    froot.Close();
  }

  std::cout << "[michel-find] output pdf: " << out_pdf << std::endl;
  std::cout << "[michel-find] output txt: " << out_txt << std::endl;
}
