// check/hist_PS_peak.C
//
// PS のピーク高さ（baseline - 最小値）の分布ヒストを作る（PS_A / PS_B）。
//
// 重要:
//  - baseline は「各 digitizer 1イベント（固定長）内の ADC 値の最頻値(mode)」で決める。
//  - PS は負パルス想定なので、振幅 A = baseline - v_min (>=0) をピーク高さとする。
//  - 1イベント内に複数パルスがあり得る前提で、全パルスを数える。
//  - パルス判定は seed/end のヒステリシス + min_sep_samples の簡易状態機械。
//  - 確認用に、先頭から5イベント・末尾から5イベントの波形プロットを
//    PDF の最後に追加し、baseline/閾値/検出パルスを可視化する。
//
// 出力: doc/mainexp/hist_PS_peak_<入力データ名>.pdf
//
// 実行例（リポジトリ直下で）:
//   root -l -q 'check/hist_PS_peak.C'

#include <iostream>
#include <fstream>
#include <vector>
#include <deque>
#include <algorithm>
#include <cmath>

#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGraph.h"
#include "TLine.h"
#include "TMarker.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TDatime.h"
#include "TStyle.h"
#include "TString.h"
#include "TSystem.h"

// ============================================================
// 手で変えるパラメータはここだけ
// ============================================================

// 入力ディレクトリ（例: data/rawdata/mainexp/60°/run17）
static const char* kInputDir = "data/rawdata/mainexp/60°/run17";

// 入力データ名（出力ファイル名に使う） 例: "60°_run17"
static const char* kInputTag = "60°_run17";

// 入力ファイル（PS のみ）
static const char* kPSFiles[2] = {
  "wave_PS_A_run17.txt",
  "wave_PS_B_run17.txt"
};

// 出力ディレクトリ（固定）
static const char* kOutputDir = "doc/mainexp";

// 波形→イベント切り出し
static const int    kSamplesPerEvent = 500;  // 1 event あたりサンプル数
static const double kDtNs = 4.0;             // 1bin の時間 [ns]（プロット用）

// パルス検出パラメータ（ピーク高さ分布を見るための暫定値。ヒストを見て手で調整）
static const double kSeedThr = 30.0;   // 開始判定 [ADC counts]（baseline - v）
static const double kEndThr  = 15.0;   // 終了判定 [ADC counts]（baseline - v）
static const int    kMinSepSamples = 100; // パルス終了後のデッドタイム [samples]

// 1D ヒスト（ピーク高さ A=baseline-vmin）設定
static const int    kNBinsPeak = 2000;
static const double kPeakMin  = 0.0;
static const double kPeakMax  = 2000.0;

// 描画設定
static const bool   kUseLogY = true;
static const int    kCanvasW = 1200;
static const int    kCanvasH = 900;

// デバッグ波形プロット（先頭/末尾から何イベント分出すか）
static const int    kNDebugFirst = 5;
static const int    kNDebugLast  = 5;

// ============================================================

struct PulseInfo {
  int start = -1;  // [sample index]
  int peak  = -1;  // [sample index]
  int end   = -1;  // [sample index]
  double amp = 0;  // peak height A = baseline - v_min
};

struct DebugEvent {
  long idx = -1;
  std::vector<double> waveA;
  std::vector<double> waveB;
  double baselineA = 0.0;
  double baselineB = 0.0;
  std::vector<PulseInfo> pulsesA;
  std::vector<PulseInfo> pulsesB;
};

static void SetupPadWave(TPad* p)
{
  p->SetFillStyle(0);
  p->SetLeftMargin(0.10);
  p->SetRightMargin(0.05);
  p->SetTopMargin(0.08);
  p->SetBottomMargin(0.14);
}

static void SetupPadHist(TPad* p)
{
  p->SetFillStyle(0);
  p->SetLeftMargin(0.12);
  p->SetRightMargin(0.05);
  p->SetTopMargin(0.08);
  p->SetBottomMargin(0.12);
}

static void DrawLabelTopLeft(const char* label)
{
  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextSize(0.06);
  lat.DrawLatex(0.14, 0.88, label);
}

// 1イベント内の ADC 値最頻値（mode）を返す。
// - 入力は double だが ADC は整数想定なので四捨五入して int にする。
// - パルスがあっても baseline 付近の出現頻度が最大になりやすい前提。
static double ModeBaselineInEvent(const std::vector<double>& samples)
{
  if (samples.empty()) return 0.0;

  int vmin = (int)std::lround(samples[0]);
  int vmax = vmin;
  for (double x : samples) {
    int v = (int)std::lround(x);
    if (v < vmin) vmin = v;
    if (v > vmax) vmax = v;
  }

  const int range = vmax - vmin + 1;
  if (range <= 0) return (double)vmin;

  // 通常は range ~ O(10^2-10^3) 程度の想定
  std::vector<int> cnt(range, 0);
  for (double x : samples) {
    int v = (int)std::lround(x);
    const int idx = v - vmin;
    if (0 <= idx && idx < range) cnt[idx]++;
  }

  int best_i = 0;
  int best_c = cnt[0];
  for (int i = 1; i < range; ++i) {
    if (cnt[i] > best_c) { best_c = cnt[i]; best_i = i; }
  }

  return (double)(vmin + best_i);
}

// 1イベント内のパルスを検出して pulse list を返す。
// a = baseline - v を振幅とし、seed/end 閾値で状態機械。
// パルス終了後は min_sep_samples 分だけ次の開始を抑制（デッドタイム）。
static std::vector<PulseInfo> FindPulsesInEvent(const std::vector<double>& samples,
                                                double baseline,
                                                double seed_thr,
                                                double end_thr,
                                                int min_sep_samples)
{
  std::vector<PulseInfo> pulses;
  if (samples.empty()) return pulses;

  bool in_pulse = false;
  PulseInfo cur;
  double amax = 0.0;

  int dead = 0;

  for (int i = 0; i < (int)samples.size(); ++i) {

    // デッドタイム中は新規開始だけ抑制（波形は読み進める）
    if (dead > 0) dead--;

    const double a = baseline - samples[i]; // 負パルス想定

    if (!in_pulse) {
      if (dead == 0 && a >= seed_thr) {
        in_pulse = true;
        cur = PulseInfo{};
        cur.start = i;
        cur.peak  = i;
        amax = a;
      }
    } else {
      if (a > amax) {
        amax = a;
        cur.peak = i;
      }
      // 終了判定
      if (a <= end_thr) {
        cur.end = i;
        cur.amp = amax;
        pulses.push_back(cur);

        in_pulse = false;
        amax = 0.0;
        dead = min_sep_samples;
      }
    }
  }

  // イベント末尾までパルスが続いた場合も一旦出す（可視化用にも必要）
  if (in_pulse) {
    cur.end = (int)samples.size() - 1;
    cur.amp = amax;
    pulses.push_back(cur);
  }

  return pulses;
}

// 1ページ（2pad）で波形 + 検出結果を描く
static void DrawDebugWavePage(const DebugEvent& ev, const TString& out_pdf, bool is_last_page)
{
  TCanvas* c = new TCanvas(Form("c_dbg_ev%ld", ev.idx),
                           Form("debug event %ld", ev.idx),
                           kCanvasW, kCanvasH);

  TPad* pA = new TPad(Form("p_dbg_A_%ld", ev.idx), "", 0.00, 0.50, 1.00, 1.00);
  TPad* pB = new TPad(Form("p_dbg_B_%ld", ev.idx), "", 0.00, 0.00, 1.00, 0.50);
  SetupPadWave(pA);
  SetupPadWave(pB);
  pA->Draw();
  pB->Draw();

  auto draw_one = [&](TPad* p,
                      const char* label,
                      const std::vector<double>& wave,
                      double baseline,
                      const std::vector<PulseInfo>& pulses)
  {
    p->cd();

    const int n = (int)wave.size();
    std::vector<double> x(n), y(n);
    for (int i = 0; i < n; ++i) {
      x[i] = i * kDtNs;
      y[i] = wave[i];
    }

    TGraph* gr = new TGraph(n, x.data(), y.data());
    gr->SetTitle(Form("%s  (event %ld);time [ns];ADC", label, ev.idx));
    gr->Draw("AL");

    // baseline / 閾値ライン
    const double y_base = baseline;
    const double y_seed = baseline - kSeedThr;
    const double y_end  = baseline - kEndThr;

    TLine* lbase = new TLine(x.front(), y_base, x.back(), y_base);
    lbase->SetLineStyle(2);
    lbase->SetLineWidth(2);
    lbase->Draw();

    TLine* lseed = new TLine(x.front(), y_seed, x.back(), y_seed);
    lseed->SetLineStyle(3);
    lseed->SetLineWidth(2);
    lseed->Draw();

    TLine* lend  = new TLine(x.front(), y_end,  x.back(), y_end);
    lend->SetLineStyle(3);
    lend->SetLineWidth(2);
    lend->Draw();

    // パルス領域とピーク位置を可視化
    for (size_t k = 0; k < pulses.size(); ++k) {
      const PulseInfo& pu = pulses[k];

      const double xs = pu.start * kDtNs;
      const double xp = pu.peak  * kDtNs;
      const double xe = pu.end   * kDtNs;

      // start/end を縦線
      TLine* ls = new TLine(xs, y_base + 50, xs, y_base - 800);
      ls->SetLineStyle(1);
      ls->SetLineWidth(1);
      ls->Draw();

      TLine* le = new TLine(xe, y_base + 50, xe, y_base - 800);
      le->SetLineStyle(1);
      le->SetLineWidth(1);
      le->Draw();

      // peak をマーカー
      TMarker* mk = new TMarker(xp, wave[pu.peak], 20);
      mk->SetMarkerSize(0.8);
      mk->Draw();
    }

    // テキスト（パラメータ表示）
    TLatex t;
    t.SetNDC(true);
    t.SetTextSize(0.045);
    t.DrawLatex(0.12, 0.86, Form("%s", label));
    t.SetTextSize(0.040);
    t.DrawLatex(0.12, 0.80, Form("baseline(mode in event) = %.2f", baseline));
    t.DrawLatex(0.12, 0.74, Form("seed_thr = %.1f   end_thr = %.1f   min_sep = %d samples",
                                 kSeedThr, kEndThr, kMinSepSamples));
    t.DrawLatex(0.12, 0.68, Form("N(pulses) = %zu", pulses.size()));

    // 凡例的にラインの意味を追記
    t.SetTextSize(0.035);
    t.DrawLatex(0.60, 0.80, "baseline: dashed");
    t.DrawLatex(0.60, 0.74, "seed/end: dotted");
    t.DrawLatex(0.60, 0.68, "peak: marker");
  };

  draw_one(pA, "PS_A", ev.waveA, ev.baselineA, ev.pulsesA);
  draw_one(pB, "PS_B", ev.waveB, ev.baselineB, ev.pulsesB);

  if (is_last_page) {
    TString out_close = out_pdf; out_close += ")";
    c->SaveAs(out_close.Data());
  } else {
    c->SaveAs(out_pdf.Data());
  }
}

// ============================================================

void hist_PS_peak()
{
  // 出力ディレクトリ作成（doc/mainexp）
  gSystem->mkdir(kOutputDir, true);

  // 出力 PDF
  TString out_pdf = Form("%s/hist_PS_peak_%s.pdf", kOutputDir, kInputTag);

  // 入力ファイル（フルパス）
  TString fA = Form("%s/%s", kInputDir, kPSFiles[0]);
  TString fB = Form("%s/%s", kInputDir, kPSFiles[1]);

  std::ifstream inA(fA.Data());
  std::ifstream inB(fB.Data());
  if (!inA.is_open() || !inB.is_open()) {
    std::cerr << "ERROR: could not open input files.\n";
    std::cerr << "  PS_A: " << fA << "\n";
    std::cerr << "  PS_B: " << fB << "\n";
    return;
  }

  // ヒスト（ピーク高さ）
  TH1D* hPeakA = new TH1D("hPeakA", "PS_A peak height;A = baseline - v_{min} [ADC];Pulses",
                          kNBinsPeak, kPeakMin, kPeakMax);
  TH1D* hPeakB = new TH1D("hPeakB", "PS_B peak height;A = baseline - v_{min} [ADC];Pulses",
                          kNBinsPeak, kPeakMin, kPeakMax);

  long n_events = 0;
  long n_pulses_A = 0;
  long n_pulses_B = 0;

  std::vector<double> bufA; bufA.reserve(kSamplesPerEvent);
  std::vector<double> bufB; bufB.reserve(kSamplesPerEvent);

  std::vector<DebugEvent> first_dbg;
  first_dbg.reserve(kNDebugFirst);

  std::deque<DebugEvent> last_dbg; // 常に最後の kNDebugLast を保持

  double vA = 0.0, vB = 0.0;
  while (true) {
    if (!(inA >> vA)) break;
    if (!(inB >> vB)) break;

    bufA.push_back(vA);
    bufB.push_back(vB);

    if ((int)bufA.size() == kSamplesPerEvent && (int)bufB.size() == kSamplesPerEvent) {

      // baseline を mode で決める（イベントごと）
      const double baseA = ModeBaselineInEvent(bufA);
      const double baseB = ModeBaselineInEvent(bufB);

      // パルス抽出（全パルス）
      std::vector<PulseInfo> pulsesA = FindPulsesInEvent(bufA, baseA, kSeedThr, kEndThr, kMinSepSamples);
      std::vector<PulseInfo> pulsesB = FindPulsesInEvent(bufB, baseB, kSeedThr, kEndThr, kMinSepSamples);

      // ヒスト fill
      for (const auto& p : pulsesA) { hPeakA->Fill(p.amp); n_pulses_A++; }
      for (const auto& p : pulsesB) { hPeakB->Fill(p.amp); n_pulses_B++; }

      // デバッグ保存（先頭/末尾）
      DebugEvent ev;
      ev.idx = n_events;
      ev.waveA = bufA;
      ev.waveB = bufB;
      ev.baselineA = baseA;
      ev.baselineB = baseB;
      ev.pulsesA = pulsesA;
      ev.pulsesB = pulsesB;

      if ((int)first_dbg.size() < kNDebugFirst) {
        first_dbg.push_back(ev);
      }
      last_dbg.push_back(ev);
      while ((int)last_dbg.size() > kNDebugLast) last_dbg.pop_front();

      // 次イベントへ
      n_events++;
      bufA.clear();
      bufB.clear();
    }
  }

  inA.close();
  inB.close();

  // ============================================================
  // PDF 1ページ目: PS_A / PS_B ピーク高さヒスト（2枚）
  // ============================================================
  gStyle->SetOptStat(1110);

  TCanvas* c1 = new TCanvas("c_peak_page1", "PS peak hist", kCanvasW, kCanvasH);

  TPad* pTop = new TPad("pTop","", 0.00, 0.50, 1.00, 1.00);
  TPad* pBot = new TPad("pBot","", 0.00, 0.00, 1.00, 0.50);
  SetupPadHist(pTop);
  SetupPadHist(pBot);
  pTop->Draw();
  pBot->Draw();

  pTop->cd();
  if (kUseLogY) gPad->SetLogy(1);
  hPeakA->Draw();
  DrawLabelTopLeft("PS_A");

  pBot->cd();
  if (kUseLogY) gPad->SetLogy(1);
  hPeakB->Draw();
  DrawLabelTopLeft("PS_B");

  TString out_open = out_pdf; out_open += "(";
  c1->SaveAs(out_open.Data());

  // ============================================================
  // PDF 2ページ目: メタ情報
  // ============================================================
  gStyle->SetOptStat(0);

  TCanvas* c2 = new TCanvas("c_peak_meta", "meta", kCanvasW, kCanvasH);
  c2->cd();

  TPaveText* pt = new TPaveText(0.04, 0.06, 0.96, 0.94, "NDC");
  pt->SetFillStyle(0);
  pt->SetBorderSize(0);
  pt->SetTextAlign(12);
  pt->SetTextSize(0.035);

  TDatime now;
  pt->AddText("=== PS peak height histogram meta ===");
  pt->AddText(Form("InputDir : %s", kInputDir));
  pt->AddText(Form("PS_A file: %s", fA.Data()));
  pt->AddText(Form("PS_B file: %s", fB.Data()));
  pt->AddText(Form("Output   : %s", out_pdf.Data()));
  pt->AddText(" ");
  pt->AddText("Definitions:");
  pt->AddText("  baseline = mode(ADC) within each digitizer event");
  pt->AddText("  peak height A = baseline - v_min (negative pulse assumed)");
  pt->AddText(" ");
  pt->AddText(Form("Event segmentation: samples_per_event = %d, dt = %.1f ns", kSamplesPerEvent, kDtNs));
  pt->AddText(Form("Pulse finding: seed_thr = %.1f, end_thr = %.1f, min_sep = %d samples",
                   kSeedThr, kEndThr, kMinSepSamples));
  pt->AddText(Form("Peak hist: nbins = %d, range = [%.1f, %.1f], logy(page1)=%d",
                   kNBinsPeak, kPeakMin, kPeakMax, (int)kUseLogY));
  pt->AddText(" ");
  pt->AddText(Form("Processed events = %ld", n_events));
  pt->AddText(Form("Total pulses: PS_A = %ld, PS_B = %ld", n_pulses_A, n_pulses_B));
  pt->AddText(Form("Generated: %04d-%02d-%02d %02d:%02d:%02d",
                   now.GetYear(), now.GetMonth(), now.GetDay(),
                   now.GetHour(), now.GetMinute(), now.GetSecond()));
  pt->Draw();

  c2->SaveAs(out_pdf.Data());

  // ============================================================
  // デバッグ波形ページ（先頭5 + 末尾5）
  // ============================================================

  // debug event list を構成
  std::vector<DebugEvent> dbg;
  dbg.reserve(kNDebugFirst + kNDebugLast);

  for (const auto& e : first_dbg) dbg.push_back(e);
  for (const auto& e : last_dbg)  dbg.push_back(e);

  // 末尾ページ判定して PDF を閉じる
  if (dbg.empty()) {
    TString out_close = out_pdf; out_close += ")";
    c2->SaveAs(out_close.Data());
    std::cout << "Saved: " << out_pdf << "\n";
    return;
  }

  for (size_t i = 0; i < dbg.size(); ++i) {
    const bool is_last = (i == dbg.size() - 1);
    DrawDebugWavePage(dbg[i], out_pdf, is_last);
  }

  std::cout << "Saved: " << out_pdf << "\n";
}
