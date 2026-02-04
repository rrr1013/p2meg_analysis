// check/hist_ADCint.C
//
// run番号だけ指定して、6ファイル（NaI 4 + PS 2）のADC積分ヒストを作って
// 1枚のPDFに配置して保存する。
// 出力: hist_ADCint_runXX.pdf
//
// 実行例（リポジトリ直下で）:
//   root -l -q 'check/hist_ADCint.C(12)'

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TDatime.h"
#include "TStyle.h"
#include "TString.h"

// 入力ファイルが置かれているディレクトリ
// 例: データがリポジトリ直下にあるなら "."、data/ にあるなら "data"。
// （run番号以外を触る必要があるとすればここだけです）
static const char* kInputDir = "data/mockdata_wave";

// ヒスト設定
static const int    kNBins   = 500;
static const double kHistMin = 3000.0;
static const double kHistMax = 6e5;

// 波形→イベント切り出し
static const int kSamplesPerEvent = 500;
static const int kBaselineSamples = 50;

// 各パッドに表示するラベル
static const char* kLabels[6] = {"NaIA1","NaIA2","NaIB1","NaIB2","PS1","PS2"};

// ============================================================

struct HistResult {
  TH1D* h = nullptr;
  long  n_events = 0;
  double vmin = 0.0;
  double vmax = 0.0;
  bool ok = false;
  TString err;
};

static HistResult MakeIntegralHist(const char* filename,
                                   const char* hname,
                                   int nbins,
                                   double h_min,
                                   double h_max,
                                   int samples_per_event,
                                   int baseline_samples)
{
  HistResult res;

  std::ifstream file(filename);
  if (!file.is_open()) {
    res.err = Form("Could not open: %s", filename);
    res.h = new TH1D(hname, Form("Spectrum (%s);Integral;Counts", filename),
                     nbins, h_min, h_max);
    res.ok = false;
    return res;
  }

  std::vector<double> event_buffer;
  event_buffer.reserve(samples_per_event);

  std::vector<double> integrals;
  integrals.reserve(10000);

  double val;
  while (file >> val) {
    event_buffer.push_back(val);

    if ((int)event_buffer.size() == samples_per_event) {
      // baseline（先頭 baseline_samples の平均）
      double baseline = 0.0;
      for (int i = 0; i < baseline_samples; ++i) baseline += event_buffer[i];
      baseline /= (double)baseline_samples;

      // integral（負パルス想定: baseline - v を全サンプルで総和）
      double integral = 0.0;
      for (double v : event_buffer) integral += (baseline - v);

      integrals.push_back(integral);
      event_buffer.clear();
    }
  }
  file.close();

  res.h = new TH1D(hname, Form("Spectrum (%s);Integral;Counts", filename),
                   nbins, h_min, h_max);

  for (double s : integrals) res.h->Fill(s);

  res.n_events = (long)integrals.size();
  if (!integrals.empty()) {
    auto it_min = std::min_element(integrals.begin(), integrals.end());
    auto it_max = std::max_element(integrals.begin(), integrals.end());
    res.vmin = *it_min;
    res.vmax = *it_max;
  }

  res.ok = true;
  return res;
}

static void SetupPad(TPad* p)
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

// run番号 → 入力ファイルパスを構成
static void BuildInputFiles(int run, TString out_files[6])
{
  // 画像で示された命名規則に合わせる
  // wave_NaI_A1_run12.txt 等
  out_files[0] = Form("%s/wave_NaI_A1_run%d.txt", kInputDir, run);
  out_files[1] = Form("%s/wave_NaI_A2_run%d.txt", kInputDir, run);
  out_files[2] = Form("%s/wave_NaI_B1_run%d.txt", kInputDir, run);
  out_files[3] = Form("%s/wave_NaI_B2_run%d.txt", kInputDir, run);
  out_files[4] = Form("%s/wave_PS_A_run%d.txt",   kInputDir, run);
  out_files[5] = Form("%s/wave_PS_B_run%d.txt",   kInputDir, run);
}

void hist_ADCint(int run = 12)
{
  gStyle->SetOptStat(1110);

  TString filesS[6];
  BuildInputFiles(run, filesS);

  HistResult R[6];
  for (int i = 0; i < 6; ++i) {
    R[i] = MakeIntegralHist(filesS[i].Data(),
                            Form("h_%s_run%d", kLabels[i], run),
                            kNBins, kHistMin, kHistMax,
                            kSamplesPerEvent, kBaselineSamples);
    R[i].h->GetXaxis()->SetRangeUser(kHistMin, kHistMax);
  }

  // 出力PDF名（末尾に _run12 形式）
  TString out_pdf = Form("doc/hist_ADCint_run%d.pdf", run);

  // レイアウト:
  // y=0.50-1.00: NaI 2x2
  // y=0.25-0.50: PS 1x2
  // y=0.00-0.25: メタ情報
  TCanvas* c = new TCanvas("c_hist", "ADC Integral Histograms", 1200, 900);

  TPad* pNaI11 = new TPad("pNaI11","", 0.00, 0.75, 0.50, 1.00);
  TPad* pNaI12 = new TPad("pNaI12","", 0.50, 0.75, 1.00, 1.00);
  TPad* pNaI21 = new TPad("pNaI21","", 0.00, 0.50, 0.50, 0.75);
  TPad* pNaI22 = new TPad("pNaI22","", 0.50, 0.50, 1.00, 0.75);

  TPad* pPS1   = new TPad("pPS1",  "", 0.00, 0.25, 0.50, 0.50);
  TPad* pPS2   = new TPad("pPS2",  "", 0.50, 0.25, 1.00, 0.50);

  TPad* pMeta  = new TPad("pMeta", "", 0.00, 0.00, 1.00, 0.25);

  SetupPad(pNaI11); SetupPad(pNaI12); SetupPad(pNaI21); SetupPad(pNaI22);
  SetupPad(pPS1);   SetupPad(pPS2);

  pMeta->SetLeftMargin(0.04);
  pMeta->SetRightMargin(0.02);
  pMeta->SetTopMargin(0.10);
  pMeta->SetBottomMargin(0.10);

  pNaI11->Draw(); pNaI12->Draw(); pNaI21->Draw(); pNaI22->Draw();
  pPS1->Draw();   pPS2->Draw();
  pMeta->Draw();

  // NaI 4枚
  pNaI11->cd(); R[0].h->Draw(); DrawLabelTopLeft(kLabels[0]);
  if (!R[0].ok) { TLatex t; t.SetNDC(true); t.SetTextSize(0.05); t.DrawLatex(0.14,0.78,R[0].err); }

  pNaI12->cd(); R[1].h->Draw(); DrawLabelTopLeft(kLabels[1]);
  if (!R[1].ok) { TLatex t; t.SetNDC(true); t.SetTextSize(0.05); t.DrawLatex(0.14,0.78,R[1].err); }

  pNaI21->cd(); R[2].h->Draw(); DrawLabelTopLeft(kLabels[2]);
  if (!R[2].ok) { TLatex t; t.SetNDC(true); t.SetTextSize(0.05); t.DrawLatex(0.14,0.78,R[2].err); }

  pNaI22->cd(); R[3].h->Draw(); DrawLabelTopLeft(kLabels[3]);
  if (!R[3].ok) { TLatex t; t.SetNDC(true); t.SetTextSize(0.05); t.DrawLatex(0.14,0.78,R[3].err); }

  // PS 2枚
  pPS1->cd(); R[4].h->Draw(); DrawLabelTopLeft(kLabels[4]);
  if (!R[4].ok) { TLatex t; t.SetNDC(true); t.SetTextSize(0.05); t.DrawLatex(0.14,0.78,R[4].err); }

  pPS2->cd(); R[5].h->Draw(); DrawLabelTopLeft(kLabels[5]);
  if (!R[5].ok) { TLatex t; t.SetNDC(true); t.SetTextSize(0.05); t.DrawLatex(0.14,0.78,R[5].err); }

  // メタ情報
  pMeta->cd();
  TPaveText* pt = new TPaveText(0.02, 0.05, 0.98, 0.95, "NDC");
  pt->SetFillStyle(0);
  pt->SetBorderSize(0);
  pt->SetTextAlign(12);
  pt->SetTextSize(0.10);

  TDatime now;
  pt->AddText(Form("Run: %d    Output: %s    Generated: %04d-%02d-%02d %02d:%02d:%02d",
                   run, out_pdf.Data(),
                   now.GetYear(), now.GetMonth(), now.GetDay(),
                   now.GetHour(), now.GetMinute(), now.GetSecond()));

  pt->AddText(Form("InputDir: %s    samples_per_event=%d, baseline_samples=%d, nbins=%d, hist_range=[%.0f, %.0f]",
                   kInputDir, kSamplesPerEvent, kBaselineSamples, kNBins, kHistMin, kHistMax));

  pt->AddText("Input files / event counts / min-max:");
  for (int i = 0; i < 6; ++i) {
    if (R[i].ok) {
      pt->AddText(Form("  %-6s : %s  | N=%ld  min=%.1f  max=%.1f",
                       kLabels[i], filesS[i].Data(), R[i].n_events, R[i].vmin, R[i].vmax));
    } else {
      pt->AddText(Form("  %-6s : %s  | ERROR (%s)",
                       kLabels[i], filesS[i].Data(), R[i].err.Data()));
    }
  }
  pt->Draw();

  c->SaveAs(out_pdf.Data());
  std::cout << "Saved: " << out_pdf << std::endl;
}
