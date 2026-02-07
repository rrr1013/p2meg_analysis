// check/hist_ADCint.C
//
// 手書きで指定した6ファイル（NaI 4 + PS 2）のADC積分ヒストを作って
// 1枚目に1Dヒスト6枚（NaI 4 + PS 2 + メタ）、2ページ目以降に全組合せ(6C2=15)の2Dヒストを
// 6枚/ページで描画してPDF保存する。
// 出力: kOutputPdfName で指定したPDF
//
// 実行例（リポジトリ直下で）:
//   root -l -q 'check/hist_ADCint.C'

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <utility>

#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TDatime.h"
#include "TStyle.h"
#include "TString.h"

// 入力ファイルが置かれているディレクトリ
static const char* kInputDir = "data/rawdata";

// 入力ファイル名（必要に応じてここを手書きで変更する）
static const char* kInputFiles[6] = {
  "wave_NaI_A1_run12.txt",
  "wave_NaI_A2_run12.txt",
  "wave_NaI_B1_run12.txt",
  "wave_NaI_B2_run12.txt",
  "wave_PS_A_run12.txt",
  "wave_PS_B_run12.txt"
};

// static const char* kInputFiles[6] = {
//   "wave_NaI_A1_run12.txt",
//   "wave_NaI_A2_run12.txt",
//   "wave_NaI_B1_run12.txt",
//   "wave_NaI_B2_run12.txt",
//   "wave_PS_A_run12.txt",
//   "wave_PS_B_run12.txt"
// };

// 出力PDF名（必要に応じて変更する）
static const char* kOutputPdfName = "doc/hist_ADCint_manual.pdf";



// 1Dヒスト設定
static const int    kNBins1D = 500;
static const double kHistMin = 3000.0;
static const double kHistMax = 6e5;

// 2Dヒスト設定（重い場合はここを下げてください）
static const int kNBins2D = 200;

// 波形→イベント切り出し
static const int kSamplesPerEvent = 500;
static const int kBaselineSamples = 50;

// 各パッドに表示するラベル
static const char* kLabels[6] = {"NaIA1","NaIA2","NaIB1","NaIB2","PS1","PS2"};

// ============================================================

struct HistResult {
  TH1D* h = nullptr;
  std::vector<double> integrals;
  long  n_events = 0;
  double vmin = 0.0;
  double vmax = 0.0;
  bool ok = false;
  TString err;
};

static void DrawLabelTopLeft(const char* label)
{
  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextSize(0.06);
  lat.DrawLatex(0.14, 0.88, label);
}

static void DrawErrorText(const TString& msg)
{
  TLatex t;
  t.SetNDC(true);
  t.SetTextSize(0.05);
  t.DrawLatex(0.14, 0.78, msg.Data());
}

static void SetupPad1D(TPad* p)
{
  p->SetFillStyle(0);
  p->SetLeftMargin(0.12);
  p->SetRightMargin(0.05);
  p->SetTopMargin(0.08);
  p->SetBottomMargin(0.12);
}

static void SetupPad2D(TPad* p)
{
  p->SetFillStyle(0);
  p->SetLeftMargin(0.12);
  p->SetRightMargin(0.13); // COLZ のカラーバー分
  p->SetTopMargin(0.08);
  p->SetBottomMargin(0.12);
}

static HistResult AnalyzeFileMake1D(const char* filename,
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

  res.integrals.swap(integrals);
  res.ok = true;
  return res;
}

// 入力ファイル名配列をフルパスに変換
static void BuildInputFilesFromManual(TString out_files[6])
{
  for (int i = 0; i < 6; ++i) {
    out_files[i] = Form("%s/%s", kInputDir, kInputFiles[i]);
  }
}

void hist_ADCint()
{
  // ---------- 6本を解析して 1D と integrals を作る ----------
  TString filesS[6];
  BuildInputFilesFromManual(filesS);

  HistResult R[6];
  for (int i = 0; i < 6; ++i) {
    R[i] = AnalyzeFileMake1D(filesS[i].Data(),
                             Form("h_%s", kLabels[i]),
                             kNBins1D, kHistMin, kHistMax,
                             kSamplesPerEvent, kBaselineSamples);
    R[i].h->GetXaxis()->SetRangeUser(kHistMin, kHistMax);
  }

  // 出力PDF名（doc/ に保存）
  TString out_pdf = kOutputPdfName;

  // ============================================================
  // 1ページ目（1D: NaI4 + PS2 + メタ）
  // ============================================================
  gStyle->SetOptStat(1110);

  TCanvas* c1 = new TCanvas("c_hist_1d", "ADC Integral 1D", 1200, 900);

  TPad* pNaI11 = new TPad("pNaI11","", 0.00, 0.75, 0.50, 1.00);
  TPad* pNaI12 = new TPad("pNaI12","", 0.50, 0.75, 1.00, 1.00);
  TPad* pNaI21 = new TPad("pNaI21","", 0.00, 0.50, 0.50, 0.75);
  TPad* pNaI22 = new TPad("pNaI22","", 0.50, 0.50, 1.00, 0.75);

  TPad* pPS1   = new TPad("pPS1",  "", 0.00, 0.25, 0.50, 0.50);
  TPad* pPS2   = new TPad("pPS2",  "", 0.50, 0.25, 1.00, 0.50);

  TPad* pMeta  = new TPad("pMeta", "", 0.00, 0.00, 1.00, 0.25);

  SetupPad1D(pNaI11); SetupPad1D(pNaI12); SetupPad1D(pNaI21); SetupPad1D(pNaI22);
  SetupPad1D(pPS1);   SetupPad1D(pPS2);

  pMeta->SetLeftMargin(0.04);
  pMeta->SetRightMargin(0.02);
  pMeta->SetTopMargin(0.10);
  pMeta->SetBottomMargin(0.10);

  pNaI11->Draw(); pNaI12->Draw(); pNaI21->Draw(); pNaI22->Draw();
  pPS1->Draw();   pPS2->Draw();
  pMeta->Draw();

  pNaI11->cd(); R[0].h->Draw(); DrawLabelTopLeft(kLabels[0]); if (!R[0].ok) DrawErrorText(R[0].err);
  pNaI12->cd(); R[1].h->Draw(); DrawLabelTopLeft(kLabels[1]); if (!R[1].ok) DrawErrorText(R[1].err);
  pNaI21->cd(); R[2].h->Draw(); DrawLabelTopLeft(kLabels[2]); if (!R[2].ok) DrawErrorText(R[2].err);
  pNaI22->cd(); R[3].h->Draw(); DrawLabelTopLeft(kLabels[3]); if (!R[3].ok) DrawErrorText(R[3].err);
  pPS1->cd();   R[4].h->Draw(); DrawLabelTopLeft(kLabels[4]); if (!R[4].ok) DrawErrorText(R[4].err);
  pPS2->cd();   R[5].h->Draw(); DrawLabelTopLeft(kLabels[5]); if (!R[5].ok) DrawErrorText(R[5].err);

  // メタ情報
  pMeta->cd();
  TPaveText* pt = new TPaveText(0.02, 0.05, 0.98, 0.95, "NDC");
  pt->SetFillStyle(0);
  pt->SetBorderSize(0);
  pt->SetTextAlign(12);
  pt->SetTextSize(0.10);

  TDatime now;
  pt->AddText(Form("Output: %s    Generated: %04d-%02d-%02d %02d:%02d:%02d",
                   out_pdf.Data(),
                   now.GetYear(), now.GetMonth(), now.GetDay(),
                   now.GetHour(), now.GetMinute(), now.GetSecond()));
  pt->AddText(Form("samples_per_event=%d, baseline_samples=%d, nbins1D=%d, nbins2D=%d, range=[%.0f, %.0f]",
                   kSamplesPerEvent, kBaselineSamples, kNBins1D, kNBins2D, kHistMin, kHistMax));

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

  // マルチページPDF開始
  TString out_pdf_open = out_pdf; out_pdf_open += "(";
  c1->SaveAs(out_pdf_open.Data());

  // ============================================================
  // 2ページ目以降（2D: 6C2 = 15通り、6枚/ページ）
  // ============================================================

  // 組合せリスト
  std::vector<std::pair<int,int>> pairs;
  pairs.reserve(15);
  for (int i = 0; i < 6; ++i) {
    for (int j = i + 1; j < 6; ++j) {
      pairs.emplace_back(i, j);
    }
  }

  // 2Dヒストを作成（必要ならここでまとめて作る）
  std::vector<TH2D*> h2s;
  h2s.reserve(pairs.size());

  for (size_t p = 0; p < pairs.size(); ++p) {
    int i = pairs[p].first;
    int j = pairs[p].second;

    TString name  = Form("h2_%s_vs_%s", kLabels[i], kLabels[j]);
    TString title = Form("%s vs %s;%s integral;%s integral",
                         kLabels[i], kLabels[j], kLabels[i], kLabels[j]);

    TH2D* h2 = new TH2D(name.Data(), title.Data(),
                        kNBins2D, kHistMin, kHistMax,
                        kNBins2D, kHistMin, kHistMax);

    // どちらかが読めていない場合は空のまま（描画時にエラー表示）
    if (R[i].ok && R[j].ok) {
      const long nfill = (long)std::min(R[i].integrals.size(), R[j].integrals.size());
      for (long k = 0; k < nfill; ++k) {
        h2->Fill(R[i].integrals[k], R[j].integrals[k]);
      }
    }

    h2s.push_back(h2);
  }

  gStyle->SetOptStat(0); // 2Dはstat boxを消す（見づらいので）

  const int per_page = 6;
  const int n_pages_2d = (int)((h2s.size() + per_page - 1) / per_page);

  for (int page = 0; page < n_pages_2d; ++page) {

    TCanvas* c2 = new TCanvas(Form("c_hist_2d_p%d", page+1),
                              Form("ADC Integral 2D page %d", page+1),
                              1200, 900);

    // 2 x 3 の6パッド（全面使用）
    TPad* pads[6];
    pads[0] = new TPad(Form("p2d_%d_0",page), "", 0.00, 0.66, 0.50, 1.00);
    pads[1] = new TPad(Form("p2d_%d_1",page), "", 0.50, 0.66, 1.00, 1.00);
    pads[2] = new TPad(Form("p2d_%d_2",page), "", 0.00, 0.33, 0.50, 0.66);
    pads[3] = new TPad(Form("p2d_%d_3",page), "", 0.50, 0.33, 1.00, 0.66);
    pads[4] = new TPad(Form("p2d_%d_4",page), "", 0.00, 0.00, 0.50, 0.33);
    pads[5] = new TPad(Form("p2d_%d_5",page), "", 0.50, 0.00, 1.00, 0.33);

    for (int k = 0; k < 6; ++k) {
      SetupPad2D(pads[k]);
      pads[k]->Draw();
    }

    for (int slot = 0; slot < per_page; ++slot) {
      const int idx = page * per_page + slot;
      pads[slot]->cd();

      if (idx >= (int)h2s.size()) {
        // 最終ページの余りスロットは空欄
        TLatex t;
        t.SetNDC(true);
        t.SetTextSize(0.07);
        t.DrawLatex(0.20, 0.50, " ");
        continue;
      }

      int i = pairs[idx].first;
      int j = pairs[idx].second;

      if (!(R[i].ok && R[j].ok)) {
        // どちらかの入力が読めていない
        TLatex t;
        t.SetNDC(true);
        t.SetTextSize(0.06);
        t.DrawLatex(0.14, 0.80, Form("%s vs %s", kLabels[i], kLabels[j]));
        t.SetTextSize(0.05);
        if (!R[i].ok) t.DrawLatex(0.14, 0.68, R[i].err.Data());
        if (!R[j].ok) t.DrawLatex(0.14, 0.58, R[j].err.Data());
        continue;
      }

      h2s[idx]->Draw("COLZ");

      // パッド内ラベル（ペア名）
      TLatex lat;
      lat.SetNDC(true);
      lat.SetTextSize(0.06);
      lat.DrawLatex(0.14, 0.88, Form("%s vs %s", kLabels[i], kLabels[j]));

      // 使ったイベント数（ズレがある場合の確認用）
      const long nfill = (long)std::min(R[i].integrals.size(), R[j].integrals.size());
      lat.SetTextSize(0.05);
      lat.DrawLatex(0.14, 0.80, Form("N(fill) = %ld", nfill));
    }

    // PDFへ追加（最後のページだけ閉じ括弧）
    const bool is_last_page = (page == n_pages_2d - 1);
    if (is_last_page) {
      TString out_pdf_close = out_pdf; out_pdf_close += ")";
      c2->SaveAs(out_pdf_close.Data());
    } else {
      c2->SaveAs(out_pdf.Data());
    }
  }

  std::cout << "Saved: " << out_pdf << std::endl;
}
