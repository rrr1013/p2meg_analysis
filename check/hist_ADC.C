// check/hist_ADC.C
//
// ADC 値そのものの 1D ヒスト（6ch: NaI4 + PS2）を作って PDF 保存する。
// baseline 推定（最頻値）を見るための step1 用。
// 1ページ目: 6枚のADCヒストのみ
// 2ページ目: メタ情報のみ
//
// 出力: doc/mainexp/hist_ADC_<入力データ名>.pdf
//
// 実行例（リポジトリ直下で）:
//   root -l -q 'check/hist_ADC.C'

#include <iostream>
#include <fstream>

#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
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
static const char* kInputDir = "data/rawdata/mainexp/180°/run1";

// 入力データ名（出力ファイル名に使う）
// 例: "60°_run17"
static const char* kInputTag = "180°_run1";

// 入力ファイル名（必要に応じて手で変更）
static const char* kInputFiles[6] = {
  "wave_NaI_A1_run1.txt",
  "wave_NaI_A2_run1.txt",
  "wave_NaI_B1_run1.txt",
  "wave_NaI_B2_run1.txt",
  "wave_PS_A_run1.txt",
  "wave_PS_B_run1.txt"
};

// 出力ディレクトリ（固定）
static const char* kOutputDir = "doc/mainexp";

// 1D ヒスト設定（必要に応じて調整）
static const int    kNBins1D = 2000;
static const double kHistMin = 12000.0;
static const double kHistMax = 17000.0;

// 描画設定
static const bool   kUseLogY = true;
static const int    kCanvasW = 1200;
static const int    kCanvasH = 900;

// 各パッドに表示するラベル
static const char* kLabels[6] = {"NaI_A1","NaI_A2","NaI_B1","NaI_B2","PS_A","PS_B"};

// ============================================================

struct HistResultADC {
  TH1D*   h = nullptr;
  long    n_samples = 0;
  double  vmin = 0.0;
  double  vmax = 0.0;
  double  mode = 0.0;   // 最頻値（最大ビン中心）
  bool    ok = false;
  TString err;
};

static void SetupPad1D(TPad* p)
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

static void DrawErrorText(const TString& msg)
{
  TLatex t;
  t.SetNDC(true);
  t.SetTextSize(0.05);
  t.DrawLatex(0.14, 0.78, msg.Data());
}

static HistResultADC AnalyzeFileMakeADCHist(const char* filename,
                                            const char* hname,
                                            int nbins,
                                            double h_min,
                                            double h_max)
{
  HistResultADC res;

  std::ifstream file(filename);
  res.h = new TH1D(hname, Form("ADC (%s);ADC counts;Samples", filename),
                   nbins, h_min, h_max);

  if (!file.is_open()) {
    res.err = Form("Could not open: %s", filename);
    res.ok = false;
    return res;
  }

  bool first = true;
  double v;
  while (file >> v) {
    res.h->Fill(v);
    res.n_samples++;

    if (first) {
      res.vmin = v;
      res.vmax = v;
      first = false;
    } else {
      if (v < res.vmin) res.vmin = v;
      if (v > res.vmax) res.vmax = v;
    }
  }
  file.close();

  // 最頻値（最大ビン中心）
  const int maxbin = res.h->GetMaximumBin();
  res.mode = res.h->GetXaxis()->GetBinCenter(maxbin);

  res.ok = true;
  return res;
}

static void BuildInputFilesFromManual(TString out_files[6])
{
  for (int i = 0; i < 6; ++i) {
    out_files[i] = Form("%s/%s", kInputDir, kInputFiles[i]);
  }
}

void hist_ADC()
{
  // 出力ディレクトリ作成（doc/mainexp）
  gSystem->mkdir(kOutputDir, true);

  // 出力 PDF 名（doc/mainexp/hist_ADC_<tag>.pdf）
  TString out_pdf = Form("%s/hist_ADC_%s.pdf", kOutputDir, kInputTag);

  // 入力ファイル（フルパス）
  TString filesS[6];
  BuildInputFilesFromManual(filesS);

  // 解析
  HistResultADC R[6];
  for (int i = 0; i < 6; ++i) {
    R[i] = AnalyzeFileMakeADCHist(filesS[i].Data(),
                                  Form("h_adc_%s", kLabels[i]),
                                  kNBins1D, kHistMin, kHistMax);
    R[i].h->GetXaxis()->SetRangeUser(kHistMin, kHistMax);
  }

  // ============================================================
  // 1ページ目: 6枚のADCヒストのみ（2x3）
  // ============================================================
  gStyle->SetOptStat(1110);

  TCanvas* c1 = new TCanvas("c_hist_adc_page1", "ADC hist page1", kCanvasW, kCanvasH);

  TPad* p11 = new TPad("p11","", 0.00, 0.66, 0.50, 1.00);
  TPad* p12 = new TPad("p12","", 0.50, 0.66, 1.00, 1.00);
  TPad* p21 = new TPad("p21","", 0.00, 0.33, 0.50, 0.66);
  TPad* p22 = new TPad("p22","", 0.50, 0.33, 1.00, 0.66);
  TPad* p31 = new TPad("p31","", 0.00, 0.00, 0.50, 0.33);
  TPad* p32 = new TPad("p32","", 0.50, 0.00, 1.00, 0.33);

  SetupPad1D(p11); SetupPad1D(p12); SetupPad1D(p21);
  SetupPad1D(p22); SetupPad1D(p31); SetupPad1D(p32);

  p11->Draw(); p12->Draw(); p21->Draw(); p22->Draw(); p31->Draw(); p32->Draw();

  TPad* pads[6] = {p11,p12,p21,p22,p31,p32};

  for (int i = 0; i < 6; ++i) {
    pads[i]->cd();
    if (kUseLogY) gPad->SetLogy(1);

    R[i].h->Draw();
    DrawLabelTopLeft(kLabels[i]);

    if (!R[i].ok) {
      DrawErrorText(R[i].err);
    } else {
      TLatex t;
      t.SetNDC(true);
      t.SetTextSize(0.05);
      t.DrawLatex(0.14, 0.80, Form("mode = %.2f", R[i].mode));
    }
  }

  // マルチページPDF開始
  TString out_open = out_pdf; out_open += "(";
  c1->SaveAs(out_open.Data());

  // ============================================================
  // 2ページ目: メタ情報のみ
  // ============================================================
  gStyle->SetOptStat(0);

  TCanvas* c2 = new TCanvas("c_hist_adc_page2", "ADC hist meta", kCanvasW, kCanvasH);
  c2->cd();

  TPaveText* pt = new TPaveText(0.04, 0.06, 0.96, 0.94, "NDC");
  pt->SetFillStyle(0);
  pt->SetBorderSize(0);
  pt->SetTextAlign(12);
  pt->SetTextSize(0.035);

  TDatime now;
  pt->AddText("=== ADC histogram meta ===");
  pt->AddText(Form("InputDir: %s", kInputDir));
  pt->AddText(Form("Tag     : %s", kInputTag));
  pt->AddText(Form("Output  : %s", out_pdf.Data()));
  pt->AddText(Form("Hist    : nbins=%d, range=[%.0f, %.0f], logy(page1)=%d",
                   kNBins1D, kHistMin, kHistMax, (int)kUseLogY));
  pt->AddText(Form("Generated: %04d-%02d-%02d %02d:%02d:%02d",
                   now.GetYear(), now.GetMonth(), now.GetDay(),
                   now.GetHour(), now.GetMinute(), now.GetSecond()));
  pt->AddText(" ");

  pt->AddText("Channel summary (Nsamples / min / max / mode):");
  for (int i = 0; i < 6; ++i) {
    if (R[i].ok) {
      pt->AddText(Form("  %-6s : N=%ld  min=%.1f  max=%.1f  mode=%.2f",
                       kLabels[i], R[i].n_samples, R[i].vmin, R[i].vmax, R[i].mode));
    } else {
      pt->AddText(Form("  %-6s : ERROR (%s)", kLabels[i], R[i].err.Data()));
    }
  }

  pt->Draw();

  // マルチページPDF終了（閉じ括弧）
  TString out_close = out_pdf; out_close += ")";
  c2->SaveAs(out_close.Data());

  std::cout << "Saved: " << out_pdf << std::endl;
}
