// check/plot_ps_threshold_validation.C
//
// 目的:
//   PS振幅ヒストグラム上に、
//   - Michel側ピーク
//   - 谷位置（自動決定しきい値）
// を重ねて表示し、しきい値設定の妥当性を可視化する。
//
// 出力:
//   doc/mainexp/plot_ps_threshold_validation_run8.pdf
//
// 実行:
//   root -l -b -q 'check/plot_ps_threshold_validation.C'

#include <fstream>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"

// PS抽出条件・しきい値探索関数を既存解析と一致させるため、
// studyマクロを読み込む。
#include "study_michel_ps_tagged.C"

void plot_ps_threshold_validation()
{
  gStyle->SetOptStat(1110);
  gSystem->mkdir("doc/mainexp", true);

  TH1D* hPSAmp[2] = {
    new TH1D("hPSAmp_A_thr", "PS_A pulse amplitude;Amplitude [ADC counts];Pulses",
             kPSAmpBins, kPSAmpMin, kPSAmpMax),
    new TH1D("hPSAmp_B_thr", "PS_B pulse amplitude;Amplitude [ADC counts];Pulses",
             kPSAmpBins, kPSAmpMin, kPSAmpMax)
  };

  long n_ps_preselected[2] = {0, 0};

  {
    TString fps[2];
    for (int i = 0; i < 2; ++i) fps[i] = Form("%s/%s", kInputDir, kPSFiles[i]);

    std::ifstream in_ps[2];
    for (int i = 0; i < 2; ++i) {
      in_ps[i].open(fps[i].Data());
      if (!in_ps[i]) {
        std::cerr << "[ps-thr] ERROR: cannot open " << fps[i] << std::endl;
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
        const double bps = ModeBaselineInEvent(buf_ps[ip]);
        std::vector<double> sps;
        BuildS(buf_ps[ip], bps, sps);
        const std::vector<Pulse> ps = FindPulsesSimple(sps, kPSSeedThr, kPSEndThr, kPSMinSep, kPSMinWidth);
        for (const auto& q : ps) {
          if (q.amp <= 0.0) continue;
          hPSAmp[ip]->Fill(q.amp);
          n_ps_preselected[ip]++;
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

  const TString out_pdf = "doc/mainexp/plot_ps_threshold_validation_run8.pdf";
  TCanvas* c = new TCanvas("c_ps_thr", "ps threshold validation", 1400, 900);
  c->Divide(2, 1);

  for (int i = 0; i < 2; ++i) {
    c->cd(i + 1);
    gPad->SetLogy();
    gPad->SetGridx();
    gPad->SetGridy();

    hPSAmp[i]->SetLineColor((i == 0) ? (kBlue + 1) : (kRed + 1));
    hPSAmp[i]->SetLineWidth(2);
    hPSAmp[i]->Draw("hist");

    const double ymax = hPSAmp[i]->GetMaximum();

    TLine* lthr = new TLine(ps_thr_auto[i], 0.8, ps_thr_auto[i], ymax);
    lthr->SetLineColor(kBlack);
    lthr->SetLineWidth(2);
    lthr->SetLineStyle(2);
    lthr->Draw();

    if (ps_peak_bin[i] > 0) {
      const double xpk = hPSAmp[i]->GetBinCenter(ps_peak_bin[i]);
      TLine* lpk = new TLine(xpk, 0.8, xpk, ymax);
      lpk->SetLineColor(kMagenta + 2);
      lpk->SetLineStyle(3);
      lpk->SetLineWidth(2);
      lpk->Draw();
    }

    if (ps_valley_bin[i] > 0) {
      const double xval = hPSAmp[i]->GetBinCenter(ps_valley_bin[i]);
      TLine* lval = new TLine(xval, 0.8, xval, ymax);
      lval->SetLineColor(kGreen + 2);
      lval->SetLineStyle(4);
      lval->SetLineWidth(2);
      lval->Draw();
    }

    const int b_thr = hPSAmp[i]->GetXaxis()->FindBin(ps_thr_auto[i]);
    const double n_all = hPSAmp[i]->Integral(1, hPSAmp[i]->GetNbinsX());
    const double n_above = hPSAmp[i]->Integral(b_thr, hPSAmp[i]->GetNbinsX());
    const double frac = (n_all > 0.0) ? (n_above / n_all) : 0.0;

    std::cout << "[ps-thr] PS_" << (i == 0 ? 'A' : 'B')
              << " threshold=" << ps_thr_auto[i]
              << " all=" << n_all
              << " above=" << n_above
              << " frac=" << frac
              << std::endl;

    TLatex t;
    t.SetNDC(true);
    t.SetTextSize(0.033);
    t.DrawLatex(0.13, 0.88, Form("PS_%c: auto threshold = %.1f", (i == 0 ? 'A' : 'B'), ps_thr_auto[i]));
    t.DrawLatex(0.13, 0.82, Form("total pulses = %ld", n_ps_preselected[i]));
    t.DrawLatex(0.13, 0.76, Form("pulses above threshold = %.0f (%.1f%%)", n_above, 100.0 * frac));

    TLegend* lg = new TLegend(0.55, 0.74, 0.92, 0.90);
    lg->SetBorderSize(0);
    lg->SetFillStyle(0);
    lg->AddEntry(hPSAmp[i], "PS amplitude histogram", "l");
    lg->AddEntry(lthr, "adopted threshold (valley)", "l");
    lg->AddEntry((TObject*)0, "peak marker: magenta", "");
    lg->AddEntry((TObject*)0, "valley marker: green", "");
    lg->Draw();
  }

  c->Print(out_pdf);
  delete c;

  std::cout << "[ps-thr] output pdf: " << out_pdf << std::endl;
}
