// macros/plot_michel_e_hist.C
//
// Michel 偏極測定用 Ee-only データのヒストグラム確認
//
// 入力: Ee-only .dat（1行1事象、1列Ee[MeV]。空行/#コメント行は無視）
// 出力: doc/michel_e_hist_<basename>.pdf（2ページ）
//   1) meta（入出力・設定・統計）
//   2) Ee ヒストグラム
//
// 実行例（リポジトリ直下から）:
//   root -l -q 'macros/plot_michel_e_hist.C("data/mockdata_michel/run1_A_plus.dat")'
//   root -l -q 'macros/plot_michel_e_hist.C("data/mockdata_michel/run1_A_plus.dat", 0, 60, 120)'
//

R__ADD_INCLUDE_PATH(./include)

#include <string>
#include <vector>
#include <cmath>

#include "TCanvas.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"
#include "TLatex.h"
#include "TLine.h"
#include "TError.h"

#include "p2meg/MichelEData.h"
#include "p2meg/MichelPolConfig.h"

// 未解決を避けるため実装を同一モジュールに取り込む（ROOTマクロ用）
#include "../src/MichelEData.cc"

static TString MakeOutputPdfPath(const char* infile)
{
    TString base = gSystem->BaseName(infile); // e.g. xxx.dat
    Ssiz_t dot = base.Last('.');
    if (dot != kNPOS) base.Remove(dot);       // e.g. xxx

    gSystem->mkdir("doc", /*recursive=*/true);
    return Form("doc/michel_e_hist_%s.pdf", base.Data());
}

static void DrawMetaPage(const char* infile,
                         const char* outpdf,
                         long long n_read,
                         long long n_skipped,
                         long long n_inrange,
                         long long n_under,
                         long long n_over,
                         double Ee_min,
                         double Ee_max,
                         int nbins,
                         double mean_in,
                         double rms_in,
                         bool file_missing)
{
    TLatex lat;
    lat.SetNDC(true);
    lat.SetTextAlign(13);

    lat.SetTextSize(0.040);
    lat.DrawLatex(0.05, 0.92, "Michel Ee histogram (p2MEG)");

    lat.SetTextSize(0.030);
    lat.DrawLatex(0.05, 0.85, Form("input  : %s", infile));
    lat.DrawLatex(0.05, 0.81, Form("output : %s", outpdf));

    lat.SetTextSize(0.028);
    lat.DrawLatex(0.05, 0.74, "I/O summary:");
    lat.DrawLatex(0.08, 0.70, Form("read    = %lld", n_read));
    lat.DrawLatex(0.08, 0.66, Form("skipped = %lld", n_skipped));

    lat.DrawLatex(0.05, 0.58, "hist range:");
    lat.DrawLatex(0.08, 0.54, Form("Ee [MeV] = [%.6g, %.6g], nbins = %d", Ee_min, Ee_max, nbins));

    lat.DrawLatex(0.05, 0.46, "counts (in range):");
    lat.DrawLatex(0.08, 0.42, Form("in range  = %lld", n_inrange));
    lat.DrawLatex(0.08, 0.38, Form("underflow = %lld", n_under));
    lat.DrawLatex(0.08, 0.34, Form("overflow  = %lld", n_over));

    lat.DrawLatex(0.05, 0.26, "stats (in range):");
    if (n_inrange > 0) {
        lat.DrawLatex(0.08, 0.22, Form("mean = %.6g MeV", mean_in));
        lat.DrawLatex(0.08, 0.18, Form("RMS  = %.6g MeV", rms_in));
    } else {
        lat.DrawLatex(0.08, 0.22, "mean = (no entries)");
        lat.DrawLatex(0.08, 0.18, "RMS  = (no entries)");
    }

    if (file_missing) {
        lat.SetTextSize(0.032);
        lat.DrawLatex(0.05, 0.08, "warning: input file not found (no data)");
    }
}

void plot_michel_e_hist(const char* path,
                        double Ee_min = michel_pol_config.Ee_min,
                        double Ee_max = michel_pol_config.Ee_max,
                        int nbins = michel_pol_config.nbins_Ee)
{
    gStyle->SetOptStat(0);

    long long n_read = 0;
    long long n_skipped = 0;
    const auto evs = ReadMichelEData(path, &n_read, &n_skipped);

    const bool file_missing = gSystem->AccessPathName(path);

    const TString outpdf = MakeOutputPdfPath(path);

    TH1D hEe("hEe", "Michel Ee spectrum;Ee [MeV];Counts", nbins, Ee_min, Ee_max);

    long long n_inrange = 0;
    long long n_under = 0;
    long long n_over = 0;
    for (const auto& ev : evs) {
        if (ev.Ee < Ee_min) {
            n_under++;
        } else if (ev.Ee >= Ee_max) {
            n_over++;
        } else {
            hEe.Fill(ev.Ee);
            n_inrange++;
        }
    }

    const double mean_in = (n_inrange > 0) ? hEe.GetMean() : 0.0;
    const double rms_in  = (n_inrange > 0) ? hEe.GetRMS()  : 0.0;

    TCanvas c0("c0", "meta", 900, 700);
    DrawMetaPage(path, outpdf.Data(), n_read, n_skipped,
                 n_inrange, n_under, n_over,
                 Ee_min, Ee_max, nbins,
                 mean_in, rms_in, file_missing);

    TCanvas c1("c1", "hist", 900, 700);
    hEe.SetLineWidth(2);
    hEe.Draw("hist");

    // fit 範囲（MichelPolConfig）を補助線として描画
    const double fit_min = michel_pol_config.fit_Ee_min;
    const double fit_max = michel_pol_config.fit_Ee_max;
    double ymax = hEe.GetMaximum();
    if (!(ymax > 0.0)) ymax = 1.0;
    if (fit_min > Ee_min && fit_min < Ee_max) {
        TLine lmin(fit_min, 0.0, fit_min, ymax * 1.05);
        lmin.SetLineColor(kRed + 1);
        lmin.SetLineStyle(2);
        lmin.Draw("same");
    }
    if (fit_max > Ee_min && fit_max < Ee_max) {
        TLine lmax(fit_max, 0.0, fit_max, ymax * 1.05);
        lmax.SetLineColor(kRed + 1);
        lmax.SetLineStyle(2);
        lmax.Draw("same");
    }

    TLatex lat;
    lat.SetNDC(true);
    lat.SetTextSize(0.030);
    lat.SetTextAlign(13);
    lat.DrawLatex(0.12, 0.86, "dashed lines: fit range (MichelPolConfig)");

    c0.Print(Form("%s[", outpdf.Data()));
    c0.Print(outpdf.Data());
    c1.Print(outpdf.Data());
    c1.Print(Form("%s]", outpdf.Data()));
}
