// macros/plot_fit_acc_time_by_eg_thesis.C
//
// 目的:
//   ACC 時間形状 fit の Eg bin ごとの 4 枚の t 分布を、
//   論文用の 2x2 レイアウトで描く。
//
// 入力:
//   5列テキスト
//     Ee Eg t phi_detector_e phi_detector_g
//
// 選別:
//   fit_acc_time_by_eg.C と同じ。
//   - phi は許可ペア
//   - Ee は解析窓内
//   - theta_eg は解析窓内
//   - Eg は bin ごと
//   - t は [-500, 500] ns
//   - fit ヒストには 4D AW 内イベントを入れない
//
// 出力:
//   doc/fig/fit_acc_time_by_eg_<input>_blind50ns_sigma30ns_thesis_2x2.pdf
//   doc/fig/fit_acc_time_by_eg_<input>_blind50ns_sigma30ns_thesis_2x2.png
//
// 実行例:
//   root -l -b -q 'macros/plot_fit_acc_time_by_eg_thesis.C()'

R__ADD_INCLUDE_PATH(./include)

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

#include "TBox.h"
#include "TCanvas.h"
#include "TError.h"
#include "TF1.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TPave.h"
#include "TPaveText.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"

#include "p2meg/AccTimeFit.h"
#include "p2meg/AnalysisWindow.h"
#include "p2meg/DetectorResolution.h"

static constexpr double kTAllMin_thesis = -500.0; // [ns]
static constexpr double kTAllMax_thesis =  500.0; // [ns]
static constexpr int    kNBinsT_thesis  = 100;
static const double kEgEdges_thesis[] = {0.0, 10.0, 20.0, 30.0, 80.0};
static constexpr int kNEgBins_thesis =
    static_cast<int>(sizeof(kEgEdges_thesis) / sizeof(kEgEdges_thesis[0])) - 1;

static bool ParseEventLine5Doubles_thesis(const std::string& line,
                                          double& Ee, double& Eg, double& t,
                                          double& phi_e, double& phi_g)
{
    if (line.empty()) return false;
    std::istringstream iss(line);
    if (!(iss >> Ee >> Eg >> t >> phi_e >> phi_g)) return false;
    if (!std::isfinite(Ee) || !std::isfinite(Eg) || !std::isfinite(t) ||
        !std::isfinite(phi_e) || !std::isfinite(phi_g)) {
        return false;
    }
    return true;
}

static int FindEgBin_thesis(double Eg)
{
    if (!std::isfinite(Eg)) return -1;
    for (int i = 0; i < kNEgBins_thesis; ++i) {
        const double lo = kEgEdges_thesis[i];
        const double hi = kEgEdges_thesis[i + 1];
        if ((Eg >= lo && Eg < hi) || (i == kNEgBins_thesis - 1 && Eg >= lo && Eg <= hi)) {
            return i;
        }
    }
    return -1;
}

static bool IsInsideFullAnalysisWindow_thesis(double Ee, double Eg,
                                              double theta_eg, double t)
{
    if (!std::isfinite(Ee) || !std::isfinite(Eg) ||
        !std::isfinite(theta_eg) || !std::isfinite(t)) {
        return false;
    }
    if (Ee < analysis_window.Ee_min || Ee > analysis_window.Ee_max) return false;
    if (Eg < analysis_window.Eg_min || Eg > analysis_window.Eg_max) return false;
    if (theta_eg < analysis_window.theta_min || theta_eg > analysis_window.theta_max) return false;
    if (t < analysis_window.t_min || t > analysis_window.t_max) return false;
    return true;
}

static void ConvertCountsToDensity_thesis(TH1D& h)
{
    for (int ib = 1; ib <= h.GetNbinsX(); ++ib) {
        const double w = h.GetXaxis()->GetBinWidth(ib);
        if (!(w > 0.0) || !std::isfinite(w)) {
            h.SetBinContent(ib, 0.0);
            h.SetBinError(ib, 0.0);
            continue;
        }
        h.SetBinContent(ib, h.GetBinContent(ib) / w);
        h.SetBinError(ib, h.GetBinError(ib) / w);
    }
}

static void StyleHist_thesis(TH1D& h)
{
    h.SetLineColor(kBlack);
    h.SetLineWidth(1);
    h.SetMarkerStyle(20);
    h.SetMarkerSize(0.45);
    h.GetXaxis()->SetTitleSize(0.070);
    h.GetYaxis()->SetTitleSize(0.070);
    h.GetXaxis()->SetLabelSize(0.055);
    h.GetYaxis()->SetLabelSize(0.055);
    h.GetXaxis()->SetTitleOffset(0.95);
    h.GetYaxis()->SetTitleOffset(1.05);
}

static TString MakeOutputBase_thesis(const char* infile,
                                     double blind_core_ns,
                                     double sigma_min_ns)
{
    TString base = gSystem->BaseName(infile);
    const Ssiz_t dot = base.Last('.');
    if (dot != kNPOS) base.Remove(dot);
    gSystem->mkdir("doc/fig", /*recursive=*/true);
    return Form("doc/fig/fit_acc_time_by_eg_%s_blind%.0fns_sigma%.0fns_thesis_2x2",
                base.Data(), blind_core_ns, sigma_min_ns);
}

static void DrawFitPanel_thesis(TH1D& h,
                                int eg_bin,
                                long long n_entries,
                                double blind_core_ns,
                                double sigma_min_ns)
{
    const double hmax = std::max(1.0, 1.20 * h.GetMaximum());

    AccTimeFitResult fitres{};
    AccTimeFit_Run(h, -blind_core_ns, blind_core_ns,
                   sigma_min_ns, 400.0, fitres);

    TF1 f_draw(Form("f_acc_time_draw_thesis_%d", eg_bin), AccTimeFit_DensityCore,
               kTAllMin_thesis, kTAllMax_thesis, 3);
    f_draw.SetParameters(fitres.A, fitres.sigma, fitres.C);
    f_draw.SetLineColor(kRed + 1);
    f_draw.SetLineWidth(1);

    gPad->SetGrid(1, 1);
    gPad->SetLeftMargin(0.16);
    gPad->SetBottomMargin(0.20);
    gPad->SetRightMargin(0.04);
    gPad->SetTopMargin(0.17);

    h.SetMaximum(hmax);
    h.SetMinimum(0.0);
    h.Draw("E");

    TBox blind_box(-blind_core_ns, 0.0, blind_core_ns, hmax);
    blind_box.SetFillColorAlpha(kGray + 1, 0.22);
    blind_box.SetLineColor(kGray + 2);
    blind_box.SetLineWidth(1);
    blind_box.DrawClone("same");

    f_draw.DrawCopy("same");

    TLatex title;
    title.SetNDC(true);
    title.SetTextAlign(22);
    title.SetTextSize(0.062);
    title.DrawLatex(0.50, 0.950,
                    Form("t, E_{#gamma} in [%.0f, %.0f]",
                         kEgEdges_thesis[eg_bin], kEgEdges_thesis[eg_bin + 1]));

    // メタ情報は白背景の箱に入れて、グリッド・ヒストと分離する。
    TPaveText* meta = new TPaveText(0.12, 0.70, 0.50, 0.90, "NDC");
    meta->SetFillColor(kWhite);
    meta->SetFillStyle(1001);
    meta->SetLineColor(kBlack);
    meta->SetLineWidth(1);
    meta->SetTextAlign(12);
    meta->SetTextFont(42);
    meta->SetTextSize(0.036);
    meta->SetMargin(0.03);
    meta->AddText(Form("E_{#gamma} in [%.0f, %.0f] MeV, entries=%lld",
                       kEgEdges_thesis[eg_bin], kEgEdges_thesis[eg_bin + 1], n_entries));
    meta->AddText(Form("fit status=%d, #chi^{2}/ndf=%.2f / %d",
                       fitres.fit_status, fitres.chi2, fitres.ndf));
    meta->AddText(Form("A=%.4g, #sigma=%.2f ns, c_{0}=%.4g",
                       fitres.A, fitres.sigma, fitres.C));
    meta->Draw();

    TPaveText* legend_box = new TPaveText(0.70, 0.67, 0.98, 0.90, "NDC");
    legend_box->SetFillColor(kWhite);
    legend_box->SetFillStyle(1001);
    legend_box->SetLineColor(kBlack);
    legend_box->SetLineWidth(1);
    legend_box->Draw();

    const double x1 = 0.73;
    const double x2 = 0.83;
    double y = 0.84;

    TLine data_line;
    data_line.SetLineColor(kBlack);
    data_line.SetLineWidth(1);
    data_line.DrawLineNDC(x1, y, x2, y);
    TLatex text;
    text.SetNDC(true);
    text.SetTextAlign(13);
    text.SetTextFont(42);
    text.SetTextSize(0.043);
    text.DrawLatex(0.84, y + 0.01, "data");

    y = 0.77;
    TLine fit_line;
    fit_line.SetLineColor(kRed + 1);
    fit_line.SetLineWidth(1);
    fit_line.DrawLineNDC(x1, y, x2, y);
    text.DrawLatex(0.84, y + 0.01, "fit");

    y = 0.70;
    TPave* blind_legend = new TPave(x1, y - 0.028, x2, y + 0.028, 0, "NDC");
    blind_legend->SetFillColor(kGray + 1);
    blind_legend->SetLineColor(kGray + 2);
    blind_legend->SetLineWidth(1);
    blind_legend->Draw();
    text.DrawLatex(0.86, y + 0.01, "blind core");
}

void plot_fit_acc_time_by_eg_thesis(
    const char* infile = "data/finaldata/step4f_run1to20_eg_doubleonly_allpatterns_500ns.txt",
    double blind_core_ns = 50.0,
    double sigma_min_ns = 30.0)
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetGridStyle(3);
    gStyle->SetGridColor(kGray + 1);

    TH1D* hT[kNEgBins_thesis];
    long long n_fit_eg_bin[kNEgBins_thesis];
    for (int i = 0; i < kNEgBins_thesis; ++i) {
        hT[i] = new TH1D(Form("hT_thesis_%d", i),
                         ";t_{e^{+}#gamma} [ns];density [counts/ns]",
                         kNBinsT_thesis, kTAllMin_thesis, kTAllMax_thesis);
        hT[i]->Sumw2();
        n_fit_eg_bin[i] = 0;
    }

    std::ifstream fin(infile);
    if (!fin) {
        Error("plot_fit_acc_time_by_eg_thesis", "failed to open input file: %s", infile);
        return;
    }

    std::string line;
    while (std::getline(fin, line)) {
        double Ee = 0.0;
        double Eg = 0.0;
        double t = 0.0;
        double phi_e = 0.0;
        double phi_g = 0.0;
        if (!ParseEventLine5Doubles_thesis(line, Ee, Eg, t, phi_e, phi_g)) continue;
        if (!std::isfinite(t) || t < kTAllMin_thesis || t > kTAllMax_thesis) continue;

        int idx_e = -1;
        int idx_g = -1;
        if (!Detector_IsAllowedPhiPairValue(phi_e, phi_g, detres, idx_e, idx_g)) continue;

        const double phi_e_disc =
            Detector_PhiGridPoint(idx_e, detres.phi_e_min, detres.phi_e_max, detres.N_phi_e);
        const double phi_g_disc =
            Detector_PhiGridPoint(idx_g, detres.phi_g_min, detres.phi_g_max, detres.N_phi_g);
        const double theta_eg = std::fabs(phi_e_disc - phi_g_disc);

        if (Ee < analysis_window.Ee_min || Ee > analysis_window.Ee_max) continue;
        if (theta_eg < analysis_window.theta_min || theta_eg > analysis_window.theta_max) continue;

        const int ib = FindEgBin_thesis(Eg);
        if (ib < 0) continue;

        if (!IsInsideFullAnalysisWindow_thesis(Ee, Eg, theta_eg, t)) {
            hT[ib]->Fill(t);
            ++n_fit_eg_bin[ib];
        }
    }

    for (int i = 0; i < kNEgBins_thesis; ++i) {
        ConvertCountsToDensity_thesis(*hT[i]);
        StyleHist_thesis(*hT[i]);
    }

    const TString outbase = MakeOutputBase_thesis(infile, blind_core_ns, sigma_min_ns);
    const TString outpdf = outbase + ".pdf";
    const TString outpng = outbase + ".png";

    TCanvas c("c_fit_acc_time_thesis", "fit acc time thesis", 1800, 980);
    c.Divide(2, 2, 0.020, 0.030);
    for (int i = 0; i < kNEgBins_thesis; ++i) {
        c.cd(i + 1);
        DrawFitPanel_thesis(*hT[i], i, n_fit_eg_bin[i], blind_core_ns, sigma_min_ns);
    }

    c.Print(outpdf.Data());
    c.Print(outpng.Data());

    Info("plot_fit_acc_time_by_eg_thesis", "wrote %s", outpdf.Data());
    Info("plot_fit_acc_time_by_eg_thesis", "wrote %s", outpng.Data());
}
