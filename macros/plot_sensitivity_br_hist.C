#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "TArrow.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLine.h"
#include "TROOT.h"
#include "TStyle.h"

// ============================================================
// 感度分布ヒストグラム
//
// 入力:
//   run_sensitivity_br_toymc の出力 txt
//
// 出力:
//   doc/finalanalysis/sensitivity_br_hist_blind50ns.pdf
//
// observed BR upper limit を赤い矢印で示す。
// ============================================================

namespace {

bool StartsWith(const std::string& s, const char* prefix)
{
    const std::string p(prefix);
    return s.rfind(p, 0) == 0;
}

}

void plot_sensitivity_br_hist(
    const char* infile = "doc/finalanalysis/sensitivity_br_scan_blind50ns.txt",
    const char* outfile = "doc/finalanalysis/sensitivity_br_hist_blind50ns.pdf",
    double observed_br90_override = 8.5e-4,
    double bin_width = 1.0e-4)
{
    std::ifstream fin(infile);
    if (!fin) {
        Error("plot_sensitivity_br_hist", "cannot open %s", infile);
        return;
    }

    double observed_br90 = -1.0;
    std::vector<double> br_values;
    bool in_toy_block = false;

    std::string line;
    while (std::getline(fin, line)) {
        if (StartsWith(line, "observed_fast_BR90")) {
            std::istringstream iss(line);
            std::string key, eq;
            iss >> key >> eq >> observed_br90;
            continue;
        }

        if (line.find("=============== Sensitivity Toys") != std::string::npos) {
            in_toy_block = true;
            continue;
        }
        if (in_toy_block && line.find("===============================================") != std::string::npos) {
            in_toy_block = false;
            continue;
        }
        if (!in_toy_block) continue;
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        int toy_index = -1;
        double br90 = 0.0;
        if (!(iss >> toy_index >> br90)) continue;
        br_values.push_back(br90);
    }

    if (br_values.empty()) {
        Error("plot_sensitivity_br_hist", "no toy BR values found in %s", infile);
        return;
    }
    if (observed_br90_override > 0.0 && std::isfinite(observed_br90_override)) {
        observed_br90 = observed_br90_override;
    }
    if (!(observed_br90 >= 0.0)) {
        Error("plot_sensitivity_br_hist", "observed_fast_BR90 not found in %s", infile);
        return;
    }

    double xmin = br_values.front();
    double xmax = br_values.front();
    for (double v : br_values) {
        if (v < xmin) xmin = v;
        if (v > xmax) xmax = v;
    }

    double step = bin_width;
    if (!(step > 0.0) || !std::isfinite(step)) step = 1.0e-4;
    xmin -= 0.5 * step;
    xmax += 0.5 * step;
    if (xmin < 0.0) xmin = 0.0;
    if (xmax <= xmin) xmax = xmin + step;

    const int nbins = static_cast<int>((xmax - xmin) / step + 0.5);
    TH1D* h = new TH1D("h_sensitivity_br",
                       ";90% C.L. upper limit on BR(#mu #rightarrow e#gamma);Toy count",
                       nbins, xmin, xmax);

    for (double v : br_values) h->Fill(v);

    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTextFont(42);

    TCanvas* c = new TCanvas("c_sensitivity_br", "c_sensitivity_br", 900, 650);
    c->SetMargin(0.12, 0.04, 0.12, 0.08);

    h->SetLineColor(TColor::GetColor("#24445c"));
    h->SetFillColor(TColor::GetColor("#8fb7c9"));
    h->SetLineWidth(2);
    h->SetMaximum(h->GetMaximum() * 1.35);
    h->Draw("hist");

    TLine* obs_line = new TLine(observed_br90, 0.0, observed_br90, h->GetMaximum() * 0.78);
    obs_line->SetLineColor(kRed + 1);
    obs_line->SetLineStyle(2);
    obs_line->SetLineWidth(3);
    obs_line->Draw();

    TArrow* obs_arrow = new TArrow(observed_br90, h->GetMaximum() * 0.95,
                                   observed_br90, h->GetMaximum() * 0.80,
                                   0.02, "|>");
    obs_arrow->SetLineColor(kRed + 1);
    obs_arrow->SetFillColor(kRed + 1);
    obs_arrow->SetLineWidth(3);
    obs_arrow->Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.036);
    latex.DrawLatex(0.14, 0.93, "p2MEG expected upper-limit distribution");
    latex.SetTextSize(0.032);
    latex.DrawLatex(0.56, 0.88, Form("Observed BR_{90} = %.2g", observed_br90));
    latex.DrawLatex(0.56, 0.83, Form("Sensitivity toys = %zu", br_values.size()));
    latex.DrawLatex(0.56, 0.78, Form("Display bin width = %.0e", step));

    c->SaveAs(outfile);
}
