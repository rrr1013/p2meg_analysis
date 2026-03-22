#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "TArrow.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TBox.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
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

double ExtractValue(const std::string& line, const char* key)
{
    if (!StartsWith(line, key)) return 0.0;
    std::istringstream iss(line);
    std::string k, eq;
    double value = 0.0;
    iss >> k >> eq >> value;
    return value;
}

double AutoBinWidth(const std::vector<double>& values)
{
    std::vector<double> uniq = values;
    std::sort(uniq.begin(), uniq.end());
    uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());
    if (uniq.size() < 2U) return 1.0e-4;

    double dx_min = 0.0;
    for (std::size_t i = 1; i < uniq.size(); ++i) {
        const double dx = uniq[i] - uniq[i - 1U];
        if (!(dx > 0.0) || !std::isfinite(dx)) continue;
        if (!(dx_min > 0.0) || dx < dx_min) dx_min = dx;
    }
    if (!(dx_min > 0.0) || !std::isfinite(dx_min)) return 1.0e-4;
    if (dx_min < 1.0e-4) return 1.0e-4;
    return dx_min;
}

}

void plot_sensitivity_br_hist(
    const char* infile = "doc/finalanalysis/sensitivity_br_scan_blind50ns_broad_100000.txt",
    const char* outfile = "doc/finalanalysis/sensitivity_br_hist_blind50ns.pdf",
    double observed_br90_override = -1.0,
    double bin_width = 0.0)
{
    std::ifstream fin(infile);
    if (!fin) {
        Error("plot_sensitivity_br_hist", "cannot open %s", infile);
        return;
    }

    double observed_br90 = -1.0;
    double br50 = -1.0;
    double br16 = -1.0;
    double br84 = -1.0;
    std::vector<double> br_values;
    bool in_toy_block = false;

    std::string line;
    while (std::getline(fin, line)) {
        if (StartsWith(line, "observed_local_exact_BR90")) {
            observed_br90 = ExtractValue(line, "observed_local_exact_BR90");
            continue;
        }

        if (StartsWith(line, "observed_fast_BR90")) {
            if (!(observed_br90 >= 0.0)) {
                observed_br90 = ExtractValue(line, "observed_fast_BR90");
            }
            continue;
        }

        if (StartsWith(line, "sensitivity_BR50")) {
            br50 = ExtractValue(line, "sensitivity_BR50");
            continue;
        }
        if (StartsWith(line, "sensitivity_BR16")) {
            br16 = ExtractValue(line, "sensitivity_BR16");
            continue;
        }
        if (StartsWith(line, "sensitivity_BR84")) {
            br84 = ExtractValue(line, "sensitivity_BR84");
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
    if (!(step > 0.0) || !std::isfinite(step)) step = AutoBinWidth(br_values);
    xmin -= 0.5 * step;
    xmax += 0.5 * step;
    if (xmin < 0.0) xmin = 0.0;
    if (xmax <= xmin) xmax = xmin + step;

    const int nbins = static_cast<int>((xmax - xmin) / step + 0.5);
    TH1D* h = new TH1D("h_sensitivity_br",
                       ";90% C.L. upper limit on BR(#mu #rightarrow e#gamma);Pseudo-experiments",
                       nbins, xmin, xmax);

    for (double v : br_values) h->Fill(v);

    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTextFont(42);

    TCanvas* c = new TCanvas("c_sensitivity_br", "c_sensitivity_br", 860, 620);
    c->SetMargin(0.145, 0.04, 0.12, 0.06);

    const Color_t hist_line = TColor::GetColor("#264653");
    const Color_t hist_fill = TColor::GetColor("#8fb7c9");
    const Color_t band_fill = TColor::GetColor("#8ecf8a");

    h->SetLineColor(hist_line);
    h->SetFillColor(hist_fill);
    h->SetLineWidth(2);
    h->SetMaximum(h->GetMaximum() * 1.20);
    h->GetXaxis()->SetTitleOffset(1.05);
    h->GetYaxis()->SetTitleOffset(1.35);

    TBox* band_box = nullptr;
    if (br16 >= 0.0 && br84 >= 0.0 && br84 >= br16) {
        band_box = new TBox(br16, 0.0, br84, h->GetMaximum());
        band_box->SetFillColorAlpha(band_fill, 0.16);
        band_box->SetLineColor(kGreen + 2);
        band_box->SetLineStyle(1);
        band_box->SetLineWidth(1);
    }

    h->Draw("hist");
    if (band_box) band_box->Draw("same");
    h->Draw("hist same");

    TLine* med_line = nullptr;
    if (br50 >= 0.0) {
        med_line = new TLine(br50, 0.0, br50, h->GetMaximum() * 0.92);
        med_line->SetLineColor(kBlack);
        med_line->SetLineStyle(1);
        med_line->SetLineWidth(3);
        med_line->Draw();
    }

    TLine* obs_line = new TLine(observed_br90, 0.0, observed_br90, h->GetMaximum() * 0.92);
    obs_line->SetLineColor(kRed + 1);
    obs_line->SetLineStyle(2);
    obs_line->SetLineWidth(3);
    obs_line->Draw();

    TLegend* leg = new TLegend(0.55, 0.76, 0.89, 0.89);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.030);
    leg->AddEntry(h, "Background-only pseudo-experiments", "f");
    if (band_box) leg->AddEntry(band_box, "Central 68% interval", "f");
    if (med_line) leg->AddEntry(med_line, Form("Median = %.5f", br50), "l");
    leg->AddEntry(obs_line, Form("Observed = %.5f", observed_br90), "l");
    leg->Draw();

    c->SaveAs(outfile);
}
