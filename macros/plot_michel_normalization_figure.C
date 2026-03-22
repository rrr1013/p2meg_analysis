// macros/plot_michel_normalization_figure.C
//
// 目的:
//   正規化節に入れるための Michel 陽電子エネルギー図を 1 ページ PDF で作る。
//   実データの run-energy 2 列ファイルを読み、E>=20 MeV の閾値を明示する。
//
// 入力:
//   run energy
//   の 2 列テキスト（先頭ヘッダ許容）
//
// 出力:
//   doc/fig/fig06_michel_normalization.pdf
//   doc/fig/fig06_michel_normalization.png

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TH1D.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"

struct RunEnergyRow {
    int run;
    double energy; // [MeV]
};

static bool ParseRunEnergyRow(const std::string& line, RunEnergyRow& row)
{
    std::istringstream iss(line);
    if (!(iss >> row.run >> row.energy)) return false;
    return true;
}

void plot_michel_normalization_figure(
    const char* infile =
        "data/finaldata/positron_run_energy_before_eg_match_run1to20_doubleonly.txt")
{
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    std::ifstream fin(infile);
    if (!fin) {
        ::Error("plot_michel_normalization_figure", "cannot open %s", infile);
        return;
    }

    std::vector<RunEnergyRow> rows;
    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') continue;
        RunEnergyRow row{};
        if (!ParseRunEnergyRow(line, row)) continue;
        rows.push_back(row);
    }

    TH1D h("h_michel_norm",
           ";E_{e^{+}} [MeV];Counts",
           120, 0.0, 60.0);
    for (std::size_t i = 0; i < rows.size(); ++i) {
        h.Fill(rows[i].energy);
    }

    gSystem->mkdir("doc/fig", /*recursive=*/true);

    TCanvas c("c_michel_norm", "Michel normalization", 920, 680);
    c.SetLeftMargin(0.12);
    c.SetBottomMargin(0.12);
    c.SetRightMargin(0.04);
    c.SetTopMargin(0.10);

    h.SetLineColor(kBlue + 1);
    h.SetLineWidth(2);
    h.Draw("hist");

    const double ymax = (h.GetMaximum() > 0.0) ? h.GetMaximum() : 1.0;

    TLine lcut(20.0, 0.0, 20.0, 1.05 * ymax);
    lcut.SetLineColor(kRed + 1);
    lcut.SetLineStyle(2);
    lcut.SetLineWidth(3);
    lcut.Draw("same");

    TLatex lat;
    lat.SetNDC(true);
    lat.SetTextFont(42);
    lat.SetTextSize(0.035);
    lat.SetTextSize(0.032);
    lat.DrawLatex(0.58, 0.86, "normalization region");
    lat.DrawLatex(0.58, 0.81, "E_{e^{+}} #geq 20 MeV");

    TLegend leg(0.58, 0.66, 0.90, 0.76);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.030);
    leg.AddEntry(&h, "data", "l");
    leg.AddEntry(&lcut, "threshold at 20 MeV", "l");
    leg.Draw();

    c.SaveAs("doc/fig/fig06_michel_normalization.pdf");
    c.SaveAs("doc/fig/fig06_michel_normalization.png");
}
