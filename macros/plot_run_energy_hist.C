// macros/plot_run_energy_hist.C
//
// run-energy 2列テキストの energy 列ヒストグラムを作る
//
// 入力:
//   - 先頭1行はヘッダを許す
//   - 2列: run, energy
//   - energy の単位は MeV を想定
//   - 空行と '#' 始まり行は無視
//
// 出力:
//   - doc/finalanalysis/run_energy_hist_<basename>.pdf（2ページ）
//     1) メタ情報
//     2) energy ヒストグラム
//
// 実行例:
//   root -l -q 'macros/plot_run_energy_hist.C("data/finaldata/positron_run_energy_before_eg_match_run1to20_doubleonly.txt")'
//   root -l -q 'macros/plot_run_energy_hist.C("data/finaldata/positron_run_energy_before_eg_match_run1to20_doubleonly.txt", 0.0, 60.0, 120)'

#include <fstream>
#include <sstream>
#include <string>

#include "TCanvas.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"
#include "TLatex.h"

struct RunEnergyEvent {
    int run;
    double energy; // [MeV]
};

static TString MakeRunEnergyOutputPdfPath(const char* infile)
{
    TString base = gSystem->BaseName(infile);
    Ssiz_t dot = base.Last('.');
    if (dot != kNPOS) base.Remove(dot);

    gSystem->mkdir("doc/finalanalysis", /*recursive=*/true);
    return Form("doc/finalanalysis/run_energy_hist_%s.pdf", base.Data());
}

static bool ParseRunEnergyLine(const std::string& line, int& run, double& energy)
{
    std::istringstream iss(line);
    if (!(iss >> run >> energy)) return false;
    return true;
}

static void ReadRunEnergyData(const char* path,
                              std::vector<RunEnergyEvent>& events,
                              long long& n_read,
                              long long& n_skipped,
                              int& run_min,
                              int& run_max)
{
    n_read = 0;
    n_skipped = 0;
    run_min = 0;
    run_max = 0;

    std::ifstream fin(path);
    if (!fin) return;

    std::string line;
    bool has_run = false;
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') continue;

        int run = 0;
        double energy = 0.0;
        if (!ParseRunEnergyLine(line, run, energy)) {
            n_skipped++;
            continue;
        }

        RunEnergyEvent ev;
        ev.run = run;
        ev.energy = energy;
        events.push_back(ev);
        n_read++;

        if (!has_run) {
            run_min = run;
            run_max = run;
            has_run = true;
        } else {
            if (run < run_min) run_min = run;
            if (run > run_max) run_max = run;
        }
    }
}

static void DrawRunEnergyMetaPage(const char* infile,
                                  const char* outpdf,
                                  long long n_read,
                                  long long n_skipped,
                                  long long n_ge20,
                                  int run_min,
                                  int run_max,
                                  double e_min,
                                  double e_max,
                                  int nbins,
                                  double mean,
                                  double rms,
                                  bool file_missing)
{
    TLatex lat;
    lat.SetNDC(true);
    lat.SetTextAlign(13);

    lat.SetTextSize(0.040);
    lat.DrawLatex(0.05, 0.92, "Run energy histogram (p2MEG)");

    lat.SetTextSize(0.030);
    lat.DrawLatex(0.05, 0.85, Form("input  : %s", infile));
    lat.DrawLatex(0.05, 0.81, Form("output : %s", outpdf));

    lat.SetTextSize(0.028);
    lat.DrawLatex(0.05, 0.73, "I/O summary:");
    lat.DrawLatex(0.08, 0.69, Form("read    = %lld", n_read));
    lat.DrawLatex(0.08, 0.65, Form("skipped = %lld", n_skipped));
    lat.DrawLatex(0.08, 0.61, Form("E >= 20 MeV = %lld", n_ge20));

    if (n_read > 0) {
        lat.DrawLatex(0.05, 0.57, "Data summary:");
        lat.DrawLatex(0.08, 0.53, Form("run range   = %d .. %d", run_min, run_max));
        lat.DrawLatex(0.08, 0.49, Form("hist range  = %.3f .. %.3f MeV", e_min, e_max));
        lat.DrawLatex(0.08, 0.45, Form("nbins       = %d", nbins));
        lat.DrawLatex(0.08, 0.41, Form("mean        = %.6f MeV", mean));
        lat.DrawLatex(0.08, 0.37, Form("RMS         = %.6f MeV", rms));
    }

    lat.SetTextSize(0.024);
    lat.DrawLatex(0.05, 0.20, "note: header/non-numeric lines are skipped");

    if (file_missing) {
        lat.SetTextColor(kRed + 1);
        lat.DrawLatex(0.05, 0.10, "warning: input file not found (no data)");
    }
}

void plot_run_energy_hist(const char* path,
                          double e_min = 0.0,
                          double e_max = 60.0,
                          int nbins = 120)
{
    gStyle->SetOptStat(0);

    std::vector<RunEnergyEvent> events;
    long long n_read = 0;
    long long n_skipped = 0;
    int run_min = 0;
    int run_max = 0;
    ReadRunEnergyData(path, events, n_read, n_skipped, run_min, run_max);

    const bool file_missing = gSystem->AccessPathName(path);
    const TString outpdf = MakeRunEnergyOutputPdfPath(path);

    TH1D hE("hE", "Positron energy before e#gamma match;Energy [MeV];Counts",
            nbins, e_min, e_max);

    long long n_ge20 = 0;
    for (size_t i = 0; i < events.size(); ++i) {
        if (events[i].energy >= 20.0) n_ge20++;
        hE.Fill(events[i].energy);
    }

    const double mean = (hE.GetEntries() > 0.0) ? hE.GetMean() : 0.0;
    const double rms  = (hE.GetEntries() > 0.0) ? hE.GetRMS()  : 0.0;

    TCanvas c0("c0", "meta", 900, 700);
    DrawRunEnergyMetaPage(path, outpdf.Data(),
                          n_read, n_skipped, n_ge20,
                          run_min, run_max,
                          e_min, e_max, nbins,
                          mean, rms, file_missing);

    TCanvas c1("c1", "hist", 900, 700);
    hE.SetLineWidth(2);
    hE.Draw("hist");

    c0.Print(Form("%s[", outpdf.Data()));
    c0.Print(outpdf.Data());
    c1.Print(outpdf.Data());
    c1.Print(Form("%s]", outpdf.Data()));
}
