// macros/plot_michelsim_epos_smeared_hist.C
//
// Michelsim 6列テキストの Epos 列に
// DetectorResolution.h の e+ エネルギー応答を乱数でかけ、
// run_energy_hist と同系統の 2 ページ PDF を作る。
//
// 入力:
//   - 6列: Epos Egam dt phi_pos phi_gam angle_deg
//   - 先頭ヘッダ、空行、'#' 始まり行は無視
//
// 出力:
//   - doc/finalanalysis/run_energy_hist_<basename>_epos_smeared.pdf
//     1) メタ情報
//     2) スメア後 Epos ヒストグラム
//
// 実行例:
//   root -l -b -q 'macros/plot_michelsim_epos_smeared_hist.C("data/finaldata/Michelsim_5e5.txt")'

R__ADD_INCLUDE_PATH(./include)

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TRandom3.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"

#include "p2meg/DetectorResolution.h"

struct MichelSimRow {
    double Epos;      // [MeV]
    double Egam;      // [MeV]
    double dt;        // [ns]
    double phi_pos;   // [rad]
    double phi_gam;   // [rad]
    double angle_deg; // [deg]
};

static TString MakeOutputPdfPathMichelSim(const char* infile)
{
    TString base = gSystem->BaseName(infile);
    Ssiz_t dot = base.Last('.');
    if (dot != kNPOS) base.Remove(dot);

    gSystem->mkdir("doc/finalanalysis", /*recursive=*/true);
    return Form("doc/finalanalysis/run_energy_hist_%s_epos_smeared.pdf", base.Data());
}

static bool ParseMichelSimLine(const std::string& line, MichelSimRow& row)
{
    std::istringstream iss(line);
    if (!(iss >> row.Epos >> row.Egam >> row.dt
              >> row.phi_pos >> row.phi_gam >> row.angle_deg)) {
        return false;
    }

    if (!std::isfinite(row.Epos) || !std::isfinite(row.Egam) ||
        !std::isfinite(row.dt) || !std::isfinite(row.phi_pos) ||
        !std::isfinite(row.phi_gam) || !std::isfinite(row.angle_deg)) {
        return false;
    }
    return true;
}

static void ReadMichelSimRows(const char* path,
                              std::vector<MichelSimRow>& rows,
                              long long& n_read,
                              long long& n_skipped)
{
    rows.clear();
    n_read = 0;
    n_skipped = 0;

    std::ifstream fin(path);
    if (!fin) return;

    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') continue;

        MichelSimRow row{};
        if (!ParseMichelSimLine(line, row)) {
            ++n_skipped;
            continue;
        }

        rows.push_back(row);
        ++n_read;
    }
}

static void DrawMichelSimMetaPage(const char* infile,
                                  const char* outpdf,
                                  long long n_read,
                                  long long n_skipped,
                                  long long n_ge20,
                                  double e_min,
                                  double e_max,
                                  int nbins,
                                  unsigned int seed,
                                  double mean,
                                  double rms,
                                  bool file_missing)
{
    TLatex lat;
    lat.SetNDC(true);
    lat.SetTextAlign(13);

    lat.SetTextSize(0.040);
    lat.DrawLatex(0.05, 0.92, "Run energy histogram from Michelsim Epos smear");

    lat.SetTextSize(0.030);
    lat.DrawLatex(0.05, 0.85, Form("input  : %s", infile));
    lat.DrawLatex(0.05, 0.81, Form("output : %s", outpdf));

    lat.SetTextSize(0.028);
    lat.DrawLatex(0.05, 0.73, "I/O summary:");
    lat.DrawLatex(0.08, 0.69, Form("read    = %lld", n_read));
    lat.DrawLatex(0.08, 0.65, Form("skipped = %lld", n_skipped));
    lat.DrawLatex(0.08, 0.61, Form("E >= 20 MeV = %lld", n_ge20));

    lat.DrawLatex(0.05, 0.55, "Smear model:");
    lat.DrawLatex(0.08, 0.51, Form("seed = %u", seed));
    lat.DrawLatex(0.08, 0.47, "Ee response = Gaussian");
    lat.DrawLatex(0.08, 0.43, "sigma(Ee) = 0.079 x Etrue");

    lat.DrawLatex(0.05, 0.37, "Histogram:");
    lat.DrawLatex(0.08, 0.33, Form("range = %.3f .. %.3f MeV", e_min, e_max));
    lat.DrawLatex(0.08, 0.29, Form("nbins = %d", nbins));
    lat.DrawLatex(0.08, 0.25, Form("mean  = %.6f MeV", mean));
    lat.DrawLatex(0.08, 0.21, Form("RMS   = %.6f MeV", rms));

    lat.SetTextSize(0.024);
    lat.DrawLatex(0.05, 0.12, "note: non-numeric/header lines are skipped");
    lat.DrawLatex(0.05, 0.09, "note: smear uses smear_energy_trandom3_e in DetectorResolution.h");

    if (file_missing) {
        lat.SetTextColor(kRed + 1);
        lat.DrawLatex(0.05, 0.04, "warning: input file not found (no data)");
    }
}

void plot_michelsim_epos_smeared_hist(const char* path,
                                      double e_min = 0.0,
                                      double e_max = 60.0,
                                      int nbins = 120,
                                      unsigned int seed = 12345u)
{
    gStyle->SetOptStat(0);

    std::vector<MichelSimRow> rows;
    long long n_read = 0;
    long long n_skipped = 0;
    ReadMichelSimRows(path, rows, n_read, n_skipped);

    const bool file_missing = gSystem->AccessPathName(path);
    const TString outpdf = MakeOutputPdfPathMichelSim(path);

    TRandom3 rng(seed);

    TH1D hE("hE", "Positron energy from Michelsim with detector smear;Energy [MeV];Counts",
            nbins, e_min, e_max);

    long long n_ge20 = 0;
    for (size_t i = 0; i < rows.size(); ++i) {
        const double Epos_smeared = smear_energy_trandom3_e(rng, rows[i].Epos);
        if (Epos_smeared >= 20.0) ++n_ge20;
        hE.Fill(Epos_smeared);
    }

    const double mean = (hE.GetEntries() > 0.0) ? hE.GetMean() : 0.0;
    const double rms  = (hE.GetEntries() > 0.0) ? hE.GetRMS()  : 0.0;

    TCanvas c0("c0", "meta", 900, 700);
    DrawMichelSimMetaPage(path, outpdf.Data(),
                          n_read, n_skipped, n_ge20,
                          e_min, e_max, nbins, seed,
                          mean, rms, file_missing);

    TCanvas c1("c1", "hist", 900, 700);
    hE.SetLineWidth(2);
    hE.Draw("hist");

    c0.Print(Form("%s[", outpdf.Data()));
    c0.Print(outpdf.Data());
    c1.Print(outpdf.Data());
    c1.Print(Form("%s]", outpdf.Data()));

    std::cout << "[plot_michelsim_epos_smeared_hist] input    : " << path << "\n";
    std::cout << "[plot_michelsim_epos_smeared_hist] output   : " << outpdf.Data() << "\n";
    std::cout << "[plot_michelsim_epos_smeared_hist] seed     : " << seed << "\n";
    std::cout << "[plot_michelsim_epos_smeared_hist] n_read   : " << n_read << "\n";
    std::cout << "[plot_michelsim_epos_smeared_hist] n_ge20   : " << n_ge20 << "\n";
    std::cout << "[plot_michelsim_epos_smeared_hist] mean/rms : "
              << mean << " / " << rms << "\n";
}
