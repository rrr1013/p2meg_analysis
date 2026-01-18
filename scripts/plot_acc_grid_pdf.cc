// scripts/plot_acc_grid_pdf.cc
//
// acc_grid.root に保存された ACC 4D 格子 PDF の形状を
// plot_data_hist と同様の 1D/2D ヒスト形式で PDF 出力する。
//
// 使い方:
//   ./build/plot_acc_grid_pdf [in.root] [key] [out.pdf]
//
// 例:
//   ./build/plot_acc_grid_pdf data/pdf_cache/acc_grid.root acc_grid doc/acc_grid_hist.pdf
//
// 出力:
//   doc/acc_grid_hist.pdf（3ページ）
//    1) メタ情報
//    2) 1D (Ee, Eg, phi_detector_e, phi_detector_g, theta_eg)
//    3) 2D (Ee,Eg), (theta_eg,Ee), (theta_eg,Eg), (phi_e,phi_g)
//
// 注意:
//  - acc_grid は Ee/Eg の密度（MeV^-2）として格納されている。
//  - 表示用ヒストは dEe dEg を掛けた「確率質量」で作る。
//  - phi は離散点（N_theta+1）であり、phi の幅は正規化に含めない。
//  - t は解析窓内一様であり、本図では表示しない。

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THn.h"
#include "TAxis.h"
#include "TNamed.h"
#include "TParameter.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TError.h"
#include "TString.h"
#include "TLatex.h"

#include "p2meg/AnalysisWindow.h"
#include "p2meg/Constants.h"

// ---- 固定ビン数（plot_data_hist に合わせる）----
static constexpr int kNBins_E   = 120; // Ee, Eg
static constexpr int kNBins_phi = 120; // phi_detector_e/g
static constexpr int kNBins_th  = 120; // theta_eg

static constexpr int kNBins2D_E   = 120;
static constexpr int kNBins2D_th  = 120;
static constexpr int kNBins2D_phi = 120;

static bool IsFinite(double x) { return std::isfinite(x); }

static std::vector<std::string> SplitLines(const std::string& s)
{
    std::vector<std::string> lines;
    std::istringstream iss(s);
    std::string line;
    while (std::getline(iss, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        lines.push_back(line);
    }
    return lines;
}

static void DrawMetaPage(const char* infile,
                         const char* key,
                         const char* outpdf,
                         const THnD& h,
                         double sum_mass,
                         long n_negative,
                         long n_nonfinite,
                         int n_theta,
                         double dphi,
                         const std::string& meta_text)
{
    TLatex lat;
    lat.SetNDC(true);

    lat.SetTextAlign(13); // left-top
    lat.SetTextSize(0.040);
    lat.DrawLatex(0.05, 0.92, "p2MEG ACC grid PDF");

    lat.SetTextSize(0.030);
    lat.DrawLatex(0.05, 0.84, Form("input  : %s", infile));
    lat.DrawLatex(0.05, 0.79, Form("key    : %s", key));
    lat.DrawLatex(0.05, 0.74, Form("output : %s", outpdf));

    lat.SetTextSize(0.030);
    lat.DrawLatex(0.05, 0.66, Form("sum(prob mass) : %.6g", sum_mass));
    lat.DrawLatex(0.05, 0.61, Form("negative bins : %ld", n_negative));
    lat.DrawLatex(0.05, 0.56, Form("nonfinite bins: %ld", n_nonfinite));

    const TAxis* axE = h.GetAxis(0);
    const TAxis* axG = h.GetAxis(1);
    const TAxis* axPe = h.GetAxis(2);
    const TAxis* axPg = h.GetAxis(3);
    lat.SetTextSize(0.030);
    lat.DrawLatex(0.05, 0.48, Form("bins: Ee=%d Eg=%d phi_e=%d phi_g=%d",
                                   axE ? axE->GetNbins() : 0,
                                   axG ? axG->GetNbins() : 0,
                                   axPe ? axPe->GetNbins() : 0,
                                   axPg ? axPg->GetNbins() : 0));
    lat.DrawLatex(0.05, 0.43, Form("phi grid: N_theta=%d  dphi=%.6g rad", n_theta, dphi));

    lat.SetTextSize(0.030);
    lat.DrawLatex(0.05, 0.36, "Analysis window (Ee, Eg, t, theta):");
    lat.SetTextSize(0.028);
    lat.DrawLatex(0.08, 0.31, Form("Ee [MeV]    : [%.6g, %.6g]", analysis_window.Ee_min, analysis_window.Ee_max));
    lat.DrawLatex(0.08, 0.27, Form("Eg [MeV]    : [%.6g, %.6g]", analysis_window.Eg_min, analysis_window.Eg_max));
    lat.DrawLatex(0.08, 0.23, Form("t  [ns]     : [%.6g, %.6g]", analysis_window.t_min,  analysis_window.t_max));
    lat.DrawLatex(0.08, 0.19, Form("theta_eg [rad] : [%.6g, %.6g]", analysis_window.theta_min, analysis_window.theta_max));

    lat.SetTextSize(0.026);
    lat.DrawLatex(0.05, 0.13, "t: uniform in window (not plotted).  phi: discrete (N_theta+1 points).");

    lat.SetTextSize(0.024);
    lat.DrawLatex(0.05, 0.08, "Pages: (1) meta  (2) 1D  (3) 2D");

    if (!meta_text.empty()) {
        const std::vector<std::string> lines = SplitLines(meta_text);
        double y = 0.94;
        lat.SetTextSize(0.020);
        lat.DrawLatex(0.60, y, "meta (excerpt):");
        y -= 0.03;
        const size_t n_show = (lines.size() > 6) ? 6 : lines.size();
        for (size_t i = 0; i < n_show; ++i) {
            lat.DrawLatex(0.60, y, lines[i].c_str());
            y -= 0.03;
        }
    }
}

static void PrintUsage(const char* prog)
{
    std::cerr << "Usage:\n";
    std::cerr << "  " << prog << " [in.root] [key] [out.pdf]\n";
}

int main(int argc, char** argv)
{
    const char* infile = "data/pdf_cache/acc_grid.root";
    const char* key = "acc_grid";
    const char* outpdf = "doc/acc_grid_hist.pdf";

    if (argc >= 2) infile = argv[1];
    if (argc >= 3) key = argv[2];
    if (argc >= 4) outpdf = argv[3];
    if (argc > 4) {
        PrintUsage(argv[0]);
        return 1;
    }

    TFile f(infile, "READ");
    if (f.IsZombie()) {
        Error("plot_acc_grid_pdf", "cannot open: %s", infile);
        return 1;
    }

    TObject* obj = f.Get(key);
    if (!obj) {
        Error("plot_acc_grid_pdf", "key not found: %s", key);
        return 1;
    }

    THnD* grid_in = dynamic_cast<THnD*>(obj);
    if (!grid_in) {
        Error("plot_acc_grid_pdf", "object is not THnD: %s", key);
        return 1;
    }
    if (grid_in->GetNdimensions() != 4) {
        Error("plot_acc_grid_pdf", "unexpected ndim=%d (expected 4)", grid_in->GetNdimensions());
        return 1;
    }

    THnD* grid = dynamic_cast<THnD*>(grid_in->Clone("acc_grid_clone"));
    if (!grid) {
        Error("plot_acc_grid_pdf", "clone failed");
        return 1;
    }

    std::string meta_text;
    const std::string meta_name = std::string(key) + "_meta";
    TNamed* meta = nullptr;
    f.GetObject(meta_name.c_str(), meta);
    if (meta) meta_text = meta->GetTitle();

    int n_theta = 0;
    double dphi = 0.0;
    const std::string ntheta_name = std::string(key) + "_N_theta";
    const std::string dphi_name = std::string(key) + "_dphi";
    TParameter<int>* par_ntheta = nullptr;
    TParameter<double>* par_dphi = nullptr;
    f.GetObject(ntheta_name.c_str(), par_ntheta);
    f.GetObject(dphi_name.c_str(), par_dphi);
    if (par_ntheta) n_theta = par_ntheta->GetVal();
    if (par_dphi) dphi = par_dphi->GetVal();
    f.Close();

    gStyle->SetOptStat(0);
    gSystem->mkdir("doc", /*recursive=*/true);

    const double Ee_min = analysis_window.Ee_min;
    const double Ee_max = analysis_window.Ee_max;
    const double Eg_min = analysis_window.Eg_min;
    const double Eg_max = analysis_window.Eg_max;
    const double th_min = analysis_window.theta_min;
    const double th_max = analysis_window.theta_max;
    const double phi_min = 0.0;
    const double phi_max = pi;

    // ---- 1D ----
    TH1D* hEe   = new TH1D("hEe",   "Ee;Ee [MeV];Entries",                 kNBins_E,  Ee_min, Ee_max);
    TH1D* hEg   = new TH1D("hEg",   "Eg;Eg [MeV];Entries",                 kNBins_E,  Eg_min, Eg_max);
    TH1D* hPhiE = new TH1D("hPhiE", "phi_{detector,e};phi_{detector,e} [rad];Entries",
                           kNBins_phi, phi_min, phi_max);
    TH1D* hPhiG = new TH1D("hPhiG", "phi_{detector,#gamma};phi_{detector,#gamma} [rad];Entries",
                           kNBins_phi, phi_min, phi_max);
    TH1D* hThEg = new TH1D("hThEg", "theta_{eg};theta_{eg} [rad];Entries",
                           kNBins_th, th_min, th_max);

    // ---- 2D ----
    TH2D* h_EeEg = new TH2D("h_EeEg", "(Ee, Eg);Ee [MeV];Eg [MeV]",
                            kNBins2D_E, Ee_min, Ee_max, kNBins2D_E, Eg_min, Eg_max);

    TH2D* h_ThEe = new TH2D("h_ThEe", "(theta_{eg}, Ee);theta_{eg} [rad];Ee [MeV]",
                            kNBins2D_th, th_min, th_max, kNBins2D_E, Ee_min, Ee_max);

    TH2D* h_ThEg = new TH2D("h_ThEg", "(theta_{eg}, Eg);theta_{eg} [rad];Eg [MeV]",
                            kNBins2D_th, th_min, th_max, kNBins2D_E, Eg_min, Eg_max);

    TH2D* h_PePg = new TH2D("h_PePg", "(phi_{detector,e}, phi_{detector,#gamma});phi_{detector,e} [rad];phi_{detector,#gamma} [rad]",
                            kNBins2D_phi, phi_min, phi_max, kNBins2D_phi, phi_min, phi_max);

    const TAxis* axE = grid->GetAxis(0);
    const TAxis* axG = grid->GetAxis(1);
    const TAxis* axPe = grid->GetAxis(2);
    const TAxis* axPg = grid->GetAxis(3);

    if (!axE || !axG || !axPe || !axPg) {
        Error("plot_acc_grid_pdf", "missing axis in grid");
        return 1;
    }

    double sum_mass = 0.0;
    long n_negative = 0;
    long n_nonfinite = 0;

    const int n0 = axE->GetNbins();
    const int n1 = axG->GetNbins();
    const int n2 = axPe->GetNbins();
    const int n3 = axPg->GetNbins();

    std::vector<int> idx(4, 1);
    for (int i0 = 1; i0 <= n0; ++i0) {
        idx[0] = i0;
        const double wEe = axE->GetBinWidth(i0);
        const double Ee = axE->GetBinCenter(i0);
        for (int i1 = 1; i1 <= n1; ++i1) {
            idx[1] = i1;
            const double wEg = axG->GetBinWidth(i1);
            const double Eg = axG->GetBinCenter(i1);
            for (int i2 = 1; i2 <= n2; ++i2) {
                idx[2] = i2;
                const double phi_e = axPe->GetBinCenter(i2);
                for (int i3 = 1; i3 <= n3; ++i3) {
                    idx[3] = i3;
                    const double phi_g = axPg->GetBinCenter(i3);

                    const Long64_t bin = grid->GetBin(idx.data());
                    const double v = grid->GetBinContent(bin);
                    if (!IsFinite(v)) {
                        ++n_nonfinite;
                        continue;
                    }
                    if (v < 0.0) {
                        ++n_negative;
                        continue;
                    }
                    if (!(wEe > 0.0) || !(wEg > 0.0)) continue;

                    const double mass = v * wEe * wEg;
                    if (!(mass > 0.0) || !IsFinite(mass)) continue;

                    sum_mass += mass;

                    const double theta_eg = std::fabs(phi_e - phi_g);
                    const bool in_theta = (theta_eg >= th_min && theta_eg <= th_max);

                    // 1D
                    hEe->Fill(Ee, mass);
                    hEg->Fill(Eg, mass);
                    hPhiE->Fill(phi_e, mass);
                    hPhiG->Fill(phi_g, mass);
                    if (in_theta) hThEg->Fill(theta_eg, mass);

                    // 2D
                    h_EeEg->Fill(Ee, Eg, mass);
                    if (in_theta) h_ThEe->Fill(theta_eg, Ee, mass);
                    if (in_theta) h_ThEg->Fill(theta_eg, Eg, mass);
                    h_PePg->Fill(phi_e, phi_g, mass);
                }
            }
        }
    }

    if (!(sum_mass > 0.0) || !IsFinite(sum_mass)) {
        Error("plot_acc_grid_pdf", "sum(prob mass) is not positive: %.6g", sum_mass);
        return 2;
    }

    // ---- ページ1：メタ ----
    TCanvas c0("c0", "meta", 1200, 800);
    c0.cd();
    DrawMetaPage(infile, key, outpdf, *grid, sum_mass,
                 n_negative, n_nonfinite, n_theta, dphi, meta_text);

    // ---- ページ2：1D ----
    TCanvas c1("c1", "1D acc grid", 1200, 800);
    c1.Divide(3, 2);

    c1.cd(1); gPad->SetGrid(); hEe->SetLineWidth(2); hEe->Draw("hist");
    c1.cd(2); gPad->SetGrid(); hEg->SetLineWidth(2); hEg->Draw("hist");
    c1.cd(3); gPad->SetGrid(); hPhiE->SetLineWidth(2); hPhiE->Draw("hist");
    c1.cd(4); gPad->SetGrid(); hPhiG->SetLineWidth(2); hPhiG->Draw("hist");
    c1.cd(5); gPad->SetGrid(); hThEg->SetLineWidth(2); hThEg->Draw("hist");

    // ---- ページ3：2D ----
    auto Draw2D = [](TH2D* h){
        gPad->SetGrid();
        gPad->SetRightMargin(0.14);
        h->Draw("colz");
    };

    TCanvas c2("c2", "2D acc grid", 1200, 800);
    c2.Divide(2, 2);
    c2.cd(1); Draw2D(h_EeEg);
    c2.cd(2); Draw2D(h_ThEe);
    c2.cd(3); Draw2D(h_ThEg);
    c2.cd(4); Draw2D(h_PePg);

    // ---- PDF（3ページ）----
    TString out = outpdf;
    c0.Print(Form("%s[", out.Data()));
    c0.Print(out.Data());
    c1.Print(out.Data());
    c2.Print(out.Data());
    c2.Print(Form("%s]", out.Data()));

    Info("plot_acc_grid_pdf", "wrote: %s (3 pages)", out.Data());
    Info("plot_acc_grid_pdf", "sum_mass=%.6g, negative=%ld, nonfinite=%ld",
         sum_mass, n_negative, n_nonfinite);

    return 0;
}
