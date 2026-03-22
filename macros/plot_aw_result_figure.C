// macros/plot_aw_result_figure.C
//
// 目的:
//   結果節に入れるための AW 内イベント分布 2 panel 図を 1 ページで作る。
//   左に (Ee, Eg)、右に (theta_eg, t) を描く。
//
// 出力:
//   doc/fig/fig07_aw_event_distribution.pdf
//   doc/fig/fig07_aw_event_distribution.png

R__ADD_INCLUDE_PATH(./include)

#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

#include "TCanvas.h"
#include "TH2D.h"
#include "TLine.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"

#include "p2meg/AnalysisWindow.h"
#include "p2meg/DetectorResolution.h"

static bool ParseEventLineAw(const std::string& line,
                             double& Ee, double& Eg, double& t,
                             double& phi_e, double& phi_g)
{
    std::istringstream iss(line);
    if (!(iss >> Ee >> Eg >> t >> phi_e >> phi_g)) return false;
    return std::isfinite(Ee) && std::isfinite(Eg) &&
           std::isfinite(t) && std::isfinite(phi_e) && std::isfinite(phi_g);
}

void plot_aw_result_figure(
    const char* infile =
        "data/finaldata/step4f_run1to20_eg_doubleonly_allpatterns_500ns.txt")
{
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    std::ifstream fin(infile);
    if (!fin) {
        ::Error("plot_aw_result_figure", "cannot open %s", infile);
        return;
    }

    TH2D hEeEg("hEeEg",
               ";E_{e^{+}} [MeV];E_{#gamma} [MeV]",
               60, analysis_window.Ee_min, analysis_window.Ee_max,
               60, analysis_window.Eg_min, analysis_window.Eg_max);
    TH2D hThetaT("hThetaT",
                 ";#theta_{e^{+}#gamma} [rad];t_{e^{+}#gamma} [ns]",
                 36, analysis_window.theta_min, analysis_window.theta_max,
                 50, analysis_window.t_min, analysis_window.t_max);

    std::string line;
    while (std::getline(fin, line)) {
        double Ee = 0.0;
        double Eg = 0.0;
        double t = 0.0;
        double phi_e = 0.0;
        double phi_g = 0.0;
        if (!ParseEventLineAw(line, Ee, Eg, t, phi_e, phi_g)) continue;

        int idx_e = -1;
        int idx_g = -1;
        if (!Detector_IsAllowedPhiPairValue(phi_e, phi_g, detres, idx_e, idx_g)) continue;

        const double phi_e_disc =
            Detector_PhiGridPoint(idx_e, detres.phi_e_min, detres.phi_e_max, detres.N_phi_e);
        const double phi_g_disc =
            Detector_PhiGridPoint(idx_g, detres.phi_g_min, detres.phi_g_max, detres.N_phi_g);
        const double theta = std::fabs(phi_e_disc - phi_g_disc);

        if (Ee < analysis_window.Ee_min || Ee > analysis_window.Ee_max) continue;
        if (Eg < analysis_window.Eg_min || Eg > analysis_window.Eg_max) continue;
        if (t < analysis_window.t_min || t > analysis_window.t_max) continue;
        if (theta < analysis_window.theta_min || theta > analysis_window.theta_max) continue;

        hEeEg.Fill(Ee, Eg);
        hThetaT.Fill(theta, t);
    }

    gSystem->mkdir("doc/fig", /*recursive=*/true);

    TCanvas c("c_aw_result", "AW result", 1180, 560);
    c.Divide(2, 1, 0.025, 0.0);

    c.cd(1);
    gPad->SetLeftMargin(0.12);
    gPad->SetBottomMargin(0.12);
    gPad->SetRightMargin(0.12);
    gPad->SetTopMargin(0.10);
    hEeEg.SetContour(80);
    hEeEg.Draw("colz");
    TLine lee(52.8, analysis_window.Eg_min, 52.8, analysis_window.Eg_max);
    TLine leg(analysis_window.Ee_min, 52.8, analysis_window.Ee_max, 52.8);
    lee.SetLineColor(kRed + 1);
    leg.SetLineColor(kRed + 1);
    lee.SetLineStyle(2);
    leg.SetLineStyle(2);
    lee.SetLineWidth(2);
    leg.SetLineWidth(2);
    lee.Draw("same");
    leg.Draw("same");

    TLatex lat;
    lat.SetNDC(true);
    lat.SetTextFont(42);
    lat.SetTextSize(0.040);
    lat.DrawLatex(0.14, 0.92, "Events in the analysis window");
    lat.SetTextSize(0.032);
    lat.DrawLatex(0.50, 0.84, "dashed lines: nominal signal point");

    c.cd(2);
    gPad->SetLeftMargin(0.12);
    gPad->SetBottomMargin(0.12);
    gPad->SetRightMargin(0.12);
    gPad->SetTopMargin(0.10);
    hThetaT.SetContour(80);
    hThetaT.Draw("colz");
    TLine ltheta(3.1415926536, analysis_window.t_min, 3.1415926536, analysis_window.t_max);
    TLine lt(analysis_window.theta_min, 0.0, analysis_window.theta_max, 0.0);
    ltheta.SetLineColor(kRed + 1);
    lt.SetLineColor(kRed + 1);
    ltheta.SetLineStyle(2);
    lt.SetLineStyle(2);
    ltheta.SetLineWidth(2);
    lt.SetLineWidth(2);
    ltheta.Draw("same");
    lt.Draw("same");

    c.SaveAs("doc/fig/fig07_aw_event_distribution.pdf");
    c.SaveAs("doc/fig/fig07_aw_event_distribution.png");
}
