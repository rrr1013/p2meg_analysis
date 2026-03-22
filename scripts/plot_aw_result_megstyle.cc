// scripts/plot_aw_result_megstyle.cc
//
// 論文用 fig7:
// MEG Fig.27 の雰囲気に合わせ、黒点の散布図に signal contour を重ねる。

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "TCanvas.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TPad.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"

#include "p2meg/AnalysisWindow.h"
#include "p2meg/DetectorResolution.h"
#include "p2meg/Constants.h"

struct AwEventPoint {
    double Ee;
    double Eg;
    double t;
    double theta;
};

static bool ParseAwLine(const std::string& line,
                        double& Ee, double& Eg, double& t,
                        double& phi_e, double& phi_g)
{
    std::istringstream iss(line);
    if (!(iss >> Ee >> Eg >> t >> phi_e >> phi_g)) return false;
    return std::isfinite(Ee) && std::isfinite(Eg) &&
           std::isfinite(t) && std::isfinite(phi_e) && std::isfinite(phi_g);
}

int main()
{
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    std::ifstream fin("data/finaldata/step4f_run1to20_eg_doubleonly_allpatterns_500ns.txt");
    if (!fin) {
        std::cerr << "[plot_aw_result_megstyle] cannot open final data\n";
        return 1;
    }

    std::vector<AwEventPoint> evs;
    std::string line;
    while (std::getline(fin, line)) {
        double Ee = 0.0, Eg = 0.0, t = 0.0, phi_e = 0.0, phi_g = 0.0;
        if (!ParseAwLine(line, Ee, Eg, t, phi_e, phi_g)) continue;
        int idx_e = -1, idx_g = -1;
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
        evs.push_back({Ee, Eg, t, theta});
    }

    TGraph grEeEg(static_cast<int>(evs.size()));
    TGraph grThetaT(static_cast<int>(evs.size()));
    for (int i = 0; i < static_cast<int>(evs.size()); ++i) {
        grEeEg.SetPoint(i, evs[i].Ee, evs[i].Eg);
        grThetaT.SetPoint(i, evs[i].theta, evs[i].t);
    }
    grEeEg.SetMarkerStyle(20);
    grEeEg.SetMarkerSize(1.1);
    grEeEg.SetMarkerColor(kBlack);
    grThetaT.SetMarkerStyle(20);
    grThetaT.SetMarkerSize(1.1);
    grThetaT.SetMarkerColor(kBlack);

    gSystem->mkdir("doc/fig", true);

    TCanvas c("c_aw_megstyle", "AW result MEG style", 1120, 600);
    TPad p1("p1", "p1", 0.06, 0.12, 0.47, 0.95);
    TPad p2("p2", "p2", 0.54, 0.12, 0.95, 0.95);
    p1.Draw();
    p2.Draw();

    const double E0 = 0.5 * kMassesPDG.m_mu;
    const double sigmaE = 0.079 * E0;
    const double sigmaTheta = pi / (2.0 * detres.N_theta);
    const double sigmaT = detres.sigma_t;
    const double thetaPlotMax = analysis_window.theta_max + 3.4 * sigmaTheta;

    auto draw_gauss_contour = [](TH2D& h, double x0, double sx, double y0, double sy,
                                 double nsigma, int line_style) {
        const int nx = h.GetNbinsX();
        const int ny = h.GetNbinsY();
        const double target = std::exp(-0.5 * nsigma * nsigma);
        for (int ix = 1; ix <= nx; ++ix) {
            const double x = h.GetXaxis()->GetBinCenter(ix);
            for (int iy = 1; iy <= ny; ++iy) {
                const double y = h.GetYaxis()->GetBinCenter(iy);
                const double dx = (x - x0) / sx;
                const double dy = (y - y0) / sy;
                h.SetBinContent(ix, iy, std::exp(-0.5 * (dx * dx + dy * dy)));
            }
        }
        const double level = target;
        h.SetContour(1, &level);
        h.SetLineColor(kRed + 1);
        h.SetLineStyle(line_style);
        h.SetLineWidth(1);
        h.Draw("CONT3 SAME");
    };

    p1.cd();
    gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.14); gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.04);
    gPad->SetFixedAspectRatio();
    TH2D frameEeEg("frameEeEg", ";E_{e^{+}} [MeV];E_{#gamma} [MeV]",
                   100, analysis_window.Ee_min, analysis_window.Ee_max,
                   100, analysis_window.Eg_min, analysis_window.Eg_max);
    frameEeEg.Draw();
    grEeEg.Draw("P same");
    TH2D hSigEeEg1("hSigEeEg1", "", 180, analysis_window.Ee_min, analysis_window.Ee_max,
                   180, analysis_window.Eg_min, analysis_window.Eg_max);
    TH2D hSigEeEg2("hSigEeEg2", "", 180, analysis_window.Ee_min, analysis_window.Ee_max,
                   180, analysis_window.Eg_min, analysis_window.Eg_max);
    TH2D hSigEeEg3("hSigEeEg3", "", 180, analysis_window.Ee_min, analysis_window.Ee_max,
                   180, analysis_window.Eg_min, analysis_window.Eg_max);
    draw_gauss_contour(hSigEeEg3, E0, sigmaE, E0, sigmaE, 3.0, 2);
    draw_gauss_contour(hSigEeEg2, E0, sigmaE, E0, sigmaE, 2.0, 1);
    draw_gauss_contour(hSigEeEg1, E0, sigmaE, E0, sigmaE, 1.0, 2);

    TLatex lat;
    lat.SetTextFont(42);
    lat.SetTextSize(0.040);
    lat.DrawLatex(E0 + 0.15 * sigmaE, E0 + 0.25 * sigmaE, "1#sigma");
    lat.DrawLatex(E0 + 0.75 * sigmaE, E0 + 0.95 * sigmaE, "2#sigma");
    lat.DrawLatex(E0 + 1.40 * sigmaE, E0 + 1.65 * sigmaE, "3#sigma");

    p2.cd();
    gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.14); gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.04);
    gPad->SetFixedAspectRatio();
    TH2D frameThetaT("frameThetaT", ";#theta_{e^{+}#gamma} [rad];t_{e^{+}#gamma} [ns]",
                     100, analysis_window.theta_min, thetaPlotMax,
                     100, analysis_window.t_min, analysis_window.t_max);
    frameThetaT.Draw();
    grThetaT.Draw("P same");
    const double theta0 = analysis_window.theta_max - 0.02;
    TH2D hSigThetaT1("hSigThetaT1", "", 180, analysis_window.theta_min, thetaPlotMax,
                     180, analysis_window.t_min, analysis_window.t_max);
    TH2D hSigThetaT2("hSigThetaT2", "", 180, analysis_window.theta_min, thetaPlotMax,
                     180, analysis_window.t_min, analysis_window.t_max);
    TH2D hSigThetaT3("hSigThetaT3", "", 180, analysis_window.theta_min, thetaPlotMax,
                     180, analysis_window.t_min, analysis_window.t_max);
    draw_gauss_contour(hSigThetaT3, theta0, sigmaTheta, detres.t_mean, sigmaT, 3.0, 2);
    draw_gauss_contour(hSigThetaT2, theta0, sigmaTheta, detres.t_mean, sigmaT, 2.0, 1);
    draw_gauss_contour(hSigThetaT1, theta0, sigmaTheta, detres.t_mean, sigmaT, 1.0, 2);

    lat.DrawLatex(theta0 + 0.10 * sigmaTheta, detres.t_mean - 0.20 * sigmaT, "1#sigma");
    lat.DrawLatex(theta0 + 0.60 * sigmaTheta, detres.t_mean + 0.85 * sigmaT, "2#sigma");
    lat.DrawLatex(theta0 + 1.10 * sigmaTheta, detres.t_mean + 1.75 * sigmaT, "3#sigma");

    c.SaveAs("doc/fig/fig07_aw_event_distribution.pdf");
    c.SaveAs("doc/fig/fig07_aw_event_distribution.png");
    return 0;
}
