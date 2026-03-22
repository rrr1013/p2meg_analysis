// scripts/plot_pdf_shapes_final.cc
//
// 論文用 fig5:
// Signal / RMD / ACC の 1D PDF 形状を 2x2 で重ね描きする。
// 連続変数は細かいビンで積分し、TGraph の折れ線で滑らかに表示する。

#include <cmath>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPad.h"
#include "TH1D.h"
#include "TSpline.h"
#include "TStyle.h"
#include "TSystem.h"

#include "p2meg/ACCGridPdf.h"
#include "p2meg/AnalysisWindow.h"
#include "p2meg/Constants.h"
#include "p2meg/DetectorResolution.h"
#include "p2meg/MathUtils.h"
#include "p2meg/RMDGridPdf.h"
#include "p2meg/SignalPdf.h"

static constexpr int kNBinsE = 80;
static constexpr int kNBinsT = 160;
static constexpr int kNBinsTheta = 36;

static constexpr int kColorSignal = kRed + 1;
static constexpr int kColorRmd = kBlue + 1;
static constexpr int kColorAcc = kGreen + 2;

static const char* kRmdRoot = "data/pdf_cache/rmd_grid.root";
static const char* kRmdKey = "rmd_grid";
static const char* kAccRoot = "data/pdf_cache/acc_grid.root";
static const char* kAccKey = "acc_grid";

static void ConvertHistToDensity(TH1D* h)
{
    if (!h) return;
    for (int ib = 1; ib <= h->GetNbinsX(); ++ib) {
        const double w = h->GetXaxis()->GetBinWidth(ib);
        if (!(w > 0.0)) continue;
        h->SetBinContent(ib, h->GetBinContent(ib) / w);
    }
}

static TGraph* MakeGraphFromHist(const TH1D* h, int color)
{
    if (!h) return nullptr;
    const int n = h->GetNbinsX();
    TGraph* gr = new TGraph(n);
    for (int i = 0; i < n; ++i) {
        const int ib = i + 1;
    gr->SetPoint(i, h->GetXaxis()->GetBinCenter(ib), h->GetBinContent(ib));
    }
    gr->SetLineColor(color);
    gr->SetLineWidth(1);
    return gr;
}

static TSpline3* MakeSplineFromHist(const TH1D* h, int color, const char* name)
{
    TGraph* gr = MakeGraphFromHist(h, color);
    if (!gr) return nullptr;
    TSpline3* sp = new TSpline3(name, gr);
    sp->SetLineColor(color);
    sp->SetLineWidth(1);
    return sp;
}

static void StyleAxis(TH1D& h)
{
    h.GetXaxis()->SetTitleSize(0.055);
    h.GetYaxis()->SetTitleSize(0.055);
    h.GetXaxis()->SetLabelSize(0.043);
    h.GetYaxis()->SetLabelSize(0.043);
    h.GetXaxis()->SetTitleOffset(1.20);
    h.GetYaxis()->SetTitleOffset(1.30);
}

int main()
{
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    if (!RMDGridPdf_Load(kRmdRoot, kRmdKey)) {
        std::cerr << "[plot_pdf_shapes_final] RMDGridPdf_Load failed\n";
        return 1;
    }
    if (!ACCGridPdf_Load(kAccRoot, kAccKey)) {
        std::cerr << "[plot_pdf_shapes_final] ACCGridPdf_Load failed\n";
        return 1;
    }

    TH1D hSigEe("hSigEe", ";E_{e^{+}} [MeV];Density [1/MeV]",
                kNBinsE, analysis_window.Ee_min, analysis_window.Ee_max);
    TH1D hRmdEe("hRmdEe", ";E_{e^{+}} [MeV];Density [1/MeV]",
                kNBinsE, analysis_window.Ee_min, analysis_window.Ee_max);
    TH1D hAccEe("hAccEe", ";E_{e^{+}} [MeV];Density [1/MeV]",
                kNBinsE, analysis_window.Ee_min, analysis_window.Ee_max);

    TH1D hSigEg("hSigEg", ";E_{#gamma} [MeV];Density [1/MeV]",
                kNBinsE, analysis_window.Eg_min, analysis_window.Eg_max);
    TH1D hRmdEg("hRmdEg", ";E_{#gamma} [MeV];Density [1/MeV]",
                kNBinsE, analysis_window.Eg_min, analysis_window.Eg_max);
    TH1D hAccEg("hAccEg", ";E_{#gamma} [MeV];Density [1/MeV]",
                kNBinsE, analysis_window.Eg_min, analysis_window.Eg_max);

    TH1D hSigT("hSigT", ";t_{e^{+}#gamma} [ns];Density [1/ns]",
               kNBinsT, analysis_window.t_min, Math_AxisMaxInclusive(analysis_window.t_max));
    TH1D hRmdT("hRmdT", ";t_{e^{+}#gamma} [ns];Density [1/ns]",
               kNBinsT, analysis_window.t_min, Math_AxisMaxInclusive(analysis_window.t_max));
    TH1D hAccT("hAccT", ";t_{e^{+}#gamma} [ns];Density [1/ns]",
               kNBinsT, analysis_window.t_min, Math_AxisMaxInclusive(analysis_window.t_max));

    TH1D hSigTheta("hSigTheta", ";#theta_{e^{+}#gamma} [rad];Density [1/rad]",
                   kNBinsTheta, analysis_window.theta_min, Math_AxisMaxInclusive(analysis_window.theta_max));
    TH1D hRmdTheta("hRmdTheta", ";#theta_{e^{+}#gamma} [rad];Density [1/rad]",
                   kNBinsTheta, analysis_window.theta_min, Math_AxisMaxInclusive(analysis_window.theta_max));
    TH1D hAccTheta("hAccTheta", ";#theta_{e^{+}#gamma} [rad];Density [1/rad]",
                   kNBinsTheta, analysis_window.theta_min, Math_AxisMaxInclusive(analysis_window.theta_max));

    StyleAxis(hSigEe); StyleAxis(hRmdEe); StyleAxis(hAccEe);
    StyleAxis(hSigEg); StyleAxis(hRmdEg); StyleAxis(hAccEg);
    StyleAxis(hSigT); StyleAxis(hRmdT); StyleAxis(hAccT);
    StyleAxis(hSigTheta); StyleAxis(hRmdTheta); StyleAxis(hAccTheta);

    const int N_phi_e = Math_GetNPhiE(detres);
    const int N_phi_g = Math_GetNPhiG(detres);

    for (int ibe = 1; ibe <= kNBinsE; ++ibe) {
        const double Ee = hSigEe.GetXaxis()->GetBinCenter(ibe);
        const double dEe = hSigEe.GetXaxis()->GetBinWidth(ibe);
        for (int ibg = 1; ibg <= kNBinsE; ++ibg) {
            const double Eg = hSigEg.GetXaxis()->GetBinCenter(ibg);
            const double dEg = hSigEg.GetXaxis()->GetBinWidth(ibg);
            for (int ibt = 1; ibt <= kNBinsT; ++ibt) {
                const double t = hSigT.GetXaxis()->GetBinCenter(ibt);
                const double dt = hSigT.GetXaxis()->GetBinWidth(ibt);
                for (int ie = 0; ie <= N_phi_e; ++ie) {
                    const double phi_e = Detector_PhiGridPoint(ie, detres.phi_e_min, detres.phi_e_max, N_phi_e);
                    const double dphi_e = Detector_PhiBinWidth(ie, detres.phi_e_min, detres.phi_e_max, N_phi_e);
                    if (!(dphi_e > 0.0)) continue;
                    for (int ig = 0; ig <= N_phi_g; ++ig) {
                        if (!Detector_IsAllowedPhiPairIndex(ie, ig, detres)) continue;
                        const double phi_g = Detector_PhiGridPoint(ig, detres.phi_g_min, detres.phi_g_max, N_phi_g);
                        const double dphi_g = Detector_PhiBinWidth(ig, detres.phi_g_min, detres.phi_g_max, N_phi_g);
                        if (!(dphi_g > 0.0)) continue;
                        const double dV = dEe * dEg * dt * dphi_e * dphi_g;
                        const double theta = std::fabs(phi_e - phi_g);

                        const double ps = SignalPdf(Ee, Eg, t, phi_e, phi_g,
                                                    analysis_window, detres, kMassesPDG);
                        const double pr = RMDGridPdf(Ee, Eg, t, phi_e, phi_g);
                        const double pa = ACCGridPdf(Ee, Eg, t, phi_e, phi_g);

                        if (ps > 0.0 && std::isfinite(ps)) {
                            const double m = ps * dV;
                            hSigEe.Fill(Ee, m); hSigEg.Fill(Eg, m); hSigT.Fill(t, m); hSigTheta.Fill(theta, m);
                        }
                        if (pr > 0.0 && std::isfinite(pr)) {
                            const double m = pr * dV;
                            hRmdEe.Fill(Ee, m); hRmdEg.Fill(Eg, m); hRmdT.Fill(t, m); hRmdTheta.Fill(theta, m);
                        }
                        if (pa > 0.0 && std::isfinite(pa)) {
                            const double m = pa * dV;
                            hAccEe.Fill(Ee, m); hAccEg.Fill(Eg, m); hAccT.Fill(t, m); hAccTheta.Fill(theta, m);
                        }
                    }
                }
            }
        }
    }

    ConvertHistToDensity(&hSigEe); ConvertHistToDensity(&hRmdEe); ConvertHistToDensity(&hAccEe);
    ConvertHistToDensity(&hSigEg); ConvertHistToDensity(&hRmdEg); ConvertHistToDensity(&hAccEg);
    ConvertHistToDensity(&hSigT);  ConvertHistToDensity(&hRmdT);  ConvertHistToDensity(&hAccT);
    ConvertHistToDensity(&hSigTheta); ConvertHistToDensity(&hRmdTheta); ConvertHistToDensity(&hAccTheta);

    TSpline3* sSigEe = MakeSplineFromHist(&hSigEe, kColorSignal, "sSigEe");
    TSpline3* sRmdEe = MakeSplineFromHist(&hRmdEe, kColorRmd, "sRmdEe");
    TSpline3* sAccEe = MakeSplineFromHist(&hAccEe, kColorAcc, "sAccEe");
    TSpline3* sSigEg = MakeSplineFromHist(&hSigEg, kColorSignal, "sSigEg");
    TSpline3* sRmdEg = MakeSplineFromHist(&hRmdEg, kColorRmd, "sRmdEg");
    TSpline3* sAccEg = MakeSplineFromHist(&hAccEg, kColorAcc, "sAccEg");
    TSpline3* sSigT = MakeSplineFromHist(&hSigT, kColorSignal, "sSigT");
    TSpline3* sRmdT = MakeSplineFromHist(&hRmdT, kColorRmd, "sRmdT");
    TSpline3* sAccT = MakeSplineFromHist(&hAccT, kColorAcc, "sAccT");
    TSpline3* sSigTheta = MakeSplineFromHist(&hSigTheta, kColorSignal, "sSigTheta");
    TSpline3* sRmdTheta = MakeSplineFromHist(&hRmdTheta, kColorRmd, "sRmdTheta");
    TSpline3* sAccTheta = MakeSplineFromHist(&hAccTheta, kColorAcc, "sAccTheta");

    gSystem->mkdir("doc/fig", true);

    TCanvas c("c_pdf_shapes_final", "PDF shapes final", 1180, 930);

    TPad pLeg("pLeg", "legend", 0.00, 0.90, 1.00, 1.00);
    TPad p1("p1", "p1", 0.06, 0.51, 0.49, 0.90);
    TPad p2("p2", "p2", 0.52, 0.51, 0.95, 0.90);
    TPad p3("p3", "p3", 0.06, 0.09, 0.49, 0.48);
    TPad p4("p4", "p4", 0.52, 0.09, 0.95, 0.48);
    pLeg.Draw(); p1.Draw(); p2.Draw(); p3.Draw(); p4.Draw();

    p1.cd();
    gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.18); gPad->SetRightMargin(0.08); gPad->SetTopMargin(0.06);
    double ymax = std::max(hSigEe.GetMaximum(), std::max(hRmdEe.GetMaximum(), hAccEe.GetMaximum()));
    hSigEe.SetMaximum(1.20 * ymax);
    hSigEe.Draw("axis");
    sSigEe->Draw("same");
    sRmdEe->Draw("same");
    sAccEe->Draw("same");

    p2.cd();
    gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.18); gPad->SetRightMargin(0.08); gPad->SetTopMargin(0.06);
    ymax = std::max(hSigEg.GetMaximum(), std::max(hRmdEg.GetMaximum(), hAccEg.GetMaximum()));
    hSigEg.SetMaximum(1.20 * ymax);
    hSigEg.Draw("axis");
    sSigEg->Draw("same");
    sRmdEg->Draw("same");
    sAccEg->Draw("same");

    p3.cd();
    gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.18); gPad->SetRightMargin(0.08); gPad->SetTopMargin(0.06);
    ymax = std::max(hSigT.GetMaximum(), std::max(hRmdT.GetMaximum(), hAccT.GetMaximum()));
    hSigT.SetMaximum(1.20 * ymax);
    hSigT.Draw("axis");
    sRmdT->Draw("same");
    sSigT->SetLineStyle(2);
    sSigT->SetLineWidth(2);
    sSigT->Draw("same");
    sAccT->Draw("same");

    p4.cd();
    gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.18); gPad->SetRightMargin(0.08); gPad->SetTopMargin(0.06);
    ymax = std::max(hSigTheta.GetMaximum(), std::max(hRmdTheta.GetMaximum(), hAccTheta.GetMaximum()));
    hSigTheta.SetMaximum(1.20 * ymax);
    hSigTheta.Draw("axis");
    sSigTheta->Draw("same");
    sRmdTheta->Draw("same");
    sAccTheta->Draw("same");

    pLeg.cd();
    gPad->SetFillStyle(0);
    TLegend leg(0.38, 0.18, 0.88, 0.82);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetNColumns(3);
    leg.SetTextSize(0.34);
    leg.AddEntry(sSigEe, "Signal PDF", "l");
    leg.AddEntry(sRmdEe, "RMD PDF", "l");
    leg.AddEntry(sAccEe, "ACC PDF", "l");
    leg.Draw();

    TLatex lat;
    lat.SetNDC(true);
    lat.SetTextFont(42);
    lat.SetTextSize(0.34);
    lat.DrawLatex(0.08, 0.33, "PDF shapes");

    c.SaveAs("doc/fig/fig05_signal_rmd_acc_pdf_shapes.pdf");
}
