// macros/plot_data_hist_thesis.C
//
// 目的:
//   論文・発表用の全データ 2D ヒストグラムを作る。
//   - 論文用: 6枚組 1ページ PDF
//   - 発表用: (Ee,Eg), (t,Eg) の 2枚組 1ページ PDF
//
// 入力:
//   5列テキスト
//     Ee  Eg  t  phi_detector_e  phi_detector_g
//   単位:
//     Ee, Eg [MeV], t [ns], phi [rad]
//
// 表示方針:
//   - 解析用のログページは出さない
//   - 2D 図のみを整理して配置する
//   - notation は論文向けに
//       E_{e^+}, E_{#gamma}, #phi_{e^+}, #phi_{#gamma},
//       #theta_{e^+#gamma}, t_{e^+#gamma}
//   - (t, E_{#gamma}) には AW / TSB の領域を重ね描きする
//
// 実行例:
//   root -l -b -q 'macros/plot_data_hist_thesis.C("data/finaldata/step4f_run1to20_eg_doubleonly_allpatterns_500ns.txt")'

R__ADD_INCLUDE_PATH(./include)

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "TBox.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TError.h"
#include "TGaxis.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"

#include "p2meg/AnalysisWindow.h"
#include "p2meg/DetectorResolution.h"
#include "p2meg/Event.h"
#include "p2meg/MathUtils.h"

// 全データ図として使う表示範囲
static constexpr double kEePlotMin = 0.0;    // [MeV]
static constexpr double kEePlotMax = 100.0;  // [MeV]
static constexpr double kEgPlotMin = 0.0;    // [MeV]
static constexpr double kEgPlotMax = 100.0;  // [MeV]
static constexpr double kTPlotMin  = -500.0; // [ns]
static constexpr double kTPlotMax  = 500.0;  // [ns]

static bool ParseEventLine5Doubles(const std::string& line,
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

static TString MakeBaseName(const char* infile)
{
    TString base = gSystem->BaseName(infile);
    const Ssiz_t dot = base.Last('.');
    if (dot != kNPOS) base.Remove(dot);
    return base;
}

static std::vector<double> BuildThetaEdgesFromAllowedPhi(const DetectorResolutionConst& res)
{
    // theta_{e^+gamma}=|phi_e-phi_g| は離散値しか取らないので、
    // 最近傍の中点で区切った可変ビンにして「線」ではなく「帯」に見せる。
    std::set<double> theta_values;
    for (int ie = 0; ie <= res.N_phi_e; ++ie) {
        for (int ig = 0; ig <= res.N_phi_g; ++ig) {
            if (!Detector_IsAllowedPhiPairIndex(ie, ig, res)) continue;
            const double phi_e = Detector_PhiGridPoint(ie, res.phi_e_min, res.phi_e_max, res.N_phi_e);
            const double phi_g = Detector_PhiGridPoint(ig, res.phi_g_min, res.phi_g_max, res.N_phi_g);
            const double theta = std::fabs(phi_e - phi_g);
            if (std::isfinite(theta)) theta_values.insert(theta);
        }
    }

    std::vector<double> vals(theta_values.begin(), theta_values.end());
    std::vector<double> edges;
    if (vals.empty()) {
        edges.push_back(analysis_window.theta_min);
        edges.push_back(Math_AxisMaxInclusive(analysis_window.theta_max));
        return edges;
    }

    edges.resize(vals.size() + 1);
    if (vals.size() == 1) {
        const double half_width = 0.15;
        edges[0] = std::max(0.0, vals[0] - half_width);
        edges[1] = Math_AxisMaxInclusive(std::min(pi, vals[0] + half_width));
        return edges;
    }

    edges[0] = std::max(0.0, vals[0] - 0.5 * (vals[1] - vals[0]));
    for (std::size_t i = 1; i < vals.size(); ++i) {
        edges[i] = 0.5 * (vals[i - 1] + vals[i]);
    }
    edges[vals.size()] = Math_AxisMaxInclusive(pi);

    for (std::size_t i = 1; i < edges.size(); ++i) {
        if (!(edges[i] > edges[i - 1])) edges[i] = Math_AxisMaxInclusive(edges[i - 1]);
    }
    return edges;
}

static void ApplyGlobalStyle()
{
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetTitleFillColor(0);
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTitleSize(0.055, "XYZ");
    gStyle->SetLabelSize(0.045, "XYZ");
    gStyle->SetTitleOffset(1.05, "X");
    gStyle->SetTitleOffset(1.15, "Y");
    gStyle->SetPalette(kBird);
    TGaxis::SetMaxDigits(3);
}

static void StyleHistogram2D(TH2D* h)
{
    if (!h) return;
    h->SetContour(80);
    h->GetXaxis()->CenterTitle(false);
    h->GetYaxis()->CenterTitle(false);
    h->GetXaxis()->SetNdivisions(505);
    h->GetYaxis()->SetNdivisions(505);
    h->GetXaxis()->SetLabelSize(0.050);
    h->GetYaxis()->SetLabelSize(0.050);
    h->GetXaxis()->SetTitleSize(0.058);
    h->GetYaxis()->SetTitleSize(0.058);
    h->GetZaxis()->SetLabelSize(0.045);
    h->GetZaxis()->SetTitleSize(0.052);
    h->GetZaxis()->SetTitle("Counts");
}

static void SetPadMarginsForPaper()
{
    gPad->SetLeftMargin(0.14);
    gPad->SetBottomMargin(0.14);
    gPad->SetRightMargin(0.16);
    gPad->SetTopMargin(0.10);
}

static void SetPadMarginsForSlides()
{
    gPad->SetLeftMargin(0.18);
    gPad->SetBottomMargin(0.13);
    gPad->SetRightMargin(0.18);
    gPad->SetTopMargin(0.12);
}

static void DrawRegionBox(double x1, double y1, double x2, double y2,
                          int line_color, int fill_color, int line_style,
                          double alpha)
{
    TBox* box = new TBox(x1, y1, x2, y2);
    box->SetFillColorAlpha(fill_color, alpha);
    box->SetLineColor(line_color);
    box->SetLineStyle(line_style);
    box->SetLineWidth(2);
    box->Draw("lfsame");
}

static void DrawLineSegment(double x1, double y1, double x2, double y2,
                            int line_color, int line_style, int line_width)
{
    TLine* line = new TLine(x1, y1, x2, y2);
    line->SetLineColor(line_color);
    line->SetLineStyle(line_style);
    line->SetLineWidth(line_width);
    line->Draw("same");
}

static void DrawAwTsbOverlayOnTEg()
{
    // Fig.22 風に、AW を中央、TSB を左右に明示する。
    // p2MEG ではブラインドしないので、塗りは薄く残す。
    const double eg1 = analysis_window.Eg_min;
    const double eg2 = analysis_window.Eg_max;
    const double t1  = analysis_window.t_min;
    const double t2  = analysis_window.t_max;

    const int aw_line  = TColor::GetColor("#b2182b");
    const int aw_fill  = TColor::GetColor("#ef8a62");
    const int tsb_line = TColor::GetColor("#2166ac");
    const int tsb_fill = TColor::GetColor("#67a9cf");

    // TSB は塗りだけ入れ、共有境界 t=t_min/max は AW 側の線だけを見せる。
    DrawRegionBox(kTPlotMin, eg1, t1, eg2, tsb_fill, tsb_fill, 1, 0.10);
    DrawRegionBox(t1, eg1, t2, eg2, aw_fill, aw_fill, 1, 0.12);
    DrawRegionBox(t2, eg1, kTPlotMax, eg2, tsb_fill, tsb_fill, 1, 0.10);

    DrawLineSegment(kTPlotMin, eg1, kTPlotMin, eg2, tsb_line, 1, 2);
    DrawLineSegment(kTPlotMin, eg2, t1, eg2, tsb_line, 1, 2);
    DrawLineSegment(kTPlotMin, eg1, t1, eg1, tsb_line, 1, 2);
    DrawLineSegment(kTPlotMax, eg1, kTPlotMax, eg2, tsb_line, 1, 2);
    DrawLineSegment(t2, eg2, kTPlotMax, eg2, tsb_line, 1, 2);
    DrawLineSegment(t2, eg1, kTPlotMax, eg1, tsb_line, 1, 2);

    // AW の赤線は最後に描いて最前面に出す。
    DrawLineSegment(t1, eg1, t1, eg2, aw_line, 1, 2);
    DrawLineSegment(t2, eg1, t2, eg2, aw_line, 1, 2);
    DrawLineSegment(t1, eg2, t2, eg2, aw_line, 1, 2);
    DrawLineSegment(t1, eg1, t2, eg1, aw_line, 1, 2);

    TLatex lat;
    lat.SetTextFont(42);
    lat.SetTextSize(0.060);
    lat.SetTextAlign(22);
    lat.SetTextColor(tsb_line);
    lat.DrawLatex(0.5 * (kTPlotMin + t1), eg2 + 3.0, "TSB");
    lat.SetTextColor(aw_line);
    lat.DrawLatex(0.5 * (t1 + t2), eg2 + 3.0, "AW");
    lat.SetTextColor(tsb_line);
    lat.DrawLatex(0.5 * (t2 + kTPlotMax), eg2 + 3.0, "TSB");

}

static void DrawAnalysisWindowOverlayOnEeEg()
{
    const double ee1 = analysis_window.Ee_min;
    const double ee2 = analysis_window.Ee_max;
    const double eg1 = analysis_window.Eg_min;
    const double eg2 = analysis_window.Eg_max;

    const int aw_line = TColor::GetColor("#b2182b");
    const int aw_fill = TColor::GetColor("#ef8a62");

    DrawRegionBox(ee1, eg1, ee2, eg2, aw_fill, aw_fill, 1, 0.10);

    DrawLineSegment(ee1, eg1, ee1, eg2, aw_line, 1, 2);
    DrawLineSegment(ee2, eg1, ee2, eg2, aw_line, 1, 2);
    DrawLineSegment(ee1, eg2, ee2, eg2, aw_line, 1, 2);
    DrawLineSegment(ee1, eg1, ee2, eg1, aw_line, 1, 2);

    TLatex lat;
    lat.SetTextFont(42);
    lat.SetTextSize(0.060);
    lat.SetTextAlign(22);
    lat.SetTextColor(aw_line);
    lat.DrawLatex(0.5 * (ee1 + ee2), eg2 + 3.0, "AW");
}

static TH2D* MakeHist2D(const char* name, const char* title,
                        int nx, double x1, double x2,
                        int ny, double y1, double y2)
{
    TH2D* h = new TH2D(name, title, nx, x1, x2, ny, y1, y2);
    StyleHistogram2D(h);
    return h;
}

static TH2D* MakeHist2DVarX(const char* name, const char* title,
                            const std::vector<double>& xedges,
                            int ny, double y1, double y2)
{
    TH2D* h = new TH2D(name, title,
                       static_cast<int>(xedges.size()) - 1, xedges.data(),
                       ny, y1, y2);
    StyleHistogram2D(h);
    return h;
}

static TH2D* MakeHist2DVarXY(const char* name, const char* title,
                             const std::vector<double>& xedges,
                             const std::vector<double>& yedges)
{
    TH2D* h = new TH2D(name, title,
                       static_cast<int>(xedges.size()) - 1, xedges.data(),
                       static_cast<int>(yedges.size()) - 1, yedges.data());
    StyleHistogram2D(h);
    return h;
}

void plot_data_hist_thesis(
    const char* infile = "data/finaldata/step4f_run1to20_eg_doubleonly_allpatterns_500ns.txt")
{
    ApplyGlobalStyle();

    gSystem->mkdir("doc/fig", /*recursive=*/true);

    const TString base = MakeBaseName(infile);
    const TString out_paper = Form("doc/fig/%s_thesis_2d_overview.pdf", base.Data());
    const TString out_slide = Form("doc/fig/%s_slides_2d_focus.pdf", base.Data());

    const std::vector<double> phi_e_edges =
        Detector_PhiEdgesFromGrid(detres.phi_e_min, detres.phi_e_max, detres.N_phi_e);
    const std::vector<double> phi_g_edges =
        Detector_PhiEdgesFromGrid(detres.phi_g_min, detres.phi_g_max, detres.N_phi_g);
    const std::vector<double> th_edges = BuildThetaEdgesFromAllowedPhi(detres);

    TH2D* hEeEg = MakeHist2D(
        "hEeEg",
        ";E_{e^{+}} [MeV];E_{#gamma} [MeV]",
        100, kEePlotMin, kEePlotMax,
        100, kEgPlotMin, kEgPlotMax);
    TH2D* hThT = MakeHist2DVarX(
        "hThT",
        ";#theta_{e^{+}#gamma} [rad];t_{e^{+}#gamma} [ns]",
        th_edges, 100, kTPlotMin, kTPlotMax);
    TH2D* hTEg = MakeHist2D(
        "hTEg",
        ";t_{e^{+}#gamma} [ns];E_{#gamma} [MeV]",
        100, kTPlotMin, kTPlotMax,
        100, kEgPlotMin, kEgPlotMax);
    TH2D* hThEe = MakeHist2DVarX(
        "hThEe",
        ";#theta_{e^{+}#gamma} [rad];E_{e^{+}} [MeV]",
        th_edges, 100, kEePlotMin, kEePlotMax);
    TH2D* hThEg = MakeHist2DVarX(
        "hThEg",
        ";#theta_{e^{+}#gamma} [rad];E_{#gamma} [MeV]",
        th_edges, 100, kEgPlotMin, kEgPlotMax);
    TH2D* hPhi = MakeHist2DVarXY(
        "hPhi",
        ";#phi_{e^{+}} [rad];#phi_{#gamma} [rad]",
        phi_e_edges, phi_g_edges);

    std::ifstream fin(infile);
    if (!fin) {
        Error("plot_data_hist_thesis", "failed to open input file: %s", infile);
        return;
    }

    long long n_lines = 0;
    long long n_parsed = 0;
    long long n_used = 0;

    std::string line;
    while (std::getline(fin, line)) {
        ++n_lines;

        double Ee = 0.0;
        double Eg = 0.0;
        double t = 0.0;
        double phi_e = 0.0;
        double phi_g = 0.0;
        if (!ParseEventLine5Doubles(line, Ee, Eg, t, phi_e, phi_g)) continue;
        ++n_parsed;

        const int idx_e = Detector_PhiIndexFromValue(phi_e, detres.phi_e_min, detres.phi_e_max, detres.N_phi_e);
        const int idx_g = Detector_PhiIndexFromValue(phi_g, detres.phi_g_min, detres.phi_g_max, detres.N_phi_g);
        if (idx_e < 0 || idx_g < 0) continue;
        if (!Detector_IsAllowedPhiPairIndex(idx_e, idx_g, detres)) continue;

        const double phi_e_disc = Detector_PhiGridPoint(idx_e, detres.phi_e_min, detres.phi_e_max, detres.N_phi_e);
        const double phi_g_disc = Detector_PhiGridPoint(idx_g, detres.phi_g_min, detres.phi_g_max, detres.N_phi_g);
        const double theta = std::fabs(phi_e_disc - phi_g_disc);

        hEeEg->Fill(Ee, Eg);
        hThT->Fill(theta, t);
        hTEg->Fill(t, Eg);
        hThEe->Fill(theta, Ee);
        hThEg->Fill(theta, Eg);
        hPhi->Fill(phi_e_disc, phi_g_disc);
        ++n_used;
    }

    if (n_used == 0) {
        Warning("plot_data_hist_thesis", "no valid events were found. Nothing to draw.");
        return;
    }

    TCanvas cpaper("cpaper", "paper overview", 1800, 1100);
    cpaper.Divide(3, 2, 0.010, 0.010);

    cpaper.cd(1); SetPadMarginsForPaper(); hEeEg->Draw("colz"); DrawAnalysisWindowOverlayOnEeEg();
    cpaper.cd(2); SetPadMarginsForPaper(); hThT->Draw("colz");
    cpaper.cd(3); SetPadMarginsForPaper(); hTEg->Draw("colz");  DrawAwTsbOverlayOnTEg();
    cpaper.cd(4); SetPadMarginsForPaper(); hThEe->Draw("colz");
    cpaper.cd(5); SetPadMarginsForPaper(); hThEg->Draw("colz");
    cpaper.cd(6); SetPadMarginsForPaper(); hPhi->Draw("colz");
    cpaper.Print(out_paper.Data());

    TCanvas cslide("cslide", "slide focus", 1800, 720);
    cslide.Divide(2, 1, 0.030, 0.010);

    cslide.cd(1); SetPadMarginsForSlides(); hEeEg->Draw("colz"); DrawAnalysisWindowOverlayOnEeEg();
    cslide.cd(2); SetPadMarginsForSlides(); hTEg->Draw("colz"); DrawAwTsbOverlayOnTEg();
    cslide.Print(out_slide.Data());

    Info("plot_data_hist_thesis", "input=%s", infile);
    Info("plot_data_hist_thesis", "paper=%s", out_paper.Data());
    Info("plot_data_hist_thesis", "slide=%s", out_slide.Data());
    Info("plot_data_hist_thesis", "lines=%lld parsed=%lld used=%lld", n_lines, n_parsed, n_used);
}
