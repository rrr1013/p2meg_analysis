// scripts/plot_fc_toy_components.cc
//
// 使い方:
//   ./build/plot_fc_toy_components
//       <datafile> <N_mu_eff_nom> <N_mu_eff_sigma> <BR_test>
//       <mean_vis_per_component> [seed] [out.pdf]
//
// 目的:
//   FC toy 生成で使っている PdfComponent / 棄却法 / AW 定義をそのまま使い、
//   Signal, RMD, ACC を別々に高統計で発生させて shape を可視化する。
//
// 注意:
//   - 図は「shape の可視化」が目的であり、1 本の toy 実験の事象数を
//     そのまま再現するものではない。
//   - ただし toy の発生アルゴリズム自体は FC と同じ
//     GenerateToyDatasetFromModel() を使う。

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TColor.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TSystem.h"

#include "p2meg/ACCGridPdf.h"
#include "p2meg/AnalysisWindow.h"
#include "p2meg/Constants.h"
#include "p2meg/DetectorResolution.h"
#include "p2meg/Event.h"
#include "p2meg/MathUtils.h"
#include "p2meg/PdfWrappers.h"
#include "p2meg/RMDGridPdf.h"
#include "p2meg/UpperLimit.h"

static const char* kDefaultRmdRoot = "data/pdf_cache/rmd_grid.root";
static const char* kDefaultRmdKey  = "rmd_grid";
static const char* kDefaultAccRoot = "data/pdf_cache/acc_grid.root";
static const char* kDefaultAccKey  = "acc_grid";

static constexpr int kColorSignal = kRed + 1;
static constexpr int kColorRmd = kBlue + 1;
static constexpr int kColorAcc = kGreen + 2;

static constexpr int kNBinsE = 40;
static constexpr int kNBinsT = 40;
static constexpr int kNBinsTheta = 24;

static bool ParseDoubleStrict(const char* s, double& out)
{
    if (!s) return false;
    char* end = nullptr;
    out = std::strtod(s, &end);
    if (end == s || *end != '\0') return false;
    return std::isfinite(out);
}

static bool ParseULLStrict(const char* s, unsigned long long& out)
{
    if (!s) return false;
    char* end = nullptr;
    out = std::strtoull(s, &end, 10);
    if (end == s || *end != '\0') return false;
    return true;
}

static bool ParseDoublesFromLine(const std::string& line, std::vector<double>& out)
{
    out.clear();
    std::istringstream iss(line);

    std::string first;
    if (!(iss >> first)) return false;
    if (!first.empty() && first[0] == '#') return false;

    char* endptr = nullptr;
    const double v0 = std::strtod(first.c_str(), &endptr);
    if (endptr == first.c_str() || *endptr != '\0') return false;

    out.push_back(v0);
    double v = 0.0;
    while (iss >> v) out.push_back(v);
    return true;
}

static bool ThetaFromPhiPair(const Event& ev, double& theta_out)
{
    int idx_e = -1;
    int idx_g = -1;
    if (!Detector_IsAllowedPhiPairValue(ev.phi_detector_e, ev.phi_detector_g,
                                        detres, idx_e, idx_g)) {
        return false;
    }

    const int N_phi_e = Math_GetNPhiE(detres);
    const int N_phi_g = Math_GetNPhiG(detres);
    const double phi_e =
        Detector_PhiGridPoint(idx_e, detres.phi_e_min, detres.phi_e_max, N_phi_e);
    const double phi_g =
        Detector_PhiGridPoint(idx_g, detres.phi_g_min, detres.phi_g_max, N_phi_g);
    theta_out = std::fabs(phi_e - phi_g);
    return std::isfinite(theta_out);
}

static bool IsInsideAnalysisWindow(const Event& ev)
{
    if (!std::isfinite(ev.Ee) || !std::isfinite(ev.Eg) ||
        !std::isfinite(ev.t) || !std::isfinite(ev.phi_detector_e) ||
        !std::isfinite(ev.phi_detector_g)) {
        return false;
    }

    double theta = 0.0;
    if (!ThetaFromPhiPair(ev, theta)) return false;

    if (ev.Ee < analysis_window.Ee_min || ev.Ee > analysis_window.Ee_max) return false;
    if (ev.Eg < analysis_window.Eg_min || ev.Eg > analysis_window.Eg_max) return false;
    if (ev.t  < analysis_window.t_min  || ev.t  > analysis_window.t_max ) return false;
    if (theta < analysis_window.theta_min || theta > analysis_window.theta_max) return false;
    return true;
}

static bool LoadEventsFromDat(const char* filepath,
                              std::vector<Event>& events,
                              long& n_skipped_non5,
                              long& n_skipped_outside)
{
    events.clear();
    n_skipped_non5 = 0;
    n_skipped_outside = 0;

    std::ifstream fin(filepath);
    if (!fin) return false;

    std::string line;
    std::vector<double> cols;
    while (std::getline(fin, line)) {
        if (!ParseDoublesFromLine(line, cols)) continue;
        if (cols.size() != 5) {
            ++n_skipped_non5;
            continue;
        }

        Event ev{};
        ev.Ee = cols[0];
        ev.Eg = cols[1];
        ev.t = cols[2];
        ev.phi_detector_e = cols[3];
        ev.phi_detector_g = cols[4];
        if (!IsInsideAnalysisWindow(ev)) {
            ++n_skipped_outside;
            continue;
        }
        events.push_back(ev);
    }

    return !events.empty();
}

static void ConvertHist1ToDensity(TH1D* h)
{
    if (!h) return;
    const double sum = h->Integral("width");
    if (!(sum > 0.0)) return;
    h->Scale(1.0 / sum);
}

static void ConvertHist2ToDensity(TH2D* h)
{
    if (!h) return;
    double sum = 0.0;
    for (int ix = 1; ix <= h->GetNbinsX(); ++ix) {
        const double wx = h->GetXaxis()->GetBinWidth(ix);
        for (int iy = 1; iy <= h->GetNbinsY(); ++iy) {
            const double wy = h->GetYaxis()->GetBinWidth(iy);
            sum += h->GetBinContent(ix, iy) * wx * wy;
        }
    }
    if (!(sum > 0.0)) return;
    h->Scale(1.0 / sum);
}

static void Style1D(TH1D* h, int color)
{
    if (!h) return;
    h->SetLineColor(color);
    h->SetMarkerColor(color);
    h->SetLineWidth(2);
    h->SetStats(0);
}

static void Style2D(TH2D* h)
{
    if (!h) return;
    h->SetStats(0);
}

static void FillToyHists(const std::vector<Event>& events,
                         TH1D* hEe,
                         TH1D* hEg,
                         TH1D* hT,
                         TH1D* hTheta,
                         TH2D* hEeEg,
                         TH2D* hTTheta,
                         TH2D* hPhi)
{
    for (const auto& ev : events) {
        double theta = 0.0;
        if (!ThetaFromPhiPair(ev, theta)) continue;
        hEe->Fill(ev.Ee);
        hEg->Fill(ev.Eg);
        hT->Fill(ev.t);
        hTheta->Fill(theta);
        hEeEg->Fill(ev.Ee, ev.Eg);
        hTTheta->Fill(ev.t, theta);
        hPhi->Fill(ev.phi_detector_e, ev.phi_detector_g);
    }

    ConvertHist1ToDensity(hEe);
    ConvertHist1ToDensity(hEg);
    ConvertHist1ToDensity(hT);
    ConvertHist1ToDensity(hTheta);
    ConvertHist2ToDensity(hEeEg);
    ConvertHist2ToDensity(hTTheta);
    ConvertHist2ToDensity(hPhi);
}

static void DrawMetaPage(const char* datafile,
                         double N_mu_eff_nom,
                         double N_mu_eff_sigma,
                         double BR_test,
                         double N_sig_nominal,
                         const FitResult& fit_prof,
                         double mean_vis_per_component,
                         const std::vector<double>& generated_counts)
{
    TLatex lat;
    lat.SetNDC(true);
    lat.SetTextAlign(13);

    lat.SetTextSize(0.040);
    lat.DrawLatex(0.06, 0.93, "FC toy component visualisation");

    lat.SetTextSize(0.030);
    lat.DrawLatex(0.06, 0.86, Form("datafile: %s", datafile));
    lat.DrawLatex(0.06, 0.81, Form("analysis window: Ee=[%.1f, %.1f] MeV, Eg=[%.1f, %.1f] MeV",
                                   analysis_window.Ee_min, analysis_window.Ee_max,
                                   analysis_window.Eg_min, analysis_window.Eg_max));
    lat.DrawLatex(0.06, 0.76, Form("analysis window: t=[%.1f, %.1f] ns, theta=[%.4f, %.4f] rad",
                                   analysis_window.t_min, analysis_window.t_max,
                                   analysis_window.theta_min, analysis_window.theta_max));
    lat.DrawLatex(0.06, 0.71, Form("detres: sigma_t=%.3f ns, t_mean=%.3f ns, P_mu=%.3f",
                                   detres.sigma_t, detres.t_mean, detres.P_mu));

    lat.DrawLatex(0.06, 0.63, Form("FC example point: BR_test = %.10g", BR_test));
    lat.DrawLatex(0.06, 0.58, Form("N_mu_eff = %.10f +/- %.10f", N_mu_eff_nom, N_mu_eff_sigma));
    lat.DrawLatex(0.06, 0.53, Form("N_sig_test_nominal = BR_test * N_mu_eff = %.6f", N_sig_nominal));
    if (fit_prof.yields_hat.size() >= 3) {
        lat.DrawLatex(0.06, 0.48, Form("fixed-BR profile mean yields: signal=%.6f, RMD=%.6f, ACC=%.6f",
                                       fit_prof.yields_hat[0], fit_prof.yields_hat[1], fit_prof.yields_hat[2]));
    }

    lat.DrawLatex(0.06, 0.39, Form("visualisation sample: mean %.0f events for each component separately",
                                   mean_vis_per_component));
    if (generated_counts.size() >= 3) {
        lat.DrawLatex(0.06, 0.34, Form("generated counts: signal=%.0f, RMD=%.0f, ACC=%.0f",
                                       generated_counts[0], generated_counts[1], generated_counts[2]));
    }

    lat.DrawLatex(0.06, 0.24, "Generation algorithm:");
    lat.DrawLatex(0.10, 0.19, "1) choose N_k from Poisson(mean_k)");
    lat.DrawLatex(0.10, 0.14, "2) sample events in AW by accept-reject from current PdfComponent");
    lat.DrawLatex(0.10, 0.09, "3) here components are generated one by one to show their shapes");
}

int main(int argc, char** argv)
{
    if (argc != 6 && argc != 7 && argc != 8) {
        std::cerr << "Usage: " << argv[0]
                  << " <datafile> <N_mu_eff_nom> <N_mu_eff_sigma> <BR_test>"
                  << " <mean_vis_per_component> [seed] [out.pdf]\n";
        return 1;
    }

    const char* datafile = argv[1];
    double N_mu_eff_nom = 0.0;
    double N_mu_eff_sigma = 0.0;
    double BR_test = 0.0;
    double mean_vis_per_component = 0.0;
    unsigned long long seed = 97531ULL;
    std::string outpdf = "doc/finalanalysis/fc_toy_components_br000085.pdf";

    if (!ParseDoubleStrict(argv[2], N_mu_eff_nom) ||
        !ParseDoubleStrict(argv[3], N_mu_eff_sigma) ||
        !ParseDoubleStrict(argv[4], BR_test) ||
        !ParseDoubleStrict(argv[5], mean_vis_per_component)) {
        std::cerr << "[plot_fc_toy_components] failed to parse arguments\n";
        return 1;
    }
    if (argc >= 7 && !ParseULLStrict(argv[6], seed)) {
        outpdf = argv[6];
    }
    if (argc == 8) {
        seed = std::strtoull(argv[6], nullptr, 10);
        outpdf = argv[7];
    }

    if (!(N_mu_eff_nom > 0.0) || !(mean_vis_per_component > 0.0) || !(BR_test >= 0.0)) {
        std::cerr << "[plot_fc_toy_components] invalid numeric input\n";
        return 1;
    }

    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);
    gSystem->mkdir("doc/finalanalysis", true);

    std::vector<Event> data_events;
    long n_skipped_non5 = 0;
    long n_skipped_outside = 0;
    if (!LoadEventsFromDat(datafile, data_events, n_skipped_non5, n_skipped_outside)) {
        std::cerr << "[plot_fc_toy_components] failed to load events from " << datafile << "\n";
        return 2;
    }

    if (!RMDGridPdf_Load(kDefaultRmdRoot, kDefaultRmdKey)) {
        std::cerr << "[plot_fc_toy_components] RMDGridPdf_Load failed\n";
        return 3;
    }
    if (!ACCGridPdf_Load(kDefaultAccRoot, kDefaultAccKey)) {
        std::cerr << "[plot_fc_toy_components] ACCGridPdf_Load failed\n";
        return 3;
    }

    static SignalPdfContext sigctx{
        analysis_window,
        detres,
        kMassesPDG
    };

    std::vector<PdfComponent> components;
    components.push_back(MakeSignalComponent(&sigctx));
    components.push_back(MakeRMDComponent());
    components.push_back(MakeACCComponent());

    const double N0 = static_cast<double>(data_events.size());
    FitConfig free_cfg;
    free_cfg.start_yields = {N0 / 3.0, N0 / 3.0, N0 / 3.0};
    free_cfg.max_calls = 20000;
    free_cfg.tol = 1e-3;

    FitConfig prof_cfg = free_cfg;

    FitResult fit_free;
    FitResult fit_prof;
    const double N_sig_nominal = BR_test * N_mu_eff_nom;
    const double q_obs = EvaluateProfileLikelihoodQ(data_events, components,
                                                    free_cfg, prof_cfg,
                                                    N_sig_nominal,
                                                    fit_free, fit_prof);
    if (!std::isfinite(q_obs) || fit_prof.status != 0) {
        std::cerr << "[plot_fc_toy_components] fixed-BR profile fit failed\n";
        return 4;
    }

    ToyGeneratorConfig toy_cfg;
    toy_cfg.seed = seed;
    // 可視化用なので pmax 推定は FC 本番より軽くする。
    toy_cfg.pmax_scan_trials = 5000;
    toy_cfg.pmax_safety = 5.0;
    toy_cfg.pmax_update = 1.2;

    const double theta_axis_max = Math_AxisMaxInclusive(analysis_window.theta_max);
    const double phi_e_axis_max = Detector_PhiAxisMaxInclusive(detres.phi_e_max);
    const double phi_g_axis_max = Detector_PhiAxisMaxInclusive(detres.phi_g_max);

    TH1D* hSigEe = new TH1D("hSigEe", "Signal;E_{e} [MeV];Density", kNBinsE,
                            analysis_window.Ee_min, analysis_window.Ee_max);
    TH1D* hRmdEe = new TH1D("hRmdEe", "RMD;E_{e} [MeV];Density", kNBinsE,
                            analysis_window.Ee_min, analysis_window.Ee_max);
    TH1D* hAccEe = new TH1D("hAccEe", "ACC;E_{e} [MeV];Density", kNBinsE,
                            analysis_window.Ee_min, analysis_window.Ee_max);

    TH1D* hSigEg = new TH1D("hSigEg", "Signal;E_{#gamma} [MeV];Density", kNBinsE,
                            analysis_window.Eg_min, analysis_window.Eg_max);
    TH1D* hRmdEg = new TH1D("hRmdEg", "RMD;E_{#gamma} [MeV];Density", kNBinsE,
                            analysis_window.Eg_min, analysis_window.Eg_max);
    TH1D* hAccEg = new TH1D("hAccEg", "ACC;E_{#gamma} [MeV];Density", kNBinsE,
                            analysis_window.Eg_min, analysis_window.Eg_max);

    TH1D* hSigT = new TH1D("hSigT", "Signal;#Deltat [ns];Density", kNBinsT,
                           analysis_window.t_min, analysis_window.t_max);
    TH1D* hRmdT = new TH1D("hRmdT", "RMD;#Deltat [ns];Density", kNBinsT,
                           analysis_window.t_min, analysis_window.t_max);
    TH1D* hAccT = new TH1D("hAccT", "ACC;#Deltat [ns];Density", kNBinsT,
                           analysis_window.t_min, analysis_window.t_max);

    TH1D* hSigTheta = new TH1D("hSigTheta", "Signal;#theta_{eg} [rad];Density", kNBinsTheta,
                               analysis_window.theta_min, theta_axis_max);
    TH1D* hRmdTheta = new TH1D("hRmdTheta", "RMD;#theta_{eg} [rad];Density", kNBinsTheta,
                               analysis_window.theta_min, theta_axis_max);
    TH1D* hAccTheta = new TH1D("hAccTheta", "ACC;#theta_{eg} [rad];Density", kNBinsTheta,
                               analysis_window.theta_min, theta_axis_max);

    TH2D* h2SigEeEg = new TH2D("h2SigEeEg", "Signal;E_{e} [MeV];E_{#gamma} [MeV]",
                               kNBinsE, analysis_window.Ee_min, analysis_window.Ee_max,
                               kNBinsE, analysis_window.Eg_min, analysis_window.Eg_max);
    TH2D* h2RmdEeEg = new TH2D("h2RmdEeEg", "RMD;E_{e} [MeV];E_{#gamma} [MeV]",
                               kNBinsE, analysis_window.Ee_min, analysis_window.Ee_max,
                               kNBinsE, analysis_window.Eg_min, analysis_window.Eg_max);
    TH2D* h2AccEeEg = new TH2D("h2AccEeEg", "ACC;E_{e} [MeV];E_{#gamma} [MeV]",
                               kNBinsE, analysis_window.Ee_min, analysis_window.Ee_max,
                               kNBinsE, analysis_window.Eg_min, analysis_window.Eg_max);

    TH2D* h2SigTTheta = new TH2D("h2SigTTheta", "Signal;#Deltat [ns];#theta_{eg} [rad]",
                                 kNBinsT, analysis_window.t_min, analysis_window.t_max,
                                 kNBinsTheta, analysis_window.theta_min, theta_axis_max);
    TH2D* h2RmdTTheta = new TH2D("h2RmdTTheta", "RMD;#Deltat [ns];#theta_{eg} [rad]",
                                 kNBinsT, analysis_window.t_min, analysis_window.t_max,
                                 kNBinsTheta, analysis_window.theta_min, theta_axis_max);
    TH2D* h2AccTTheta = new TH2D("h2AccTTheta", "ACC;#Deltat [ns];#theta_{eg} [rad]",
                                 kNBinsT, analysis_window.t_min, analysis_window.t_max,
                                 kNBinsTheta, analysis_window.theta_min, theta_axis_max);

    TH2D* h2SigPhi = new TH2D("h2SigPhi", "Signal;#phi_{e} [rad];#phi_{#gamma} [rad]",
                              detres.N_phi_e + 1, detres.phi_e_min, phi_e_axis_max,
                              detres.N_phi_g + 1, detres.phi_g_min, phi_g_axis_max);
    TH2D* h2RmdPhi = new TH2D("h2RmdPhi", "RMD;#phi_{e} [rad];#phi_{#gamma} [rad]",
                              detres.N_phi_e + 1, detres.phi_e_min, phi_e_axis_max,
                              detres.N_phi_g + 1, detres.phi_g_min, phi_g_axis_max);
    TH2D* h2AccPhi = new TH2D("h2AccPhi", "ACC;#phi_{e} [rad];#phi_{#gamma} [rad]",
                              detres.N_phi_e + 1, detres.phi_e_min, phi_e_axis_max,
                              detres.N_phi_g + 1, detres.phi_g_min, phi_g_axis_max);

    std::vector<double> generated_counts(3, 0.0);
    std::vector<std::vector<Event>> vis_events(3);
    std::vector<double> pmax_cache(3, 0.0);
    for (std::size_t k = 0; k < 3; ++k) {
        std::vector<double> mean_yields(3, 0.0);
        mean_yields[k] = mean_vis_per_component;
        std::vector<double> gen_yields;
        if (!GenerateToyDatasetFromModel(components, mean_yields, toy_cfg,
                                         static_cast<unsigned long long>(1000 + k),
                                         vis_events[k], &pmax_cache, &gen_yields)) {
            std::cerr << "[plot_fc_toy_components] toy generation failed for component " << k << "\n";
            return 5;
        }
        if (gen_yields.size() == 3) generated_counts[k] = gen_yields[k];
    }

    FillToyHists(vis_events[0], hSigEe, hSigEg, hSigT, hSigTheta, h2SigEeEg, h2SigTTheta, h2SigPhi);
    FillToyHists(vis_events[1], hRmdEe, hRmdEg, hRmdT, hRmdTheta, h2RmdEeEg, h2RmdTTheta, h2RmdPhi);
    FillToyHists(vis_events[2], hAccEe, hAccEg, hAccT, hAccTheta, h2AccEeEg, h2AccTTheta, h2AccPhi);

    Style1D(hSigEe, kColorSignal); Style1D(hRmdEe, kColorRmd); Style1D(hAccEe, kColorAcc);
    Style1D(hSigEg, kColorSignal); Style1D(hRmdEg, kColorRmd); Style1D(hAccEg, kColorAcc);
    Style1D(hSigT, kColorSignal);  Style1D(hRmdT, kColorRmd);  Style1D(hAccT, kColorAcc);
    Style1D(hSigTheta, kColorSignal); Style1D(hRmdTheta, kColorRmd); Style1D(hAccTheta, kColorAcc);
    Style2D(h2SigEeEg); Style2D(h2RmdEeEg); Style2D(h2AccEeEg);
    Style2D(h2SigTTheta); Style2D(h2RmdTTheta); Style2D(h2AccTTheta);
    Style2D(h2SigPhi); Style2D(h2RmdPhi); Style2D(h2AccPhi);

    TCanvas c("c", "fc toy components", 1200, 900);
    c.Print((outpdf + "[").c_str());

    c.Clear();
    DrawMetaPage(datafile, N_mu_eff_nom, N_mu_eff_sigma, BR_test,
                 N_sig_nominal, fit_prof, mean_vis_per_component, generated_counts);
    c.Print(outpdf.c_str());

    c.Clear();
    c.Divide(2, 2);
    c.cd(1);
    hSigEe->Draw("hist");
    hRmdEe->Draw("hist same");
    hAccEe->Draw("hist same");
    {
        TLegend leg(0.58, 0.68, 0.88, 0.88);
        leg.AddEntry(hSigEe, "signal toy", "l");
        leg.AddEntry(hRmdEe, "RMD toy", "l");
        leg.AddEntry(hAccEe, "ACC toy", "l");
        leg.Draw();
    }
    c.cd(2);
    hSigEg->Draw("hist");
    hRmdEg->Draw("hist same");
    hAccEg->Draw("hist same");
    c.cd(3);
    hSigT->Draw("hist");
    hRmdT->Draw("hist same");
    hAccT->Draw("hist same");
    c.cd(4);
    hSigTheta->Draw("hist");
    hRmdTheta->Draw("hist same");
    hAccTheta->Draw("hist same");
    c.Print(outpdf.c_str());

    c.Clear();
    c.Divide(2, 2);
    c.cd(1); h2SigEeEg->Draw("colz");
    c.cd(2); h2SigTTheta->Draw("colz");
    c.cd(3); h2SigPhi->Draw("colz text");
    c.cd(4);
    {
        TLatex lat;
        lat.SetNDC(true);
        lat.SetTextAlign(13);
        lat.SetTextSize(0.038);
        lat.DrawLatex(0.08, 0.90, "Signal toy");
        lat.SetTextSize(0.030);
        lat.DrawLatex(0.08, 0.78, "Ee-Eg: endpoint peak near the kinematic edge");
        lat.DrawLatex(0.08, 0.70, "t-theta: narrow in time, concentrated near back-to-back");
        lat.DrawLatex(0.08, 0.62, "phi pairs: detector mask gives only allowed discrete cells");
    }
    c.Print(outpdf.c_str());

    c.Clear();
    c.Divide(2, 2);
    c.cd(1); h2RmdEeEg->Draw("colz");
    c.cd(2); h2RmdTTheta->Draw("colz");
    c.cd(3); h2RmdPhi->Draw("colz text");
    c.cd(4);
    {
        TLatex lat;
        lat.SetNDC(true);
        lat.SetTextAlign(13);
        lat.SetTextSize(0.038);
        lat.DrawLatex(0.08, 0.90, "RMD toy");
        lat.SetTextSize(0.030);
        lat.DrawLatex(0.08, 0.78, "Ee-Eg: broad along the allowed phase space");
        lat.DrawLatex(0.08, 0.70, "t-theta: prompt in time, but angle/energy broader than signal");
        lat.DrawLatex(0.08, 0.62, "phi pairs: same detector acceptance mask as signal");
    }
    c.Print(outpdf.c_str());

    c.Clear();
    c.Divide(2, 2);
    c.cd(1); h2AccEeEg->Draw("colz");
    c.cd(2); h2AccTTheta->Draw("colz");
    c.cd(3); h2AccPhi->Draw("colz text");
    c.cd(4);
    {
        TLatex lat;
        lat.SetNDC(true);
        lat.SetTextAlign(13);
        lat.SetTextSize(0.038);
        lat.DrawLatex(0.08, 0.90, "ACC toy");
        lat.SetTextSize(0.030);
        lat.DrawLatex(0.08, 0.78, "Ee-Eg: current ACC grid shape in AW");
        lat.DrawLatex(0.08, 0.70, "t-theta: fit-based time shape times angular acceptance");
        lat.DrawLatex(0.08, 0.62, "phi pairs: acceptance-weighted discrete detector cells");
    }
    c.Print(outpdf.c_str());

    c.Print((outpdf + "]").c_str());

    std::cout << std::setprecision(10);
    std::cout << "[plot_fc_toy_components] wrote " << outpdf << "\n";
    std::cout << "[plot_fc_toy_components] q_obs = " << q_obs << "\n";
    std::cout << "[plot_fc_toy_components] N_sig_test_nominal = " << N_sig_nominal << "\n";
    if (fit_prof.yields_hat.size() >= 3) {
        std::cout << "[plot_fc_toy_components] fixed-BR profile means: "
                  << fit_prof.yields_hat[0] << " "
                  << fit_prof.yields_hat[1] << " "
                  << fit_prof.yields_hat[2] << "\n";
    }
    std::cout << "[plot_fc_toy_components] generated counts: "
              << generated_counts[0] << " "
              << generated_counts[1] << " "
              << generated_counts[2] << "\n";

    return 0;
}
