// scripts/plot_rmd_theory.cc
//
// RMD の理論式（スメア前）を解析窓内でモンテカルロ積分し、
// plot_data_hist と同様の 1D/2D ヒストを PDF に出力する。
//
// 使い方:
//   ./build/plot_rmd_theory 10000000 doc/rmd_theory.pdf
//
// 出力:
//   doc/rmd_theory_hist.pdf（3ページ）
//    1) メタ情報
//    2) 1D (Ee, Eg, phi_detector_e, phi_detector_g, theta_eg)
//    3) 2D (Ee,Eg), (theta_eg,Ee), (theta_eg,Eg), (phi_e,phi_g)
//
// 注意:
//  - 理論式は RMD_d6B_dEe_dEg_dOmegae_dOmegag を使用。
//  - エネルギーのスメアは行わない（真値のまま）。
//  - t は理論式に含まれないため解析窓で積分し、ヒスト表示はしない。
//  - phi は 0..pi の連続角として扱う（角度離散化はしない）。

#include <cerrno>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TError.h"
#include "TString.h"
#include "TLatex.h"
#include "TRandom3.h"

#include "p2meg/AnalysisWindow.h"
#include "p2meg/Constants.h"
#include "p2meg/DetectorResolution.h"
#include "p2meg/RMDSpectrum.h"

// ---- 固定ビン数（plot_data_hist に合わせる）----
static constexpr int kNBins_E   = 120; // Ee, Eg
static constexpr int kNBins_phi = 120; // phi_detector_e/g
static constexpr int kNBins_th  = 120; // theta_eg

static constexpr int kNBins2D_E   = 120;
static constexpr int kNBins2D_th  = 120;
static constexpr int kNBins2D_phi = 120;

// ---- MC 設定 ----
static constexpr long long kDefaultSamples = 2000000LL;
static constexpr unsigned long kSeed = 20260201UL;
static constexpr double kDMin = 1e-6;

static constexpr int kProgressBarWidth = 40;
static constexpr int kProgressIntervalMs = 200;

static bool IsFinite(double x) { return std::isfinite(x); }

static double Clamp(double x, double lo, double hi)
{
    if (x < lo) return lo;
    if (x > hi) return hi;
    return x;
}

static bool ParseLL(const char* s, long long& out)
{
    if (!s) return false;
    char* end = nullptr;
    errno = 0;
    const long long v = std::strtoll(s, &end, 10);
    if (errno != 0) return false;
    if (end == s || *end != '\0') return false;
    if (v <= 0) return false;
    out = v;
    return true;
}

static double ThetaFromPhi(double phi_e, double phi_g)
{
    if (!IsFinite(phi_e) || !IsFinite(phi_g)) return 0.0;
    const double pe = Clamp(phi_e, 0.0, pi);
    const double pg = Clamp(phi_g, 0.0, pi);
    return std::fabs(pe - pg);
}

// 理論式（d^6B/dEe dEg dOmegae dOmegag）から (Ee,Eg,phi_e,phi_g) への重み
static double RmdWeightNoSmear(double Ee, double Eg,
                               double phi_e, double phi_g,
                               double Pmu, double d_min)
{
    if (!(Ee > 0.0) || !(Eg > 0.0)) return 0.0;
    if (!IsFinite(phi_e) || !IsFinite(phi_g)) return 0.0;

    const double cosThetaE  = std::cos(phi_e);
    const double cosThetaG  = std::cos(phi_g);
    const double cosThetaEG = Clamp(std::cos(phi_e - phi_g), -1.0, 1.0);

    const double w0 = RMD_d6B_dEe_dEg_dOmegae_dOmegag(
        Ee, Eg, cosThetaEG, cosThetaE, cosThetaG, Pmu, d_min);

    if (!(w0 > 0.0) || !IsFinite(w0)) return 0.0;

    return w0;
}

static void PrintProgressBar(const char* label, long long done, long long total)
{
    if (total <= 0) return;
    const double frac = (done <= 0) ? 0.0
                                    : (done >= total ? 1.0
                                                     : (static_cast<double>(done) /
                                                        static_cast<double>(total)));
    const int filled = static_cast<int>(std::round(frac * kProgressBarWidth));
    std::ostringstream oss;
    oss << "\r" << label << " [";
    for (int i = 0; i < kProgressBarWidth; ++i) {
        oss << (i < filled ? '#' : '-');
    }
    oss << "] " << std::fixed << std::setprecision(1)
        << (frac * 100.0) << "% (" << done << "/" << total << ")";
    std::cerr << oss.str() << std::flush;
    if (done >= total) {
        std::cerr << "\n";
    }
}

static void DrawMetaPage(const char* outpdf,
                         long long n_samples,
                         long long n_theta_ok,
                         long long n_wpos,
                         double scale,
                         double sum_w_scaled)
{
    TLatex lat;
    lat.SetNDC(true);

    lat.SetTextAlign(13); // left-top
    lat.SetTextSize(0.040);
    lat.DrawLatex(0.05, 0.92, "p2MEG RMD theory (no smearing)");

    lat.SetTextSize(0.030);
    lat.DrawLatex(0.05, 0.84, Form("output : %s", outpdf));
    lat.DrawLatex(0.05, 0.79, Form("samples (MC) : %lld", n_samples));
    lat.DrawLatex(0.05, 0.74, Form("seed         : %lu", kSeed));
    lat.DrawLatex(0.05, 0.69, Form("P_mu         : %.6g", detres.P_mu));
    lat.DrawLatex(0.05, 0.64, Form("d_min        : %.3g", kDMin));
    lat.DrawLatex(0.05, 0.59, Form("theta in window : %lld", n_theta_ok));
    lat.DrawLatex(0.05, 0.54, Form("positive weight : %lld", n_wpos));
    lat.DrawLatex(0.05, 0.49, Form("MC scale (V/N)   : %.6g", scale));
    lat.DrawLatex(0.05, 0.44, Form("sum(weight*scale): %.6g", sum_w_scaled));

    lat.SetTextSize(0.032);
    lat.DrawLatex(0.05, 0.35, "Analysis window:");
    lat.SetTextSize(0.030);
    lat.DrawLatex(0.08, 0.30, Form("Ee [MeV]    : [%.6g, %.6g]", analysis_window.Ee_min, analysis_window.Ee_max));
    lat.DrawLatex(0.08, 0.26, Form("Eg [MeV]    : [%.6g, %.6g]", analysis_window.Eg_min, analysis_window.Eg_max));
    lat.DrawLatex(0.08, 0.22, Form("t  [ns]     : [%.6g, %.6g]", analysis_window.t_min,  analysis_window.t_max));
    lat.DrawLatex(0.08, 0.18, Form("theta_eg [rad] : [%.6g, %.6g]", analysis_window.theta_min, analysis_window.theta_max));

    lat.SetTextSize(0.028);
    lat.DrawLatex(0.05, 0.12, "t: uniform in window (integrated; not plotted).  phi: continuous [0, pi].");

    lat.SetTextSize(0.032);
    lat.DrawLatex(0.05, 0.06, "Pages: (1) meta  (2) 1D  (3) 2D");
}

int main(int argc, char** argv)
{
    long long n_samples = kDefaultSamples;
    std::string outpdf = "doc/rmd_theory_hist.pdf";

    if (argc == 2) {
        long long tmp = 0;
        if (ParseLL(argv[1], tmp)) {
            n_samples = tmp;
        } else {
            outpdf = argv[1];
        }
    } else if (argc == 3) {
        long long tmp = 0;
        if (!ParseLL(argv[1], tmp)) {
            std::cerr << "Usage: " << argv[0] << " [n_samples] [out_pdf]\n";
            return 1;
        }
        n_samples = tmp;
        outpdf = argv[2];
    } else if (argc > 3) {
        std::cerr << "Usage: " << argv[0] << " [n_samples] [out_pdf]\n";
        return 1;
    }

    gStyle->SetOptStat(0);
    gSystem->mkdir("doc", /*recursive=*/true);

    const double Ee_min = analysis_window.Ee_min;
    const double Ee_max = analysis_window.Ee_max;
    const double Eg_min = analysis_window.Eg_min;
    const double Eg_max = analysis_window.Eg_max;
    const double t_min  = analysis_window.t_min;
    const double t_max  = analysis_window.t_max;
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

    const double V_Ee  = Ee_max - Ee_min;
    const double V_Eg  = Eg_max - Eg_min;
    const double V_t   = t_max - t_min;
    const double V_phi = phi_max - phi_min;
    const double V_total = V_Ee * V_Eg * V_t * V_phi * V_phi;
    const double scale = V_total / static_cast<double>(n_samples);

    TRandom3 rng(kSeed);

    long long n_theta_ok = 0;
    long long n_wpos = 0;
    double sum_w = 0.0;

    auto last_progress = std::chrono::steady_clock::now();
    PrintProgressBar("rmd-theory", 0, n_samples);

    for (long long i = 0; i < n_samples; ++i) {
        const double Ee = rng.Uniform(Ee_min, Ee_max);
        const double Eg = rng.Uniform(Eg_min, Eg_max);
        const double phi_e = rng.Uniform(phi_min, phi_max);
        const double phi_g = rng.Uniform(phi_min, phi_max);

        const double theta_eg = ThetaFromPhi(phi_e, phi_g);
        if (theta_eg < th_min || theta_eg > th_max) {
            // 解析窓カット
        } else {
            ++n_theta_ok;
            const double w = RmdWeightNoSmear(Ee, Eg, phi_e, phi_g, detres.P_mu, kDMin);
            if (w > 0.0 && IsFinite(w)) {
                ++n_wpos;
                sum_w += w;
                const double wt = w * scale;

                // 1D
                hEe->Fill(Ee, wt);
                hEg->Fill(Eg, wt);
                hPhiE->Fill(phi_e, wt);
                hPhiG->Fill(phi_g, wt);
                hThEg->Fill(theta_eg, wt);

                // 2D
                h_EeEg->Fill(Ee, Eg, wt);
                h_ThEe->Fill(theta_eg, Ee, wt);
                h_ThEg->Fill(theta_eg, Eg, wt);
                h_PePg->Fill(phi_e, phi_g, wt);
            }
        }

        const auto now = std::chrono::steady_clock::now();
        const auto dt_ms = std::chrono::duration_cast<std::chrono::milliseconds>(now - last_progress).count();
        if (dt_ms >= kProgressIntervalMs || i + 1 == n_samples) {
            PrintProgressBar("rmd-theory", i + 1, n_samples);
            last_progress = now;
        }
    }

    if (n_wpos == 0) {
        Warning("plot_rmd_theory", "no positive-weight samples. Nothing to plot.");
        return 2;
    }

    const double sum_w_scaled = sum_w * scale;

    // ---- ページ1：メタ ----
    TCanvas c0("c0", "meta", 1200, 800);
    c0.cd();
    DrawMetaPage(outpdf.c_str(), n_samples, n_theta_ok, n_wpos, scale, sum_w_scaled);

    // ---- ページ2：1D ----
    TCanvas c1("c1", "1D theory", 1200, 800);
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

    TCanvas c2("c2", "2D theory", 1200, 800);
    c2.Divide(2, 2);
    c2.cd(1); Draw2D(h_EeEg);
    c2.cd(2); Draw2D(h_ThEe);
    c2.cd(3); Draw2D(h_ThEg);
    c2.cd(4); Draw2D(h_PePg);

    // ---- PDF（3ページ）----
    TString out = outpdf.c_str();
    c0.Print(Form("%s[", out.Data()));
    c0.Print(out.Data());
    c1.Print(out.Data());
    c2.Print(out.Data());
    c2.Print(Form("%s]", out.Data()));

    Info("plot_rmd_theory", "wrote: %s (3 pages)", out.Data());
    Info("plot_rmd_theory", "samples=%lld, theta_ok=%lld, wpos=%lld, sum_w_scaled=%.6g",
         n_samples, n_theta_ok, n_wpos, sum_w_scaled);

    return 0;
}
