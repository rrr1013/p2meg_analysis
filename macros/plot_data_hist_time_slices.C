// macros/plot_data_hist_time_slices.C
//
// 入力: .txt/.dat テキスト。
//      「5列すべてが double として読める行」だけを Event として扱う。
//      列順: Ee Eg t phi_detector_e phi_detector_g
//
// 目的:
//   - 全時間範囲 [-500, 500] ns を 10 分割し、
//     各時間窓ごとに plot_data_hist と同様の 1D ヒストを作る。
//   - 出力は 1 つの PDF にまとめる。
//
// 出力:
//   doc/finalanalysis/data_hist_time_slices_<入力ファイル名(拡張子除く)>.pdf
//
// 実行例:
//   root -l -q 'macros/plot_data_hist_time_slices.C("data/finaldata/step4f_run1to20_eg_doubleonly_allpatterns_500ns.txt")'
//

R__ADD_INCLUDE_PATH(./include)

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

#include "TCanvas.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TError.h"
#include "TString.h"
#include "TLatex.h"
#include "TPad.h"

#include "p2meg/Event.h"
#include "p2meg/DetectorResolution.h"
#include "p2meg/MathUtils.h"

// ---- 固定ビン数（plot_data_hist に合わせる）----
static constexpr int kNBins_E   = 120; // Ee, Eg
static constexpr int kNBins_t   = 160; // t
static constexpr int kNBins_phi = 120; // phi_detector_e/g
static constexpr int kNBins_th  = 120; // theta_eg

// ---- 時間分割設定 ----
static constexpr int kNTimeSlices = 10;
static constexpr double kTimeAllMin = -500.0; // [ns]
static constexpr double kTimeAllMax =  500.0; // [ns]

static bool ParseEventLine5Doubles(const std::string& line,
                                   double& Epos, double& Egam, double& dt,
                                   double& phi_pos, double& phi_gam)
{
    if (line.empty()) return false;

    std::istringstream iss(line);
    if (!(iss >> Epos >> Egam >> dt >> phi_pos >> phi_gam)) {
        return false;
    }
    if (!std::isfinite(Epos) || !std::isfinite(Egam) || !std::isfinite(dt) ||
        !std::isfinite(phi_pos) || !std::isfinite(phi_gam)) {
        return false;
    }
    return true;
}

static TString MakeOutputPdfPath(const char* infile)
{
    TString base = gSystem->BaseName(infile);
    const Ssiz_t dot = base.Last('.');
    if (dot != kNPOS) base.Remove(dot);

    gSystem->mkdir("doc/finalanalysis", /*recursive=*/true);
    return Form("doc/finalanalysis/data_hist_time_slices_%s.pdf", base.Data());
}

static int FindTimeSliceIndex(double t)
{
    if (!std::isfinite(t)) return -1;
    if (t < kTimeAllMin || t > kTimeAllMax) return -1;

    const double width = (kTimeAllMax - kTimeAllMin) / static_cast<double>(kNTimeSlices);
    if (!(width > 0.0)) return -1;

    if (t == kTimeAllMax) return kNTimeSlices - 1;

    const int idx = static_cast<int>((t - kTimeAllMin) / width);
    if (idx < 0 || idx >= kNTimeSlices) return -1;
    return idx;
}

static double SliceTimeMin(int idx)
{
    const double width = (kTimeAllMax - kTimeAllMin) / static_cast<double>(kNTimeSlices);
    return kTimeAllMin + width * static_cast<double>(idx);
}

static double SliceTimeMax(int idx)
{
    const double width = (kTimeAllMax - kTimeAllMin) / static_cast<double>(kNTimeSlices);
    if (idx == kNTimeSlices - 1) return kTimeAllMax;
    return kTimeAllMin + width * static_cast<double>(idx + 1);
}

static void DrawMetaPage(const char* infile,
                         const char* outpdf,
                         long long n_lines,
                         long long n_parsed,
                         long long n_phi_ok,
                         long long n_time_ok,
                         const std::vector<long long>& slice_counts)
{
    TLatex lat;
    lat.SetNDC(true);
    lat.SetTextAlign(13);

    lat.SetTextSize(0.040);
    lat.DrawLatex(0.05, 0.94, "p2MEG time-sliced histograms");

    lat.SetTextSize(0.030);
    lat.DrawLatex(0.05, 0.87, Form("input  : %s", infile));
    lat.DrawLatex(0.05, 0.82, Form("output : %s", outpdf));

    lat.DrawLatex(0.05, 0.74, Form("lines read         : %lld", n_lines));
    lat.DrawLatex(0.05, 0.69, Form("parsed (5 doubles) : %lld", n_parsed));
    lat.DrawLatex(0.05, 0.64, Form("phi-accepted       : %lld", n_phi_ok));
    lat.DrawLatex(0.05, 0.59, Form("time in [-500,500] : %lld", n_time_ok));

    lat.DrawLatex(0.05, 0.51, "Time slices:");
    lat.SetTextSize(0.027);
    for (int i = 0; i < kNTimeSlices; ++i) {
        const double x = (i < 5) ? 0.08 : 0.53;
        const double y = (i < 5)
            ? (0.46 - 0.07 * static_cast<double>(i))
            : (0.46 - 0.07 * static_cast<double>(i - 5));
        const double tmin = SliceTimeMin(i);
        const double tmax = SliceTimeMax(i);
        const char* bracket = (i == kNTimeSlices - 1) ? "]" : ")";
        lat.DrawLatex(x, y,
                      Form("%2d: [%.0f, %.0f%s  entries = %lld",
                           i + 1, tmin, tmax, bracket, slice_counts[static_cast<size_t>(i)]));
    }

    lat.SetTextSize(0.026);
    lat.DrawLatex(0.05, 0.08, "Each following page shows Ee, Eg, t, phi_e, phi_g, theta_eg for one time slice.");
}

void plot_data_hist_time_slices(
    const char* infile = "data/finaldata/step4f_run1to20_eg_doubleonly_allpatterns_500ns.txt")
{
    const TString outpdf = MakeOutputPdfPath(infile);

    gStyle->SetOptStat(0);

    const double pi_val = 3.14159265358979323846;
    const double Ee_min = 0.0;
    const double Ee_max = 70.0;
    const double Eg_min = 0.0;
    const double Eg_max = 70.0;
    const double t_min  = kTimeAllMin;
    const double t_max  = kTimeAllMax;

    double th_plot_min = 0.0;
    double th_plot_max = pi_val;
    double th_det_min = 0.0;
    double th_det_max = 0.0;
    if (Detector_ThetaRangeFromAllowedPhi(detres, th_det_min, th_det_max)) {
        th_plot_min = th_det_min;
        th_plot_max = th_det_max;
    }
    const double th_plot_min_axis = -Math_AxisMaxInclusive(-th_plot_min);
    const double th_plot_max_axis = Math_AxisMaxInclusive(th_plot_max);

    const double phi_e_min = detres.phi_e_min;
    const double phi_e_max = detres.phi_e_max;
    const double phi_g_min = detres.phi_g_min;
    const double phi_g_max = detres.phi_g_max;
    const double phi_e_max_plot = Detector_PhiAxisMaxInclusive(phi_e_max);
    const double phi_g_max_plot = Detector_PhiAxisMaxInclusive(phi_g_max);

    std::vector<TH1D*> hEe;
    std::vector<TH1D*> hEg;
    std::vector<TH1D*> ht;
    std::vector<TH1D*> hPhiE;
    std::vector<TH1D*> hPhiG;
    std::vector<TH1D*> hThEg;
    std::vector<long long> slice_counts(static_cast<size_t>(kNTimeSlices), 0LL);

    hEe.reserve(kNTimeSlices);
    hEg.reserve(kNTimeSlices);
    ht.reserve(kNTimeSlices);
    hPhiE.reserve(kNTimeSlices);
    hPhiG.reserve(kNTimeSlices);
    hThEg.reserve(kNTimeSlices);

    for (int i = 0; i < kNTimeSlices; ++i) {
        hEe.push_back(new TH1D(Form("hEe_%d", i), "Ee;Ee [MeV];Entries", kNBins_E, Ee_min, Ee_max));
        hEg.push_back(new TH1D(Form("hEg_%d", i), "Eg;Eg [MeV];Entries", kNBins_E, Eg_min, Eg_max));
        ht.push_back(new TH1D(Form("ht_%d", i), "t;t [ns];Entries", kNBins_t, t_min, t_max));
        hPhiE.push_back(new TH1D(Form("hPhiE_%d", i),
                                 "phi_{detector,e};phi_{detector,e} [rad];Entries",
                                 kNBins_phi, phi_e_min, phi_e_max_plot));
        hPhiG.push_back(new TH1D(Form("hPhiG_%d", i),
                                 "phi_{detector,#gamma};phi_{detector,#gamma} [rad];Entries",
                                 kNBins_phi, phi_g_min, phi_g_max_plot));
        hThEg.push_back(new TH1D(Form("hThEg_%d", i),
                                 "theta_{eg};theta_{eg} [rad];Entries",
                                 kNBins_th, th_plot_min_axis, th_plot_max_axis));
    }

    std::ifstream fin(infile);
    if (!fin) {
        Error("plot_data_hist_time_slices", "failed to open input file: %s", infile);
        return;
    }

    long long n_lines  = 0;
    long long n_parsed = 0;
    long long n_phi_ok = 0;
    long long n_time_ok = 0;
    long long n_phi_out = 0;
    long long n_time_out = 0;

    std::string line;
    while (std::getline(fin, line)) {
        ++n_lines;

        double Epos = 0.0;
        double Egam = 0.0;
        double dt = 0.0;
        double phi_pos = 0.0;
        double phi_gam = 0.0;
        if (!ParseEventLine5Doubles(line, Epos, Egam, dt, phi_pos, phi_gam)) {
            continue;
        }
        ++n_parsed;

        Event ev;
        ev.Ee = Epos;
        ev.Eg = Egam;
        ev.t = dt;
        ev.phi_detector_e = phi_pos;
        ev.phi_detector_g = phi_gam;

        const int idx_e = Detector_PhiIndexFromValue(ev.phi_detector_e,
                                                     detres.phi_e_min, detres.phi_e_max,
                                                     Math_GetNPhiE(detres));
        const int idx_g = Detector_PhiIndexFromValue(ev.phi_detector_g,
                                                     detres.phi_g_min, detres.phi_g_max,
                                                     Math_GetNPhiG(detres));
        if (idx_e < 0 || idx_g < 0 || !Detector_IsAllowedPhiPairIndex(idx_e, idx_g, detres)) {
            ++n_phi_out;
            continue;
        }
        ++n_phi_ok;

        const int slice_idx = FindTimeSliceIndex(ev.t);
        if (slice_idx < 0) {
            ++n_time_out;
            continue;
        }
        ++n_time_ok;

        const double phi_e_disc = Detector_PhiGridPoint(idx_e, detres.phi_e_min, detres.phi_e_max,
                                                        Math_GetNPhiE(detres));
        const double phi_g_disc = Detector_PhiGridPoint(idx_g, detres.phi_g_min, detres.phi_g_max,
                                                        Math_GetNPhiG(detres));
        const double theta_eg = std::fabs(phi_e_disc - phi_g_disc);

        hEe[static_cast<size_t>(slice_idx)]->Fill(ev.Ee);
        hEg[static_cast<size_t>(slice_idx)]->Fill(ev.Eg);
        ht[static_cast<size_t>(slice_idx)]->Fill(ev.t);
        hPhiE[static_cast<size_t>(slice_idx)]->Fill(ev.phi_detector_e);
        hPhiG[static_cast<size_t>(slice_idx)]->Fill(ev.phi_detector_g);
        hThEg[static_cast<size_t>(slice_idx)]->Fill(theta_eg);
        ++slice_counts[static_cast<size_t>(slice_idx)];
    }

    if (n_time_ok == 0) {
        Warning("plot_data_hist_time_slices", "no events in the requested time range. Nothing to plot.");
        return;
    }

    TCanvas cmeta("cmeta", "meta", 1200, 800);
    cmeta.cd();
    DrawMetaPage(infile, outpdf.Data(), n_lines, n_parsed, n_phi_ok, n_time_ok, slice_counts);

    cmeta.Print(Form("%s[", outpdf.Data()));
    cmeta.Print(outpdf.Data());

    TLatex page_label;
    page_label.SetNDC(true);
    page_label.SetTextAlign(13);
    page_label.SetTextSize(0.032);

    for (int i = 0; i < kNTimeSlices; ++i) {
        TCanvas c(Form("c_slice_%d", i), Form("slice_%d", i), 1200, 820);
        c.Divide(3, 2);

        c.cd(1); gPad->SetGrid(); hEe[static_cast<size_t>(i)]->SetLineWidth(2); hEe[static_cast<size_t>(i)]->Draw("hist");
        c.cd(2); gPad->SetGrid(); hEg[static_cast<size_t>(i)]->SetLineWidth(2); hEg[static_cast<size_t>(i)]->Draw("hist");
        c.cd(3); gPad->SetGrid(); ht[static_cast<size_t>(i)]->SetLineWidth(2);  ht[static_cast<size_t>(i)]->SetMinimum(0.0); ht[static_cast<size_t>(i)]->Draw("hist");
        c.cd(4); gPad->SetGrid(); hPhiE[static_cast<size_t>(i)]->SetLineWidth(2); hPhiE[static_cast<size_t>(i)]->Draw("hist");
        c.cd(5); gPad->SetGrid(); hPhiG[static_cast<size_t>(i)]->SetLineWidth(2); hPhiG[static_cast<size_t>(i)]->Draw("hist");
        c.cd(6); gPad->SetGrid(); hThEg[static_cast<size_t>(i)]->SetLineWidth(2); hThEg[static_cast<size_t>(i)]->Draw("hist");

        c.cd();
        page_label.DrawLatex(0.05, 0.98,
                             Form("slice %d: [%.0f, %.0f%s   entries = %lld",
                                  i + 1,
                                  SliceTimeMin(i),
                                  SliceTimeMax(i),
                                  (i == kNTimeSlices - 1) ? "]" : ")",
                                  slice_counts[static_cast<size_t>(i)]));

        c.Print(outpdf.Data());
    }

    cmeta.Print(Form("%s]", outpdf.Data()));

    Info("plot_data_hist_time_slices",
         "wrote: %s (%d slice pages + 1 meta page)", outpdf.Data(), kNTimeSlices);
    Info("plot_data_hist_time_slices",
         "lines=%lld, parsed=%lld, phi_ok=%lld, phi_out=%lld, time_ok=%lld, time_out=%lld",
         n_lines, n_parsed, n_phi_ok, n_phi_out, n_time_ok, n_time_out);
}
