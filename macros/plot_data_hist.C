// macros/plot_data_hist.C
//
// 入力: .dat テキスト。先頭/途中にメタ情報行が混ざっていてもよい。
//      「5列すべてが double として読める行」だけを Event として扱う。
//      列順: Ee Eg t phi_detector_e phi_detector_g
//
// 旧フォーマット(6列: Epos Egam dt theta cos_pos cos_gam)からの移行手順:
//   1) Ee, Eg, t は 1-3 列をそのまま使う。
//   2) cos_pos, cos_gam から角度を作る場合は
//      phi_detector_e = acos(cos_pos), phi_detector_g = acos(cos_gam)
//      として 4-5 列に書き出す（符号情報が失われるので注意）。
//   3) theta_eg は新フォーマット側で |phi_e - phi_g| として計算する。
//
// 解析窓: include/p2meg/AnalysisWindow.h の analysis_window を使用。
//        解析窓内に入ったイベントのみをヒストに入れる。
//
// 出力: doc/data_hist_<入力ファイル名(拡張子除く)>.pdf（3ページ）
//   1ページ目: メタ情報
//   2ページ目: 1D (Ee, Eg, t, phi_detector_e, phi_detector_g, theta_eg)
//   3ページ目: 2D (Ee,Eg), (theta_eg,t), (t,Ee), (theta_eg,Ee), (theta_eg,Eg), (phi_e,phi_g)
//
// 実行例（リポジトリ直下から）:
//   root -l -q 'macros/plot_data_hist.C("data/mockdata/testdata1.dat")'
//

R__ADD_INCLUDE_PATH(./include)

#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TError.h"
#include "TString.h"
#include "TLatex.h"

#include "p2meg/AnalysisWindow.h"
#include "p2meg/AnalysisWindowUtils.h"
#include "p2meg/Event.h"
#include "p2meg/DetectorResolution.h"
#include "p2meg/MathUtils.h"

// ---- 解析窓カットを外して素の分布を見る場合は、下のコメントアウトを外す ----
// #define P2MEG_PLOT_ALLDATA

// ---- 固定ビン数（必要ならここだけ調整）----
static constexpr int kNBins_E   = 120; // Ee, Eg
static constexpr int kNBins_t   = 160; // t
static constexpr int kNBins_phi = 120; // phi_detector_e/g
static constexpr int kNBins_th  = 120; // theta_eg

static constexpr int kNBins2D_E   = 120;
static constexpr int kNBins2D_t   = 160;
static constexpr int kNBins2D_th  = 120;
static constexpr int kNBins2D_phi = 120;

static bool ParseEventLine5Doubles(const std::string& line,
                                  double& Epos, double& Egam, double& dt,
                                  double& phi_pos, double& phi_gam)
{
    if (line.empty()) return false;

    std::istringstream iss(line);
    if (!(iss >> Epos >> Egam >> dt >> phi_pos >> phi_gam)) {
        return false; // 5列 double でなければメタ行としてスキップ
    }
    if (!std::isfinite(Epos) || !std::isfinite(Egam) || !std::isfinite(dt) ||
        !std::isfinite(phi_pos) || !std::isfinite(phi_gam)) {
        return false;
    }
    return true;
}


static TString MakeOutputPdfPath(const char* infile)
{
    TString base = gSystem->BaseName(infile); // e.g. xxx.dat
    Ssiz_t dot = base.Last('.');
    if (dot != kNPOS) base.Remove(dot);       // e.g. xxx

    gSystem->mkdir("doc", /*recursive=*/true);
#ifdef P2MEG_PLOT_ALLDATA
    return Form("doc/data_hist_%s_alldata.pdf", base.Data());
#else
    return Form("doc/data_hist_%s.pdf", base.Data());
#endif
}

static void DrawMetaPage(const char* infile,
                         const char* outpdf,
                         long long n_lines,
                         long long n_parsed,
                         long long n_inwin,
                         double th_plot_min,
                         double th_plot_max)
{
    TLatex lat;
    lat.SetNDC(true);

    lat.SetTextAlign(13); // left-top
    lat.SetTextSize(0.040);
    lat.DrawLatex(0.05, 0.92, "p2MEG histograms");

    lat.SetTextSize(0.030);
    lat.DrawLatex(0.05, 0.84, Form("input  : %s", infile));
    lat.DrawLatex(0.05, 0.79, Form("output : %s", outpdf));

    lat.DrawLatex(0.05, 0.71, Form("lines read           : %lld", n_lines));
    lat.DrawLatex(0.05, 0.66, Form("parsed (5 doubles)   : %lld", n_parsed));
#ifdef P2MEG_PLOT_ALLDATA
    lat.DrawLatex(0.05, 0.61, Form("events plotted       : %lld", n_inwin));
#else
    lat.DrawLatex(0.05, 0.61, Form("events in window     : %lld", n_inwin));
#endif

    lat.SetTextSize(0.032);
#ifdef P2MEG_PLOT_ALLDATA
    lat.DrawLatex(0.05, 0.52, "Plot range (alldata):");
    lat.SetTextSize(0.030);
    lat.DrawLatex(0.08, 0.47, "Ee [MeV]    : [0, 70]");
    lat.DrawLatex(0.08, 0.42, "Eg [MeV]    : [0, 70]");
    lat.DrawLatex(0.08, 0.37, "t  [ns]     : [-10, 10]");
    lat.DrawLatex(0.08, 0.32, Form("theta_eg [rad] : [%.6g, %.6g]", th_plot_min, th_plot_max));
#else
    lat.DrawLatex(0.05, 0.52, "Analysis window:");
    lat.SetTextSize(0.030);
    lat.DrawLatex(0.08, 0.47, Form("Ee [MeV]    : [%.6g, %.6g]", analysis_window.Ee_min, analysis_window.Ee_max));
    lat.DrawLatex(0.08, 0.42, Form("Eg [MeV]    : [%.6g, %.6g]", analysis_window.Eg_min, analysis_window.Eg_max));
    lat.DrawLatex(0.08, 0.37, Form("t  [ns]     : [%.6g, %.6g]", analysis_window.t_min,  analysis_window.t_max));
    lat.DrawLatex(0.08, 0.32, Form("theta_eg [rad] : [%.6g, %.6g]", analysis_window.theta_min, analysis_window.theta_max));
#endif

    lat.SetTextSize(0.032);
    lat.DrawLatex(0.05, 0.23, "Fixed bins:");
    lat.SetTextSize(0.030);
    lat.DrawLatex(0.08, 0.18, Form("1D: Ee/Eg=%d  t=%d  phi=%d  theta_eg=%d",
                                  kNBins_E, kNBins_t, kNBins_phi, kNBins_th));
    lat.DrawLatex(0.08, 0.13, Form("2D: E=%d  t=%d  theta_eg=%d  phi=%d",
                                  kNBins2D_E, kNBins2D_t, kNBins2D_th, kNBins2D_phi));

    lat.SetTextSize(0.026);
    lat.DrawLatex(0.05, 0.06, "Pages: (1) meta  (2) 1D  (3) 2D");
}

void plot_data_hist(const char* infile = "data/data.dat")
{
    const TString outpdf = MakeOutputPdfPath(infile);

    gStyle->SetOptStat(0);

    // 表示範囲
    const double pi_val = 3.14159265358979323846;
#ifdef P2MEG_PLOT_ALLDATA
    const double Ee_min = 0.0;
    const double Ee_max = 70.0;
    const double Eg_min = 0.0;
    const double Eg_max = 70.0;
    const double t_min  = -10.0;
    const double t_max  = 10.0;
#else
    const double Ee_min = analysis_window.Ee_min;
    const double Ee_max = analysis_window.Ee_max;
    const double Eg_min = analysis_window.Eg_min;
    const double Eg_max = analysis_window.Eg_max;
    const double t_min  = analysis_window.t_min;
    const double t_max  = analysis_window.t_max;
    const double th_win_min = analysis_window.theta_min;
    const double th_win_max = analysis_window.theta_max;
#endif
    double th_plot_min = 0.0;
    double th_plot_max = pi_val;
#ifndef P2MEG_PLOT_ALLDATA
    th_plot_min = th_win_min;
    th_plot_max = th_win_max;
#endif
    const double phi_e_min = detres.phi_e_min;
    const double phi_e_max = detres.phi_e_max;
    const double phi_g_min = detres.phi_g_min;
    const double phi_g_max = detres.phi_g_max;
    const double phi_e_max_plot = Detector_PhiAxisMaxInclusive(phi_e_max);
    const double phi_g_max_plot = Detector_PhiAxisMaxInclusive(phi_g_max);

    double th_det_min = 0.0;
    double th_det_max = 0.0;
    if (Detector_ThetaRangeFromAllowedPhi(detres, th_det_min, th_det_max)) {
#ifdef P2MEG_PLOT_ALLDATA
        th_plot_min = th_det_min;
        th_plot_max = th_det_max;
#else
        if (th_det_min > th_plot_min) th_plot_min = th_det_min;
        if (th_det_max < th_plot_max) th_plot_max = th_det_max;
        if (!(th_plot_min < th_plot_max)) {
            th_plot_min = th_win_min;
            th_plot_max = th_win_max;
        }
#endif
    }
    const double th_plot_max_axis = Math_AxisMaxInclusive(th_plot_max);

    // ---- 1D ----
    TH1D* hEe   = new TH1D("hEe",   "Ee;Ee [MeV];Entries",                 kNBins_E,  Ee_min, Ee_max);
    TH1D* hEg   = new TH1D("hEg",   "Eg;Eg [MeV];Entries",                 kNBins_E,  Eg_min, Eg_max);
    TH1D* ht    = new TH1D("ht",    "t;t [ns];Entries",                    kNBins_t,  t_min,  t_max);
    TH1D* hPhiE = new TH1D("hPhiE", "phi_{detector,e};phi_{detector,e} [rad];Entries",
                           kNBins_phi, phi_e_min, phi_e_max_plot);
    TH1D* hPhiG = new TH1D("hPhiG", "phi_{detector,#gamma};phi_{detector,#gamma} [rad];Entries",
                           kNBins_phi, phi_g_min, phi_g_max_plot);
    TH1D* hThEg = new TH1D("hThEg", "theta_{eg};theta_{eg} [rad];Entries",
                           kNBins_th, th_plot_min, th_plot_max_axis);

    // ---- 2D ----
    TH2D* h_EeEg = new TH2D("h_EeEg", "(Ee, Eg);Ee [MeV];Eg [MeV]",
                            kNBins2D_E, Ee_min, Ee_max, kNBins2D_E, Eg_min, Eg_max);

    TH2D* h_ThT  = new TH2D("h_ThT", "(theta_{eg}, t);theta_{eg} [rad];t [ns]",
                            kNBins2D_th, th_plot_min, th_plot_max_axis, kNBins2D_t,  t_min,  t_max);

    TH2D* h_TEe  = new TH2D("h_TEe", "(t, Ee);t [ns];Ee [MeV]",
                            kNBins2D_t,  t_min,  t_max, kNBins2D_E, Ee_min, Ee_max);

    TH2D* h_ThEe = new TH2D("h_ThEe", "(theta_{eg}, Ee);theta_{eg} [rad];Ee [MeV]",
                            kNBins2D_th, th_plot_min, th_plot_max_axis, kNBins2D_E, Ee_min, Ee_max);

    TH2D* h_ThEg = new TH2D("h_ThEg", "(theta_{eg}, Eg);theta_{eg} [rad];Eg [MeV]",
                            kNBins2D_th, th_plot_min, th_plot_max_axis, kNBins2D_E, Eg_min, Eg_max);

    TH2D* h_PePg = new TH2D("h_PePg", "(phi_{detector,e}, phi_{detector,#gamma});phi_{detector,e} [rad];phi_{detector,#gamma} [rad]",
                            kNBins2D_phi, phi_e_min, phi_e_max_plot, kNBins2D_phi, phi_g_min, phi_g_max_plot);

    // 読み込み
    std::ifstream fin(infile);
    if (!fin) {
        Error("plot_data_hist", "failed to open input file: %s", infile);
        return;
    }

    long long n_lines  = 0;
    long long n_parsed = 0;
    long long n_inwin  = 0;

    // 参考用に数えるが PDF には出さない
    long long n_outwin = 0;
    long long n_phi_out = 0;

    std::string line;
    while (std::getline(fin, line)) {
        ++n_lines;

        double Epos=0.0, Egam=0.0, dt=0.0, phi_pos=0.0, phi_gam=0.0;
        if (!ParseEventLine5Doubles(line, Epos, Egam, dt, phi_pos, phi_gam)) {
            continue;
        }
        ++n_parsed;

        Event ev;
        ev.Ee = Epos;
        ev.Eg = Egam;
        ev.t  = dt;
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

        const double phi_e_disc = Detector_PhiGridPoint(idx_e, detres.phi_e_min, detres.phi_e_max,
                                                        Math_GetNPhiE(detres));
        const double phi_g_disc = Detector_PhiGridPoint(idx_g, detres.phi_g_min, detres.phi_g_max,
                                                        Math_GetNPhiG(detres));
        const double phi_e_plot = phi_e_disc;
        const double phi_g_plot = phi_g_disc;
        double theta_eg = std::fabs(phi_e_plot - phi_g_plot);

#ifndef P2MEG_PLOT_ALLDATA
        if (!AnalysisWindow_In4D(analysis_window, ev.Ee, ev.Eg, ev.t, theta_eg)) {
            ++n_outwin;
            continue;
        }
        ++n_inwin;
#else
        ++n_inwin;
#endif

        // 1D
        hEe->Fill(ev.Ee);
        hEg->Fill(ev.Eg);
        ht->Fill(ev.t);
        hPhiE->Fill(phi_e_plot);
        hPhiG->Fill(phi_g_plot);
        hThEg->Fill(theta_eg);

        // 2D
        h_EeEg->Fill(ev.Ee, ev.Eg);
        h_ThT->Fill(theta_eg, ev.t);
        h_TEe->Fill(ev.t, ev.Ee);
        h_ThEe->Fill(theta_eg, ev.Ee);
        h_ThEg->Fill(theta_eg, ev.Eg);
        h_PePg->Fill(phi_e_plot, phi_g_plot);
    }

    if (n_inwin == 0) {
        Warning("plot_data_hist", "no events in analysis window. Nothing to plot.");
        return;
    }

    // ---- ページ1：メタ ----
    TCanvas c0("c0", "meta", 1200, 800);
    c0.cd();
    DrawMetaPage(infile, outpdf.Data(), n_lines, n_parsed, n_inwin, th_plot_min, th_plot_max);

    // ---- ページ2：1D ----
    TCanvas c1("c1", "1D sanity", 1200, 800);
    c1.Divide(3, 2);

    c1.cd(1); gPad->SetGrid(); hEe->SetLineWidth(2); hEe->Draw("hist");
    c1.cd(2); gPad->SetGrid(); hEg->SetLineWidth(2); hEg->Draw("hist");
    c1.cd(3); gPad->SetGrid(); ht->SetLineWidth(2); ht->SetMinimum(0);  ht->Draw("hist");
    c1.cd(4); gPad->SetGrid(); hPhiE->SetLineWidth(2); hPhiE->Draw("hist");
    c1.cd(5); gPad->SetGrid(); hPhiG->SetLineWidth(2); hPhiG->Draw("hist");
    c1.cd(6); gPad->SetGrid(); hThEg->SetLineWidth(2); hThEg->Draw("hist");

    // ---- ページ3：2D ----
    auto Draw2D = [](TH2D* h){
        gPad->SetGrid();
        gPad->SetRightMargin(0.14);
        h->Draw("colz");
    };

    TCanvas c2("c2", "2D sanity", 1200, 800);
    c2.Divide(3, 2);
    c2.cd(1); Draw2D(h_EeEg);
    c2.cd(2); Draw2D(h_ThT);
    c2.cd(3); Draw2D(h_TEe);
    c2.cd(4); Draw2D(h_ThEe);
    c2.cd(5); Draw2D(h_ThEg);
    c2.cd(6); Draw2D(h_PePg);

    // ---- PDF（3ページ）----
    c0.Print(Form("%s[", outpdf.Data()));
    c0.Print(outpdf.Data());
    c1.Print(outpdf.Data());
    c2.Print(outpdf.Data());
    c2.Print(Form("%s]", outpdf.Data()));

    Info("plot_data_hist", "wrote: %s (3 pages)", outpdf.Data());
    Info("plot_data_hist", "lines=%lld, parsed=%lld, inwin=%lld, outwin=%lld, phi_out=%lld",
         n_lines, n_parsed, n_inwin, n_outwin, n_phi_out);
}
