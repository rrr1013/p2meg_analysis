// macros/plot_data_hist.C
//
// 入力: .dat テキスト。先頭/途中にメタ情報行が混ざっていてもよい。
//      「6列すべてが double として読める行」だけを Event として扱う。
//      列順: Epos Egam dt theta cos_pos cos_gam
//
// 解析窓: include/p2meg/AnalysisWindow.h の analysis_window を使用。
//        解析窓内に入ったイベントのみをヒストに入れる。
//
// 出力: doc/data_hist_<入力ファイル名(拡張子除く)>.pdf（3ページ）
//   1ページ目: メタ情報
//   2ページ目: 1D (Ee, Eg, t, theta, cos_detector_e, cos_detector_g)
//   3ページ目: 2D (Ee,Eg), (theta,t), (t,Ee), (t,Eg), (theta,Ee), (cos_e,cos_g)
//
// 実行例（リポジトリ直下から）:
//   root -l -q 'macros/plot_data_hist.C("data/mockdata/MEGonly_simulation_dataset.dat")'
//

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

#include "../include/p2meg/AnalysisWindow.h"
#include "../include/p2meg/Event.h"

// ---- 固定ビン数（必要ならここだけ調整）----
static constexpr int kNBins_E   = 120; // Ee, Eg
static constexpr int kNBins_t   = 160; // t
static constexpr int kNBins_th  = 120; // theta
static constexpr int kNBins_cos = 120; // cos_detector_e/g

static constexpr int kNBins2D_E   = 120;
static constexpr int kNBins2D_t   = 160;
static constexpr int kNBins2D_th  = 120;
static constexpr int kNBins2D_cos = 120;

static bool ParseEventLine6Doubles(const std::string& line,
                                  double& Epos, double& Egam, double& dt,
                                  double& theta, double& cos_pos, double& cos_gam)
{
    if (line.empty()) return false;

    std::istringstream iss(line);
    if (!(iss >> Epos >> Egam >> dt >> theta >> cos_pos >> cos_gam)) {
        return false; // 6列 double でなければメタ行としてスキップ
    }
    if (!std::isfinite(Epos) || !std::isfinite(Egam) || !std::isfinite(dt) ||
        !std::isfinite(theta) || !std::isfinite(cos_pos) || !std::isfinite(cos_gam)) {
        return false;
    }
    return true;
}

static bool InAnalysisWindow(const Event& ev)
{
    if (ev.Ee    < analysis_window.Ee_min    || ev.Ee    > analysis_window.Ee_max)    return false;
    if (ev.Eg    < analysis_window.Eg_min    || ev.Eg    > analysis_window.Eg_max)    return false;
    if (ev.t     < analysis_window.t_min     || ev.t     > analysis_window.t_max)     return false;
    if (ev.theta < analysis_window.theta_min || ev.theta > analysis_window.theta_max) return false;
    return true;
}

static TString MakeOutputPdfPath(const char* infile)
{
    TString base = gSystem->BaseName(infile); // e.g. xxx.dat
    Ssiz_t dot = base.Last('.');
    if (dot != kNPOS) base.Remove(dot);       // e.g. xxx

    gSystem->mkdir("doc", /*recursive=*/true);
    return Form("doc/data_hist_%s.pdf", base.Data());
}

static void DrawMetaPage(const char* infile,
                         const char* outpdf,
                         long long n_lines,
                         long long n_parsed,
                         long long n_inwin)
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
    lat.DrawLatex(0.05, 0.66, Form("parsed (6 doubles)   : %lld", n_parsed));
    lat.DrawLatex(0.05, 0.61, Form("events in window     : %lld", n_inwin));

    lat.SetTextSize(0.032);
    lat.DrawLatex(0.05, 0.52, "Analysis window:");
    lat.SetTextSize(0.030);
    lat.DrawLatex(0.08, 0.47, Form("Ee [MeV]    : [%.6g, %.6g]", analysis_window.Ee_min, analysis_window.Ee_max));
    lat.DrawLatex(0.08, 0.42, Form("Eg [MeV]    : [%.6g, %.6g]", analysis_window.Eg_min, analysis_window.Eg_max));
    lat.DrawLatex(0.08, 0.37, Form("t  [ns]     : [%.6g, %.6g]", analysis_window.t_min,  analysis_window.t_max));
    lat.DrawLatex(0.08, 0.32, Form("theta [rad] : [%.6g, %.6g]", analysis_window.theta_min, analysis_window.theta_max));

    lat.SetTextSize(0.032);
    lat.DrawLatex(0.05, 0.23, "Fixed bins:");
    lat.SetTextSize(0.030);
    lat.DrawLatex(0.08, 0.18, Form("1D: Ee/Eg=%d  t=%d  theta=%d  cos=%d",
                                  kNBins_E, kNBins_t, kNBins_th, kNBins_cos));
    lat.DrawLatex(0.08, 0.13, Form("2D: E=%d  t=%d  theta=%d  cos=%d",
                                  kNBins2D_E, kNBins2D_t, kNBins2D_th, kNBins2D_cos));

    lat.SetTextSize(0.026);
    lat.DrawLatex(0.05, 0.06, "Pages: (1) meta  (2) 1D  (3) 2D");
}

void plot_data_hist(const char* infile = "data/data.dat")
{
    const TString outpdf = MakeOutputPdfPath(infile);

    gStyle->SetOptStat(0);

    // 表示範囲＝解析窓
    const double Ee_min = analysis_window.Ee_min;
    const double Ee_max = analysis_window.Ee_max;
    const double Eg_min = analysis_window.Eg_min;
    const double Eg_max = analysis_window.Eg_max;
    const double t_min  = analysis_window.t_min;
    const double t_max  = analysis_window.t_max;
    const double th_min = analysis_window.theta_min;
    const double th_max = analysis_window.theta_max;

    // ---- 1D ----
    TH1D* hEe   = new TH1D("hEe",   "Ee;Ee [MeV];Entries",                 kNBins_E,  Ee_min, Ee_max);
    TH1D* hEg   = new TH1D("hEg",   "Eg;Eg [MeV];Entries",                 kNBins_E,  Eg_min, Eg_max);
    TH1D* ht    = new TH1D("ht",    "t;t [ns];Entries",                    kNBins_t,  t_min,  t_max);
    TH1D* hTh   = new TH1D("hTh",   "theta;theta [rad];Entries",           kNBins_th, th_min, th_max);
    TH1D* hCose = new TH1D("hCose", "cos_{detector,e};cos_{detector,e};Entries",
                           kNBins_cos, -1.0, 1.0);
    TH1D* hCosg = new TH1D("hCosg", "cos_{detector,#gamma};cos_{detector,#gamma};Entries",
                           kNBins_cos, -1.0, 1.0);

    // ---- 2D ----
    TH2D* h_EeEg = new TH2D("h_EeEg", "(Ee, Eg);Ee [MeV];Eg [MeV]",
                            kNBins2D_E, Ee_min, Ee_max, kNBins2D_E, Eg_min, Eg_max);

    TH2D* h_ThT  = new TH2D("h_ThT", "(theta, t);theta [rad];t [ns]",
                            kNBins2D_th, th_min, th_max, kNBins2D_t,  t_min,  t_max);

    TH2D* h_TEe  = new TH2D("h_TEe", "(t, Ee);t [ns];Ee [MeV]",
                            kNBins2D_t,  t_min,  t_max, kNBins2D_E, Ee_min, Ee_max);

    TH2D* h_TEg  = new TH2D("h_TEg", "(t, Eg);t [ns];Eg [MeV]",
                            kNBins2D_t,  t_min,  t_max, kNBins2D_E, Eg_min, Eg_max);

    TH2D* h_ThEe = new TH2D("h_ThEe", "(theta, Ee);theta [rad];Ee [MeV]",
                            kNBins2D_th, th_min, th_max, kNBins2D_E, Ee_min, Ee_max);

    TH2D* h_CeCg = new TH2D("h_CeCg", "(cos_detector_e, cos_detector_g);cos_{detector,e};cos_{detector,#gamma}",
                            kNBins2D_cos, -1.0, 1.0, kNBins2D_cos, -1.0, 1.0);

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
    long long n_cos_out = 0;

    std::string line;
    while (std::getline(fin, line)) {
        ++n_lines;

        double Epos=0.0, Egam=0.0, dt=0.0, theta=0.0, cos_pos=0.0, cos_gam=0.0;
        if (!ParseEventLine6Doubles(line, Epos, Egam, dt, theta, cos_pos, cos_gam)) {
            continue;
        }
        ++n_parsed;

        Event ev;
        ev.Ee = Epos;
        ev.Eg = Egam;
        ev.t  = dt;
        ev.theta = theta;
        ev.cos_detector_e = cos_pos;
        ev.cos_detector_g = cos_gam;

        if (!InAnalysisWindow(ev)) {
            ++n_outwin;
            continue;
        }
        ++n_inwin;

        if (ev.cos_detector_e < -1.0 || ev.cos_detector_e > 1.0 ||
            ev.cos_detector_g < -1.0 || ev.cos_detector_g > 1.0) {
            ++n_cos_out;
        }

        // 1D
        hEe->Fill(ev.Ee);
        hEg->Fill(ev.Eg);
        ht->Fill(ev.t);
        hTh->Fill(ev.theta);
        hCose->Fill(ev.cos_detector_e);
        hCosg->Fill(ev.cos_detector_g);

        // 2D
        h_EeEg->Fill(ev.Ee, ev.Eg);
        h_ThT->Fill(ev.theta, ev.t);
        h_TEe->Fill(ev.t, ev.Ee);
        h_TEg->Fill(ev.t, ev.Eg);
        h_ThEe->Fill(ev.theta, ev.Ee);
        h_CeCg->Fill(ev.cos_detector_e, ev.cos_detector_g);
    }

    if (n_inwin == 0) {
        Warning("plot_data_hist", "no events in analysis window. Nothing to plot.");
        return;
    }

    // ---- ページ1：メタ ----
    TCanvas c0("c0", "meta", 1200, 800);
    c0.cd();
    DrawMetaPage(infile, outpdf.Data(), n_lines, n_parsed, n_inwin);

    // ---- ページ2：1D ----
    TCanvas c1("c1", "1D sanity", 1200, 800);
    c1.Divide(3, 2);

    c1.cd(1); gPad->SetGrid(); hEe->SetLineWidth(2); hEe->Draw("hist");
    c1.cd(2); gPad->SetGrid(); hEg->SetLineWidth(2); hEg->Draw("hist");
    c1.cd(3); gPad->SetGrid(); ht->SetLineWidth(2);  ht->Draw("hist");
    c1.cd(4); gPad->SetGrid(); hTh->SetLineWidth(2); hTh->Draw("hist");
    c1.cd(5); gPad->SetGrid(); hCose->SetLineWidth(2); hCose->Draw("hist");
    c1.cd(6); gPad->SetGrid(); hCosg->SetLineWidth(2); hCosg->Draw("hist");

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
    c2.cd(4); Draw2D(h_TEg);
    c2.cd(5); Draw2D(h_ThEe);
    c2.cd(6); Draw2D(h_CeCg);

    // ---- PDF（3ページ）----
    c0.Print(Form("%s[", outpdf.Data()));
    c0.Print(outpdf.Data());
    c1.Print(outpdf.Data());
    c2.Print(outpdf.Data());
    c2.Print(Form("%s]", outpdf.Data()));

    Info("plot_data_hist", "wrote: %s (3 pages)", outpdf.Data());
    Info("plot_data_hist", "lines=%lld, parsed=%lld, inwin=%lld, outwin=%lld, cos_out=%lld",
         n_lines, n_parsed, n_inwin, n_outwin, n_cos_out);
}
