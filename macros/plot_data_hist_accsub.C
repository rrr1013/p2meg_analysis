// macros/plot_data_hist_accsub.C
//
// 入力: .dat テキスト。先頭/途中にメタ情報行が混ざっていてもよい。
//      「5列すべてが double として読める行」だけを Event として扱う。
//      列順: Ee Eg t phi_detector_e phi_detector_g
//
// 解析窓: include/p2meg/AnalysisWindow.h の analysis_window を使用。
//        AW : 解析窓4D (Ee, Eg, t, theta_eg) 内
//        TSB: (Ee, Eg, theta_eg) は解析窓内、t は解析窓外（全時間範囲内）
//
// 出力: doc/data_hist_<入力ファイル名(拡張子除く)>_acc_subtracted.pdf（3ページ）
//   1ページ目: メタ情報（TSB幅、スケール、Nacc予測など）
//   2ページ目: 1D（青:差し引き後=折れ線 + 水色の誤差帯、赤:raw=薄い点のみ）
//   3ページ目: 2D（差し引き後のみ）
//
// 実行例（リポジトリ直下から）:
//   root -l -q 'macros/plot_data_hist_accsub.C("data/mockdata/testdata1.dat")'
//

R__ADD_INCLUDE_PATH(./include)

#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>

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

// ---- 固定ビン数（必要ならここだけ調整）----
static constexpr int kNBins_E   = 120; // Ee, Eg
static constexpr int kNBins_t   = 160; // t
static constexpr int kNBins_phi = 120; // phi_detector_e/g
static constexpr int kNBins_th  = 120; // theta_eg

static constexpr int kNBins2D_E   = 120;
static constexpr int kNBins2D_t   = 160;
static constexpr int kNBins2D_th  = 120;
static constexpr int kNBins2D_phi = 120;

// ---- TSB の全時間範囲（基本案）----
// ここは実データの取得レンジに合わせて変更してください。
static constexpr double kTAllMin = -10.0;
static constexpr double kTAllMax =  10.0;

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
    return Form("doc/data_hist_%s_acc_subtracted.pdf", base.Data());
}

static void DrawMetaPage(const char* infile,
                         const char* outpdf,
                         long long n_lines,
                         long long n_parsed,
                         long long n_aw,
                         long long n_tsb,
                         double t_all_min,
                         double t_all_max,
                         double w_win,
                         double w_tsb,
                         double scale_s,
                         double nacc_pred,
                         double nacc_err)
{
    TLatex lat;
    lat.SetNDC(true);

    lat.SetTextAlign(13); // left-top
    lat.SetTextSize(0.040);
    lat.DrawLatex(0.05, 0.92, "p2MEG histograms (ACC-subtracted using TSB)");

    lat.SetTextSize(0.030);
    lat.DrawLatex(0.05, 0.84, Form("input  : %s", infile));
    lat.DrawLatex(0.05, 0.79, Form("output : %s", outpdf));

    lat.DrawLatex(0.05, 0.71, Form("lines read           : %lld", n_lines));
    lat.DrawLatex(0.05, 0.66, Form("parsed (5 doubles)   : %lld", n_parsed));
    lat.DrawLatex(0.05, 0.61, Form("AW events (in window): %lld", n_aw));
    lat.DrawLatex(0.05, 0.56, Form("TSB events (E/theta in win, t in SB): %lld", n_tsb));

    lat.SetTextSize(0.032);
    lat.DrawLatex(0.05, 0.48, "Time definition:");
    lat.SetTextSize(0.030);
    lat.DrawLatex(0.08, 0.43, Form("t_all [ns]   : [%.6g, %.6g]", t_all_min, t_all_max));
    lat.DrawLatex(0.08, 0.38, Form("t_win [ns]   : [%.6g, %.6g]", analysis_window.t_min, analysis_window.t_max));
    lat.DrawLatex(0.08, 0.33, Form("W_win [ns]   : %.6g", w_win));
    lat.DrawLatex(0.08, 0.28, Form("W_TSB [ns]   : %.6g", w_tsb));
    lat.DrawLatex(0.08, 0.23, Form("scale s=W_win/W_TSB : %.6g", scale_s));

    lat.SetTextSize(0.032);
    lat.DrawLatex(0.05, 0.15, "ACC prediction in AW (stat only):");
    lat.SetTextSize(0.030);
    lat.DrawLatex(0.08, 0.10, Form("N_acc_pred = %.6g  +/-  %.6g", nacc_pred, nacc_err));

    lat.SetTextSize(0.026);
    lat.DrawLatex(0.05, 0.05, "Pages: (1) meta  (2) 1D (sub line + error band, raw points)  (3) 2D (sub only)");
}

static void StyleRaw(TH1* h)
{
    if (!h) return;

    // 赤：折れ線（線のみ、エラーバー無し）
    h->SetLineWidth(2);
    h->SetLineStyle(1);
    h->SetLineColorAlpha(kRed+1, 0.55); // 少し濃い赤

    // 点は出さない
    h->SetMarkerSize(0.0);
}

static void StyleSub(TH1* h)
{
    if (!h) return;

    // 青：折れ線（線のみ）
    h->SetLineWidth(2);
    h->SetLineStyle(1);
    h->SetLineColor(kBlue+1);

    // 点は出さない（誤差は帯で出す）
    h->SetMarkerSize(0.0);
}

static TH1D* MakeAccPredT_FromTSB(const TH1D* hRawT, long long n_tsb, double w_tsb)
{
    // 解析窓内で一様 (pedestal) を仮定:
    // rate = n_tsb / w_tsb
    // bin content = rate * bin_width
    if (!hRawT) return nullptr;

    TH1D* h = (TH1D*)hRawT->Clone("hAccPred_t");
    h->Reset("ICES");
    // Clone元がSumw2済みならSumw2は既に確保されている場合がある（警告が出ても問題ない）
    h->Sumw2();

    if (w_tsb <= 0.0 || n_tsb <= 0) return h;

    const double rate = (double)n_tsb / w_tsb; // [counts/ns] in TSB (E/theta in win)
    for (int ib = 1; ib <= h->GetNbinsX(); ++ib) {
        const double w = h->GetXaxis()->GetBinWidth(ib);
        const double mu = rate * w;

        // 統計誤差（TSB総数からレート推定）:
        // Var(mu) = (w/w_tsb)^2 * n_tsb
        const double err = std::sqrt((double)n_tsb) * (w / w_tsb);

        h->SetBinContent(ib, mu);
        h->SetBinError(ib, err);
    }
    return h;
}

static TH2D* MakeAccPred2D_ThT_FromTSB(const TH2D* hRawThT, const TH1D* hTSBTh, double w_tsb)
{
    // ACC: t一様。theta分布はTSBのhTSBThから。
    // content(i_th, j_t) = (M_th(i)/w_tsb) * dt_bin(j)
    if (!hRawThT || !hTSBTh) return nullptr;

    TH2D* h = (TH2D*)hRawThT->Clone("hAccPred_ThT");
    h->Reset("ICES");

    if (w_tsb <= 0.0) return h;

    for (int ith = 1; ith <= h->GetNbinsX(); ++ith) {
        const double M = hTSBTh->GetBinContent(ith);
        const double rate_th = M / w_tsb;
        for (int jt = 1; jt <= h->GetNbinsY(); ++jt) {
            const double dt = h->GetYaxis()->GetBinWidth(jt);
            h->SetBinContent(ith, jt, rate_th * dt);
        }
    }
    return h;
}

static TH2D* MakeAccPred2D_TEe_FromTSB(const TH2D* hRawTEe, const TH1D* hTSBEe, double w_tsb)
{
    // ACC: t一様。Ee分布はTSBのhTSBEeから。
    // content(j_t, k_Ee) = (M_Ee(k)/w_tsb) * dt_bin(j)
    if (!hRawTEe || !hTSBEe) return nullptr;

    TH2D* h = (TH2D*)hRawTEe->Clone("hAccPred_TEe");
    h->Reset("ICES");

    if (w_tsb <= 0.0) return h;

    for (int jt = 1; jt <= h->GetNbinsX(); ++jt) {
        const double dt = h->GetXaxis()->GetBinWidth(jt);
        for (int k = 1; k <= h->GetNbinsY(); ++k) {
            const double M = hTSBEe->GetBinContent(k);
            const double rate_Ee = M / w_tsb;
            h->SetBinContent(jt, k, rate_Ee * dt);
        }
    }
    return h;
}

void plot_data_hist_accsub(const char* infile = "data/data.dat")
{
    const TString outpdf = MakeOutputPdfPath(infile);

    gStyle->SetOptStat(0);

    // 表示範囲は解析窓（AW）に合わせる
    const double Ee_min = analysis_window.Ee_min;
    const double Ee_max = analysis_window.Ee_max;
    const double Eg_min = analysis_window.Eg_min;
    const double Eg_max = analysis_window.Eg_max;
    const double t_min  = analysis_window.t_min;
    const double t_max  = analysis_window.t_max;
    const double th_min = analysis_window.theta_min;
    const double th_max = analysis_window.theta_max;

    const double phi_e_min = detres.phi_e_min;
    const double phi_e_max = detres.phi_e_max;
    const double phi_g_min = detres.phi_g_min;
    const double phi_g_max = detres.phi_g_max;

    // ---- raw (AW, ACC込み) ----
    TH1D* hRawEe   = new TH1D("hRawEe",   "Ee;Ee [MeV];Entries", kNBins_E,  Ee_min, Ee_max);
    TH1D* hRawEg   = new TH1D("hRawEg",   "Eg;Eg [MeV];Entries", kNBins_E,  Eg_min, Eg_max);
    TH1D* hRawt    = new TH1D("hRawt",    "t;t [ns];Entries",    kNBins_t,  t_min,  t_max);
    TH1D* hRawPhiE = new TH1D("hRawPhiE", "phi_{detector,e};phi_{detector,e} [rad];Entries",
                              kNBins_phi, phi_e_min, phi_e_max);
    TH1D* hRawPhiG = new TH1D("hRawPhiG", "phi_{detector,#gamma};phi_{detector,#gamma} [rad];Entries",
                              kNBins_phi, phi_g_min, phi_g_max);
    TH1D* hRawThEg = new TH1D("hRawThEg", "theta_{eg};theta_{eg} [rad];Entries",
                              kNBins_th, th_min, th_max);

    TH2D* hRaw_EeEg = new TH2D("hRaw_EeEg", "(Ee, Eg);Ee [MeV];Eg [MeV]",
                               kNBins2D_E, Ee_min, Ee_max, kNBins2D_E, Eg_min, Eg_max);

    TH2D* hRaw_ThT  = new TH2D("hRaw_ThT", "(theta_{eg}, t);theta_{eg} [rad];t [ns]",
                               kNBins2D_th, th_min, th_max, kNBins2D_t,  t_min,  t_max);

    TH2D* hRaw_TEe  = new TH2D("hRaw_TEe", "(t, Ee);t [ns];Ee [MeV]",
                               kNBins2D_t,  t_min,  t_max, kNBins2D_E, Ee_min, Ee_max);

    TH2D* hRaw_ThEe = new TH2D("hRaw_ThEe", "(theta_{eg}, Ee);theta_{eg} [rad];Ee [MeV]",
                               kNBins2D_th, th_min, th_max, kNBins2D_E, Ee_min, Ee_max);

    TH2D* hRaw_ThEg = new TH2D("hRaw_ThEg", "(theta_{eg}, Eg);theta_{eg} [rad];Eg [MeV]",
                               kNBins2D_th, th_min, th_max, kNBins2D_E, Eg_min, Eg_max);

    TH2D* hRaw_PePg = new TH2D("hRaw_PePg", "(phi_{detector,e}, phi_{detector,#gamma});phi_{detector,e} [rad];phi_{detector,#gamma} [rad]",
                               kNBins2D_phi, phi_e_min, phi_e_max, kNBins2D_phi, phi_g_min, phi_g_max);

    // ---- TSB (E/thetaは窓内, tはSB) ----
    TH1D* hTSBEe   = new TH1D("hTSBEe",   "Ee(TSB);Ee [MeV];Entries", kNBins_E,  Ee_min, Ee_max);
    TH1D* hTSBEg   = new TH1D("hTSBEg",   "Eg(TSB);Eg [MeV];Entries", kNBins_E,  Eg_min, Eg_max);
    TH1D* hTSBPhiE = new TH1D("hTSBPhiE", "phi_{detector,e}(TSB);phi_{detector,e} [rad];Entries",
                              kNBins_phi, phi_e_min, phi_e_max);
    TH1D* hTSBPhiG = new TH1D("hTSBPhiG", "phi_{detector,#gamma}(TSB);phi_{detector,#gamma} [rad];Entries",
                              kNBins_phi, phi_g_min, phi_g_max);
    TH1D* hTSBThEg = new TH1D("hTSBThEg", "theta_{eg}(TSB);theta_{eg} [rad];Entries",
                              kNBins_th, th_min, th_max);

    TH2D* hTSB_EeEg = new TH2D("hTSB_EeEg", "(Ee, Eg) (TSB);Ee [MeV];Eg [MeV]",
                               kNBins2D_E, Ee_min, Ee_max, kNBins2D_E, Eg_min, Eg_max);

    TH2D* hTSB_ThEe = new TH2D("hTSB_ThEe", "(theta_{eg}, Ee) (TSB);theta_{eg} [rad];Ee [MeV]",
                               kNBins2D_th, th_min, th_max, kNBins2D_E, Ee_min, Ee_max);

    TH2D* hTSB_ThEg = new TH2D("hTSB_ThEg", "(theta_{eg}, Eg) (TSB);theta_{eg} [rad];Eg [MeV]",
                               kNBins2D_th, th_min, th_max, kNBins2D_E, Eg_min, Eg_max);

    TH2D* hTSB_PePg = new TH2D("hTSB_PePg", "(phi_{detector,e}, phi_{detector,#gamma}) (TSB);phi_{detector,e} [rad];phi_{detector,#gamma} [rad]",
                               kNBins2D_phi, phi_e_min, phi_e_max, kNBins2D_phi, phi_g_min, phi_g_max);

    // 誤差伝播のため（raw/TSB）
    hRawEe->Sumw2();   hRawEg->Sumw2();   hRawt->Sumw2();    hRawPhiE->Sumw2(); hRawPhiG->Sumw2(); hRawThEg->Sumw2();
    hTSBEe->Sumw2();   hTSBEg->Sumw2();   hTSBPhiE->Sumw2(); hTSBPhiG->Sumw2(); hTSBThEg->Sumw2();

    // 読み込み
    std::ifstream fin(infile);
    if (!fin) {
        Error("plot_data_hist_accsub", "failed to open input file: %s", infile);
        return;
    }

    long long n_lines  = 0;
    long long n_parsed = 0;
    long long n_aw     = 0;
    long long n_tsb    = 0;

    std::string line;
    while (std::getline(fin, line)) {
        ++n_lines;

        double Epos=0.0, Egam=0.0, dt=0.0, phi_pos=0.0, phi_gam=0.0;
        if (!ParseEventLine5Doubles(line, Epos, Egam, dt, phi_pos, phi_gam)) continue;
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
            continue;
        }

        const double phi_e_disc = Detector_PhiGridPoint(idx_e, detres.phi_e_min, detres.phi_e_max,
                                                        Math_GetNPhiE(detres));
        const double phi_g_disc = Detector_PhiGridPoint(idx_g, detres.phi_g_min, detres.phi_g_max,
                                                        Math_GetNPhiG(detres));
        const double theta_eg = std::fabs(phi_e_disc - phi_g_disc);

        // AW (4D window)
        if (AnalysisWindow_In4D(analysis_window, ev.Ee, ev.Eg, ev.t, theta_eg)) {
            ++n_aw;

            hRawEe->Fill(ev.Ee);
            hRawEg->Fill(ev.Eg);
            hRawt->Fill(ev.t);
            hRawPhiE->Fill(ev.phi_detector_e);
            hRawPhiG->Fill(ev.phi_detector_g);
            hRawThEg->Fill(theta_eg);

            hRaw_EeEg->Fill(ev.Ee, ev.Eg);
            hRaw_ThT->Fill(theta_eg, ev.t);
            hRaw_TEe->Fill(ev.t, ev.Ee);
            hRaw_ThEe->Fill(theta_eg, ev.Ee);
            hRaw_ThEg->Fill(theta_eg, ev.Eg);
            hRaw_PePg->Fill(ev.phi_detector_e, ev.phi_detector_g);
        }

        // TSB: (Ee,Eg,theta) in window, time in sideband
        if (AnalysisWindow_In3D(analysis_window, ev.Ee, ev.Eg, theta_eg) &&
            AnalysisWindow_InTimeSideband(analysis_window, ev.t, kTAllMin, kTAllMax)) {
            ++n_tsb;

            hTSBEe->Fill(ev.Ee);
            hTSBEg->Fill(ev.Eg);
            hTSBPhiE->Fill(ev.phi_detector_e);
            hTSBPhiG->Fill(ev.phi_detector_g);
            hTSBThEg->Fill(theta_eg);

            hTSB_EeEg->Fill(ev.Ee, ev.Eg);
            hTSB_ThEe->Fill(theta_eg, ev.Ee);
            hTSB_ThEg->Fill(theta_eg, ev.Eg);
            hTSB_PePg->Fill(ev.phi_detector_e, ev.phi_detector_g);
        }
    }

    if (n_aw == 0) {
        Warning("plot_data_hist_accsub", "no AW events in analysis window. Nothing to plot.");
        return;
    }

    const double w_win = analysis_window.t_max - analysis_window.t_min;
    const double w_tsb = AnalysisWindow_TimeSidebandWidth(analysis_window, kTAllMin, kTAllMax);
    const double s = (w_tsb > 0.0) ? (w_win / w_tsb) : 0.0;

    const double nacc_pred = s * (double)n_tsb;
    const double nacc_err  = (n_tsb > 0 && s > 0.0) ? (s * std::sqrt((double)n_tsb)) : 0.0;

    // ターミナル出力（要求事項）
    Info("plot_data_hist_accsub", "ACC prediction in AW : N_acc_pred = %.6g +/- %.6g", nacc_pred, nacc_err);

    // ---- ACC予測ヒスト（TSBからスケール or 解析的に作る）----
    TH1D* hAccEe   = (TH1D*)hTSBEe->Clone("hAccEe");     hAccEe->Scale(s);
    TH1D* hAccEg   = (TH1D*)hTSBEg->Clone("hAccEg");     hAccEg->Scale(s);
    TH1D* hAccPhiE = (TH1D*)hTSBPhiE->Clone("hAccPhiE"); hAccPhiE->Scale(s);
    TH1D* hAccPhiG = (TH1D*)hTSBPhiG->Clone("hAccPhiG"); hAccPhiG->Scale(s);
    TH1D* hAccThEg = (TH1D*)hTSBThEg->Clone("hAccThEg"); hAccThEg->Scale(s);

    // t は pedestal一様で作る
    TH1D* hAcct = MakeAccPredT_FromTSB(hRawt, n_tsb, w_tsb);

    TH2D* hAcc_EeEg = (TH2D*)hTSB_EeEg->Clone("hAcc_EeEg"); hAcc_EeEg->Scale(s);
    TH2D* hAcc_ThEe = (TH2D*)hTSB_ThEe->Clone("hAcc_ThEe"); hAcc_ThEe->Scale(s);
    TH2D* hAcc_ThEg = (TH2D*)hTSB_ThEg->Clone("hAcc_ThEg"); hAcc_ThEg->Scale(s);
    TH2D* hAcc_PePg = (TH2D*)hTSB_PePg->Clone("hAcc_PePg"); hAcc_PePg->Scale(s);

    // tを含む2Dは外積で作る
    TH2D* hAcc_ThT = MakeAccPred2D_ThT_FromTSB(hRaw_ThT, hTSBThEg, w_tsb);
    TH2D* hAcc_TEe = MakeAccPred2D_TEe_FromTSB(hRaw_TEe, hTSBEe, w_tsb);

    // ---- 差し引き後（sig+rmd）----
    TH1D* hSubEe   = (TH1D*)hRawEe->Clone("hSubEe");       hSubEe->Add(hAccEe,   -1.0);
    TH1D* hSubEg   = (TH1D*)hRawEg->Clone("hSubEg");       hSubEg->Add(hAccEg,   -1.0);
    TH1D* hSubt    = (TH1D*)hRawt->Clone("hSubt");         hSubt->Add(hAcct,     -1.0);
    TH1D* hSubPhiE = (TH1D*)hRawPhiE->Clone("hSubPhiE");   hSubPhiE->Add(hAccPhiE, -1.0);
    TH1D* hSubPhiG = (TH1D*)hRawPhiG->Clone("hSubPhiG");   hSubPhiG->Add(hAccPhiG, -1.0);
    TH1D* hSubThEg = (TH1D*)hRawThEg->Clone("hSubThEg");   hSubThEg->Add(hAccThEg, -1.0);

    TH2D* hSub_EeEg = (TH2D*)hRaw_EeEg->Clone("hSub_EeEg"); hSub_EeEg->Add(hAcc_EeEg, -1.0);
    TH2D* hSub_ThT  = (TH2D*)hRaw_ThT->Clone("hSub_ThT");   if (hAcc_ThT) hSub_ThT->Add(hAcc_ThT, -1.0);
    TH2D* hSub_TEe  = (TH2D*)hRaw_TEe->Clone("hSub_TEe");   if (hAcc_TEe) hSub_TEe->Add(hAcc_TEe, -1.0);
    TH2D* hSub_ThEe = (TH2D*)hRaw_ThEe->Clone("hSub_ThEe"); hSub_ThEe->Add(hAcc_ThEe, -1.0);
    TH2D* hSub_ThEg = (TH2D*)hRaw_ThEg->Clone("hSub_ThEg"); hSub_ThEg->Add(hAcc_ThEg, -1.0);
    TH2D* hSub_PePg = (TH2D*)hRaw_PePg->Clone("hSub_PePg"); hSub_PePg->Add(hAcc_PePg, -1.0);

    // ---- スタイル設定（1D）----
    StyleRaw(hRawEe);   StyleRaw(hRawEg);   StyleRaw(hRawt);    StyleRaw(hRawPhiE); StyleRaw(hRawPhiG); StyleRaw(hRawThEg);
    StyleSub(hSubEe);   StyleSub(hSubEg);   StyleSub(hSubt);    StyleSub(hSubPhiE); StyleSub(hSubPhiG); StyleSub(hSubThEg);

    // ---- ページ1：メタ ----
    TCanvas c0("c0", "meta", 1200, 800);
    c0.cd();
    DrawMetaPage(infile, outpdf.Data(),
                 n_lines, n_parsed, n_aw, n_tsb,
                 kTAllMin, kTAllMax, w_win, w_tsb, s, nacc_pred, nacc_err);

    // ---- ページ2：1D（青:折れ線 + 水色誤差帯、赤:薄い点のみ）----
    auto Draw1DOverlay = [](TH1D* hSub, TH1D* hRaw){
        gPad->SetGrid();
        if (!hSub || !hRaw) return;

        // 1) 誤差帯（水色, 透明）
        const Color_t fill_save  = hSub->GetFillColor();
        const Style_t style_save = hSub->GetFillStyle();

            hSub->SetFillColorAlpha(kCyan+1, 0.25);
        hSub->SetFillStyle(1001);
        hSub->Draw("E2");  // 軸確定

        // 2) 赤の折れ線（先に描く）
        hRaw->Draw("hist same");

        // 3) 青の折れ線（最後に描いて手前に）
        hSub->SetFillColor(fill_save);
        hSub->SetFillStyle(style_save);
        hSub->Draw("hist same");
    };

    TCanvas c1("c1", "1D (sub line + error band, raw points)", 1200, 800);
    c1.Divide(3, 2);
    c1.cd(1); Draw1DOverlay(hSubEe,   hRawEe);
    c1.cd(2); Draw1DOverlay(hSubEg,   hRawEg);
    c1.cd(3); Draw1DOverlay(hSubt,    hRawt);
    c1.cd(4); Draw1DOverlay(hSubPhiE, hRawPhiE);
    c1.cd(5); Draw1DOverlay(hSubPhiG, hRawPhiG);
    c1.cd(6); Draw1DOverlay(hSubThEg, hRawThEg);

    // ---- ページ3：2D（差し引き後のみ）----
    auto Draw2D = [](TH2D* h){
        gPad->SetGrid();
        gPad->SetRightMargin(0.14);
        if (h) h->Draw("colz");
    };

    TCanvas c2("c2", "2D (sub)", 1200, 800);
    c2.Divide(3, 2);
    c2.cd(1); Draw2D(hSub_EeEg);
    c2.cd(2); Draw2D(hSub_ThT);
    c2.cd(3); Draw2D(hSub_TEe);
    c2.cd(4); Draw2D(hSub_ThEe);
    c2.cd(5); Draw2D(hSub_ThEg);
    c2.cd(6); Draw2D(hSub_PePg);

    // ---- PDF（3ページ）----
    c0.Print(Form("%s[", outpdf.Data()));
    c0.Print(outpdf.Data());
    c1.Print(outpdf.Data());
    c2.Print(outpdf.Data());
    c2.Print(Form("%s]", outpdf.Data()));

    Info("plot_data_hist_accsub", "wrote: %s (3 pages)", outpdf.Data());
    Info("plot_data_hist_accsub", "lines=%lld, parsed=%lld, AW=%lld, TSB=%lld, s=%.6g",
         n_lines, n_parsed, n_aw, n_tsb, s);
}
