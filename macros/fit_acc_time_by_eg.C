// macros/fit_acc_time_by_eg.C
//
// 目的:
//   - 最終5Dデータから、ACC 優勢とみなす全データの t 分布を
//     Eg bin ごとに調べる。
//   - |t| < blind_core_ns は blind にして、対称な smooth 関数
//     （単一ガウス + 定数 pedestal）で broad な時間構造を記述する。
//   - まずは「この関数形が見た目に無理なく broad ACC shape を再現できるか」
//     を確認するための PDF を出す。
//
// 入力:
//   - 5列テキスト: Ee Eg t phi_detector_e phi_detector_g
//
// 選別:
//   - phi は DetectorResolution の許可ペアに丸めて判定
//   - Ee は解析窓内
//   - theta_eg = |phi_e - phi_g| は解析窓内
//   - Eg は bin ごとに選別（解析窓内に限らない）
//   - t は [-500, 500] ns に制限
//
// 注意:
//   - 時間形状 fit に使うヒストグラムには、
//     4D 解析窓 (Ee, Eg, theta_eg, t) の内側に入る事象を入れない。
//   - ただし R_data = N_AW/N_TSB の生カウントは、従来通り
//     実データの AW/TSB 計数をそのまま使う。
//
// フィット:
//   - p(t) = A * exp(-t^2/(2*sigma^2)) + C
//   - 偶関数（中心 0 固定）を仮定
//   - blind 領域 |t| < blind_core_ns は TF1::RejectPoint() で除外
//   - これは RMD/signal の狭い prompt 成分を見ないための解析上の blind であり、
//     物理カットではない
//
// 出力:
//   - doc/finalanalysis/fit_acc_time_by_eg_<input>.pdf
//
// 実行例:
//   root -l -q 'macros/fit_acc_time_by_eg.C("data/finaldata/step4f_run1to20_eg_doubleonly_allpatterns_500ns.txt")'

R__ADD_INCLUDE_PATH(./include)

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdlib>

#include "TBox.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFitResultPtr.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TError.h"

#include "p2meg/AnalysisWindow.h"
#include "p2meg/AccTimeFit.h"
#include "p2meg/DetectorResolution.h"
#include "p2meg/Event.h"
#include "p2meg/MathUtils.h"

static constexpr double kTAllMin = -500.0;   // [ns]
static constexpr double kTAllMax =  500.0;   // [ns]
static constexpr int    kNBinsT  = 100;
static constexpr double kBlindCoreNs = 20.0; // [ns]

// Eg の依存性を見るため、低エネルギー側を細かく切り、
// 高エネルギー側は統計確保のため 30-80 MeV でまとめる。
static const double kEgEdges[] = {0.0, 10.0, 20.0, 30.0, 80.0};
static constexpr int kNEgBins = static_cast<int>(sizeof(kEgEdges) / sizeof(kEgEdges[0])) - 1;

static double gBlindCoreNs = kBlindCoreNs;

static double ReadSigmaMinNsFromEnv()
{
    const char* s = std::getenv("P2MEG_ACC_SIGMA_MIN_NS");
    if (!s) return 20.0;
    char* endptr = nullptr;
    const double v = std::strtod(s, &endptr);
    if (endptr == s || !std::isfinite(v) || !(v > 0.0)) return 20.0;
    return v;
}

static double ReadBlindCoreNsFromEnv()
{
    const char* s = std::getenv("P2MEG_ACC_BLIND_CORE_NS");
    if (!s) return kBlindCoreNs;
    char* endptr = nullptr;
    const double v = std::strtod(s, &endptr);
    if (endptr == s || !std::isfinite(v) || !(v > 0.0)) return kBlindCoreNs;
    return v;
}

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

static TString MakeOutputPdfPath(const char* infile)
{
    TString base = gSystem->BaseName(infile);
    const Ssiz_t dot = base.Last('.');
    if (dot != kNPOS) base.Remove(dot);
    gSystem->mkdir("doc/finalanalysis", /*recursive=*/true);
    return Form("doc/finalanalysis/fit_acc_time_by_eg_%s.pdf", base.Data());
}

static int FindEgBin(double Eg)
{
    if (!std::isfinite(Eg)) return -1;
    for (int i = 0; i < kNEgBins; ++i) {
        const double lo = kEgEdges[i];
        const double hi = kEgEdges[i + 1];
        if ((Eg >= lo && Eg < hi) || (i == kNEgBins - 1 && Eg >= lo && Eg <= hi)) {
            return i;
        }
    }
    return -1;
}

static bool IsInsideFullAnalysisWindow(double Ee, double Eg, double theta_eg, double t)
{
    if (!std::isfinite(Ee) || !std::isfinite(Eg) ||
        !std::isfinite(theta_eg) || !std::isfinite(t)) {
        return false;
    }
    if (Ee < analysis_window.Ee_min || Ee > analysis_window.Ee_max) return false;
    if (Eg < analysis_window.Eg_min || Eg > analysis_window.Eg_max) return false;
    if (theta_eg < analysis_window.theta_min || theta_eg > analysis_window.theta_max) return false;
    if (t < analysis_window.t_min || t > analysis_window.t_max) return false;
    return true;
}

static double HistIntegralRange(const TH1D& h, double x_min, double x_max)
{
    if (!(x_max > x_min)) return 0.0;
    double sum = 0.0;
    for (int ib = 1; ib <= h.GetNbinsX(); ++ib) {
        const double lo = h.GetXaxis()->GetBinLowEdge(ib);
        const double hi = lo + h.GetXaxis()->GetBinWidth(ib);
        const double ov = std::min(hi, x_max) - std::max(lo, x_min);
        if (!(ov > 0.0)) continue;
        const double dens = h.GetBinContent(ib);
        if (std::isfinite(dens) && dens >= 0.0) sum += dens * ov;
    }
    return sum;
}

static void ConvertCountsToDensity(TH1D& h)
{
    for (int ib = 1; ib <= h.GetNbinsX(); ++ib) {
        const double w = h.GetXaxis()->GetBinWidth(ib);
        if (!(w > 0.0) || !std::isfinite(w)) {
            h.SetBinContent(ib, 0.0);
            h.SetBinError(ib, 0.0);
            continue;
        }
        h.SetBinContent(ib, h.GetBinContent(ib) / w);
        h.SetBinError(ib, h.GetBinError(ib) / w);
    }
}

static void StyleHist(TH1D& h)
{
    h.SetLineColor(kBlack);
    h.SetLineWidth(2);
    h.SetMarkerStyle(20);
    h.SetMarkerSize(0.7);
}

static void DrawMetaPage(const char* infile,
                         const char* outpdf,
                         long long n_lines,
                         long long n_parsed,
                         long long n_selected,
                         long long n_fit_used,
                         const long long n_eg_bin[kNEgBins])
{
    TLatex lat;
    lat.SetNDC(true);
    lat.SetTextAlign(13);

    lat.SetTextSize(0.040);
    lat.DrawLatex(0.05, 0.93, "ACC time-shape fit by Eg bin");

    lat.SetTextSize(0.030);
    lat.DrawLatex(0.05, 0.85, Form("input  : %s", infile));
    lat.DrawLatex(0.05, 0.80, Form("output : %s", outpdf));

    lat.DrawLatex(0.05, 0.70, Form("lines read           : %lld", n_lines));
    lat.DrawLatex(0.05, 0.65, Form("parsed (5 doubles)   : %lld", n_parsed));
    lat.DrawLatex(0.05, 0.60, Form("selected (Ee/theta, |t|<500) : %lld", n_selected));
    lat.DrawLatex(0.05, 0.55, Form("used in fit hist (4D AW excluded): %lld", n_fit_used));

    lat.DrawLatex(0.05, 0.50, "Selection:");
    lat.DrawLatex(0.08, 0.45, Form("Ee in [%.1f, %.1f] MeV", analysis_window.Ee_min, analysis_window.Ee_max));
    lat.DrawLatex(0.08, 0.40, Form("theta_eg in [%.3f, %.3f] rad", analysis_window.theta_min, analysis_window.theta_max));
    lat.DrawLatex(0.08, 0.35, Form("t in [%.0f, %.0f] ns", kTAllMin, kTAllMax));
    lat.DrawLatex(0.08, 0.30, Form("blind core: |t| < %.1f ns", gBlindCoreNs));

    lat.DrawLatex(0.05, 0.22, "Eg bins:");
    for (int i = 0; i < kNEgBins; ++i) {
        const double x = (i < 3) ? 0.08 : 0.45;
        const double y = (i < 3) ? (0.17 - 0.05 * i) : (0.17 - 0.05 * (i - 3));
        lat.DrawLatex(x, y,
                      Form("[%.0f, %.0f]%s MeV : %lld entries",
                           kEgEdges[i], kEgEdges[i + 1],
                           (i == kNEgBins - 1) ? "" : "",
                           n_eg_bin[i]));
    }
}

static TString MakeSummaryTxtPath(const char* infile)
{
    TString base = gSystem->BaseName(infile);
    const Ssiz_t dot = base.Last('.');
    if (dot != kNPOS) base.Remove(dot);
    gSystem->mkdir("doc/finalanalysis", /*recursive=*/true);
    return Form("doc/finalanalysis/fit_acc_time_by_eg_%s.txt", base.Data());
}

void fit_acc_time_by_eg(
    const char* infile = "data/finaldata/step4f_run1to20_eg_doubleonly_allpatterns_500ns.txt")
{
    gBlindCoreNs = ReadBlindCoreNsFromEnv();
    const double sigma_min_ns = ReadSigmaMinNsFromEnv();
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    TH1D* hT[kNEgBins];
        long long n_eg_bin[kNEgBins];
        long long n_fit_eg_bin[kNEgBins];
    long long n_aw_bin[kNEgBins];
    long long n_tsb_bin[kNEgBins];
    for (int i = 0; i < kNEgBins; ++i) {
        hT[i] = new TH1D(Form("hT_egbin_%d", i),
                         Form("t, Eg in [%.0f, %.0f];t [ns];density [counts/ns]",
                              kEgEdges[i], kEgEdges[i + 1]),
                         kNBinsT, kTAllMin, kTAllMax);
        hT[i]->Sumw2();
        n_eg_bin[i] = 0;
        n_fit_eg_bin[i] = 0;
        n_aw_bin[i] = 0;
        n_tsb_bin[i] = 0;
    }

    std::ifstream fin(infile);
    if (!fin) {
        Error("fit_acc_time_by_eg", "failed to open input file: %s", infile);
        return;
    }

    long long n_lines = 0;
    long long n_parsed = 0;
    long long n_selected = 0;
    long long n_fit_used = 0;

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

        if (!std::isfinite(t) || t < kTAllMin || t > kTAllMax) continue;

        int idx_e = -1;
        int idx_g = -1;
        if (!Detector_IsAllowedPhiPairValue(phi_e, phi_g, detres, idx_e, idx_g)) continue;

        const double phi_e_disc =
            Detector_PhiGridPoint(idx_e, detres.phi_e_min, detres.phi_e_max, detres.N_phi_e);
        const double phi_g_disc =
            Detector_PhiGridPoint(idx_g, detres.phi_g_min, detres.phi_g_max, detres.N_phi_g);
        const double theta_eg = std::fabs(phi_e_disc - phi_g_disc);

        if (Ee < analysis_window.Ee_min || Ee > analysis_window.Ee_max) continue;
        if (theta_eg < analysis_window.theta_min || theta_eg > analysis_window.theta_max) continue;

        const int ib = FindEgBin(Eg);
        if (ib < 0) continue;

        ++n_eg_bin[ib];
        if (t >= analysis_window.t_min && t <= analysis_window.t_max) {
            ++n_aw_bin[ib];
        } else {
            ++n_tsb_bin[ib];
        }

        if (!IsInsideFullAnalysisWindow(Ee, Eg, theta_eg, t)) {
            hT[ib]->Fill(t);
            ++n_fit_eg_bin[ib];
            ++n_fit_used;
        }
        ++n_selected;
    }

    for (int i = 0; i < kNEgBins; ++i) {
        ConvertCountsToDensity(*hT[i]);
        StyleHist(*hT[i]);
    }

    const TString outpdf = MakeOutputPdfPath(infile);
    const TString outtxt = MakeSummaryTxtPath(infile);

    double r_fit[kNEgBins];
    double r_data[kNEgBins];
    double n_aw_fit[kNEgBins];
    double n_tsb_fit[kNEgBins];
    int fit_status[kNEgBins];
    double chi2_over_ndf[kNEgBins];
    for (int i = 0; i < kNEgBins; ++i) {
        r_fit[i] = 0.0;
        r_data[i] = 0.0;
        n_aw_fit[i] = 0.0;
        n_tsb_fit[i] = 0.0;
        fit_status[i] = -999;
        chi2_over_ndf[i] = 0.0;
    }

    TCanvas c0("c0_fit_acc_time_meta", "meta", 1200, 800);
    c0.cd();
    DrawMetaPage(infile, outpdf.Data(), n_lines, n_parsed, n_selected, n_fit_used, n_fit_eg_bin);

    c0.Print(Form("%s[", outpdf.Data()));
    c0.Print(outpdf.Data());

    for (int i = 0; i < kNEgBins; ++i) {
        TCanvas c(Form("c_fit_acc_time_%d", i), "fit", 1200, 900);

        TH1D* h = hT[i];
        const double hmax = std::max(1.0, 1.25 * h->GetMaximum());

        AccTimeFitResult fitres{};
        const int fit_ret =
            AccTimeFit_Run(*h, -gBlindCoreNs, gBlindCoreNs,
                           sigma_min_ns, 400.0, fitres);

        TF1 f_draw(Form("f_acc_time_draw_egbin_%d", i), AccTimeFit_DensityCore,
                   kTAllMin, kTAllMax, 3);
        f_draw.SetParameters(fitres.A, fitres.sigma, fitres.C);
        f_draw.SetLineColor(kRed + 1);
        f_draw.SetLineWidth(2);

        gPad->SetGrid();
        h->SetMaximum(hmax);
        h->Draw("E");

        TBox blind_box(-gBlindCoreNs, 0.0, gBlindCoreNs, hmax);
        blind_box.SetFillColorAlpha(kGray + 1, 0.25);
        blind_box.SetLineColor(kGray + 2);
        blind_box.Draw("same");

        f_draw.Draw("same");

        TLegend lg(0.58, 0.68, 0.88, 0.88);
        lg.SetBorderSize(0);
        lg.SetFillStyle(0);
        lg.AddEntry(h, "data", "lep");
        lg.AddEntry(&f_draw, "fit", "l");
        lg.AddEntry(&blind_box, "blind core", "f");
        lg.Draw();

        TLatex lat;
        lat.SetNDC(true);
        lat.SetTextAlign(13);
        lat.SetTextSize(0.030);
        lat.DrawLatex(0.12, 0.92,
                      Form("Eg in [%.0f, %.0f] MeV, entries=%lld",
                           kEgEdges[i], kEgEdges[i + 1], n_fit_eg_bin[i]));
        lat.DrawLatex(0.12, 0.86,
                      Form("fit status=%d, chi2/ndf=%.2f / %d",
                           fitres.fit_status,
                           fitres.chi2, fitres.ndf));
        lat.DrawLatex(0.12, 0.80,
                      Form("A=%.4g, sigma=%.2f ns, c0=%.4g",
                           fitres.A, fitres.sigma, fitres.C));

        c.Print(outpdf.Data());

        const double fit_aw =
            f_draw.Integral(analysis_window.t_min, analysis_window.t_max);
        const double fit_tsb =
            f_draw.Integral(kTAllMin, analysis_window.t_min) +
            f_draw.Integral(analysis_window.t_max, kTAllMax);
        const double data_side =
            HistIntegralRange(*h, -kTAllMax, -gBlindCoreNs) +
            HistIntegralRange(*h, gBlindCoreNs, kTAllMax);

        Info("fit_acc_time_by_eg",
             "Eg=[%.0f,%.0f] entries=%lld fit_status=%d chi2/ndf=%.3f/%d fit_AW=%.6g fit_TSB=%.6g side_mass(>|blind|)=%.6g",
             kEgEdges[i], kEgEdges[i + 1], n_fit_eg_bin[i], fitres.fit_status,
             fitres.chi2, fitres.ndf, fit_aw, fit_tsb, data_side);

        n_aw_fit[i] = fit_aw;
        n_tsb_fit[i] = fit_tsb;
        r_fit[i] = (fit_tsb > 0.0) ? (fit_aw / fit_tsb) : 0.0;
        r_data[i] = (n_tsb_bin[i] > 0) ? (static_cast<double>(n_aw_bin[i]) / static_cast<double>(n_tsb_bin[i])) : 0.0;
        fit_status[i] = (fit_ret == 0) ? fitres.fit_status : fit_ret;
        chi2_over_ndf[i] = fitres.chi2_ndf;
    }

    {
        std::ofstream fout(outtxt.Data());
        if (fout) {
            fout << "# Eg_low Eg_high entries N_AW_data N_TSB_data R_data fit_status chi2_ndf N_AW_fit N_TSB_fit R_fit\n";
            for (int i = 0; i < kNEgBins; ++i) {
                fout << kEgEdges[i] << " "
                     << kEgEdges[i + 1] << " "
                     << n_fit_eg_bin[i] << " "
                     << n_aw_bin[i] << " "
                     << n_tsb_bin[i] << " "
                     << r_data[i] << " "
                     << fit_status[i] << " "
                     << chi2_over_ndf[i] << " "
                     << n_aw_fit[i] << " "
                     << n_tsb_fit[i] << " "
                     << r_fit[i] << "\n";
            }
        }
    }

    TCanvas csum("c_fit_acc_time_summary", "summary", 1200, 800);
    csum.Divide(1, 2);

    TH1D hRFit("hRFit", "R = AW / TSB by Eg bin;Eg bin;R", kNEgBins, 0.0, static_cast<double>(kNEgBins));
    TH1D hRData("hRData", "R = AW / TSB by Eg bin;Eg bin;R", kNEgBins, 0.0, static_cast<double>(kNEgBins));
    TH1D hCounts("hCounts", "Counts by Eg bin;Eg bin;Entries", kNEgBins, 0.0, static_cast<double>(kNEgBins));

    hRFit.SetLineColor(kRed + 1);
    hRFit.SetMarkerColor(kRed + 1);
    hRFit.SetMarkerStyle(20);
    hRFit.SetLineWidth(2);

    hRData.SetLineColor(kBlue + 1);
    hRData.SetMarkerColor(kBlue + 1);
    hRData.SetMarkerStyle(21);
    hRData.SetLineWidth(2);

    hCounts.SetFillColorAlpha(kGray + 1, 0.35);
    hCounts.SetLineColor(kBlack);

    for (int i = 0; i < kNEgBins; ++i) {
        const int ib = i + 1;
        hRFit.GetXaxis()->SetBinLabel(ib, Form("%.0f-%.0f", kEgEdges[i], kEgEdges[i + 1]));
        hRData.GetXaxis()->SetBinLabel(ib, Form("%.0f-%.0f", kEgEdges[i], kEgEdges[i + 1]));
        hCounts.GetXaxis()->SetBinLabel(ib, Form("%.0f-%.0f", kEgEdges[i], kEgEdges[i + 1]));
        hRFit.SetBinContent(ib, r_fit[i]);
        hRData.SetBinContent(ib, r_data[i]);
        hCounts.SetBinContent(ib, n_fit_eg_bin[i]);
    }

    csum.cd(1);
    gPad->SetGrid();
    hRFit.SetMinimum(0.0);
    hRFit.SetMaximum(std::max(0.7, 1.15 * std::max(hRFit.GetMaximum(), hRData.GetMaximum())));
    hRFit.Draw("hist p");
    hRData.Draw("hist p same");
    TLegend lgR(0.62, 0.74, 0.88, 0.88);
    lgR.SetBorderSize(0);
    lgR.SetFillStyle(0);
    lgR.AddEntry(&hRFit, "R_fit from smooth model", "lp");
    lgR.AddEntry(&hRData, "R_data = N_AW/N_TSB", "lp");
    lgR.Draw();

    csum.cd(2);
    gPad->SetGrid();
    hCounts.Draw("hist");

    csum.Print(outpdf.Data());

    c0.Print(Form("%s]", outpdf.Data()));
    Info("fit_acc_time_by_eg", "wrote %s", outpdf.Data());
    Info("fit_acc_time_by_eg", "wrote %s", outtxt.Data());
}
