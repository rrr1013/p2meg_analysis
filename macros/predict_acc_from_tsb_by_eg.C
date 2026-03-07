// macros/predict_acc_from_tsb_by_eg.C
//
// 目的:
//   - fit_acc_time_by_eg.C で求めた Eg 依存の R_i = AW/TSB を使い、
//     TSB から ACC の AW 事象数を予測する。
//   - ここでは新しい fit はせず、既に目視確認した fit 結果テキストを読む。
//
// 入力:
//   1) 最終5Dデータ
//   2) fit_acc_time_by_eg.C の summary txt
//
// 出力:
//   - doc/finalanalysis/predict_acc_from_tsb_by_eg_<input>.pdf
//
// 実行例:
//   root -l -q 'macros/predict_acc_from_tsb_by_eg.C("data/finaldata/step4f_run1to20_eg_doubleonly_allpatterns_500ns.txt")'

R__ADD_INCLUDE_PATH(./include)

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

#include "TCanvas.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TError.h"

#include "p2meg/AnalysisWindow.h"
#include "p2meg/DetectorResolution.h"

static constexpr double kTAllMin = -500.0;   // [ns]
static constexpr double kTAllMax =  500.0;   // [ns]

struct FitSummaryRow {
    double eg_low = 0.0;
    double eg_high = 0.0;
    long long entries = 0;
    long long n_aw_data = 0;
    long long n_tsb_data = 0;
    double r_data = 0.0;
    int fit_status = -999;
    double chi2_ndf = 0.0;
    double n_aw_fit = 0.0;
    double n_tsb_fit = 0.0;
    double r_fit = 0.0;
};

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

static TString MakeDefaultFitTxtPath(const char* infile)
{
    TString base = gSystem->BaseName(infile);
    const Ssiz_t dot = base.Last('.');
    if (dot != kNPOS) base.Remove(dot);
    return Form("doc/finalanalysis/fit_acc_time_by_eg_%s.txt", base.Data());
}

static TString MakeOutputPdfPath(const char* infile)
{
    TString base = gSystem->BaseName(infile);
    const Ssiz_t dot = base.Last('.');
    if (dot != kNPOS) base.Remove(dot);
    gSystem->mkdir("doc/finalanalysis", /*recursive=*/true);
    return Form("doc/finalanalysis/predict_acc_from_tsb_by_eg_%s.pdf", base.Data());
}

static bool LoadFitSummary(const char* path, std::vector<FitSummaryRow>& rows)
{
    rows.clear();

    std::ifstream fin(path);
    if (!fin) return false;

    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') continue;

        std::istringstream iss(line);
        FitSummaryRow row;
        if (!(iss >> row.eg_low
                  >> row.eg_high
                  >> row.entries
                  >> row.n_aw_data
                  >> row.n_tsb_data
                  >> row.r_data
                  >> row.fit_status
                  >> row.chi2_ndf
                  >> row.n_aw_fit
                  >> row.n_tsb_fit
                  >> row.r_fit)) {
            continue;
        }
        rows.push_back(row);
    }
    return !rows.empty();
}

static int FindEgBin(const std::vector<FitSummaryRow>& rows, double Eg)
{
    if (!std::isfinite(Eg)) return -1;
    for (int i = 0; i < static_cast<int>(rows.size()); ++i) {
        const double lo = rows[i].eg_low;
        const double hi = rows[i].eg_high;
        if ((Eg >= lo && Eg < hi) ||
            (i == static_cast<int>(rows.size()) - 1 && Eg >= lo && Eg <= hi)) {
            return i;
        }
    }
    return -1;
}

static void DrawMetaPage(const char* infile,
                         const char* fit_txt,
                         const char* outpdf,
                         long long n_lines,
                         long long n_parsed,
                         long long n_selected,
                         double n_acc_pred_total,
                         double n_acc_pred_err_total,
                         double n_acc_excess_total)
{
    TLatex lat;
    lat.SetNDC(true);
    lat.SetTextAlign(13);

    lat.SetTextSize(0.040);
    lat.DrawLatex(0.05, 0.92, "ACC prediction from TSB using Eg-dependent R_i");

    lat.SetTextSize(0.030);
    lat.DrawLatex(0.05, 0.84, Form("input data : %s", infile));
    lat.DrawLatex(0.05, 0.79, Form("fit summary: %s", fit_txt));
    lat.DrawLatex(0.05, 0.74, Form("output     : %s", outpdf));

    lat.DrawLatex(0.05, 0.66, Form("lines read         : %lld", n_lines));
    lat.DrawLatex(0.05, 0.61, Form("parsed (5 doubles) : %lld", n_parsed));
    lat.DrawLatex(0.05, 0.56, Form("selected events    : %lld", n_selected));

    lat.DrawLatex(0.05, 0.46, Form("predicted ACC in AW (sum over Eg bins) : %.6g #pm %.6g", n_acc_pred_total, n_acc_pred_err_total));
    lat.DrawLatex(0.05, 0.41, Form("observed AW - predicted ACC            : %.6g", n_acc_excess_total));

    lat.DrawLatex(0.05, 0.33, "Error model:");
    lat.DrawLatex(0.08, 0.28, "#sigma^{2}(N_{acc,pred,i}) = (R_{fit,i}#sqrt{N_{TSB,i}})^{2} + (N_{TSB,i}|R_{fit,i}-R_{data,i}|)^{2}");

    lat.DrawLatex(0.05, 0.20, "Selection:");
    lat.DrawLatex(0.08, 0.15, Form("Ee in [%.1f, %.1f] MeV", analysis_window.Ee_min, analysis_window.Ee_max));
    lat.DrawLatex(0.08, 0.10, Form("Eg in [%.1f, %.1f] MeV", analysis_window.Eg_min, analysis_window.Eg_max));
    lat.DrawLatex(0.50, 0.15, Form("theta_eg in [%.3f, %.3f] rad", analysis_window.theta_min, analysis_window.theta_max));
    lat.DrawLatex(0.50, 0.10, Form("TSB: t in [%.0f, %.0f] ns and outside [%.0f, %.0f] ns",
                                   kTAllMin, kTAllMax, analysis_window.t_min, analysis_window.t_max));
}

void predict_acc_from_tsb_by_eg(
    const char* infile = "data/finaldata/step4f_run1to20_eg_doubleonly_allpatterns_500ns.txt",
    const char* fit_summary_txt = "")
{
    gStyle->SetOptStat(0);

    TString fit_txt = fit_summary_txt;
    if (fit_txt.IsNull() || fit_txt.Length() == 0) {
        fit_txt = MakeDefaultFitTxtPath(infile);
    }

    std::vector<FitSummaryRow> rows;
    if (!LoadFitSummary(fit_txt.Data(), rows)) {
        Error("predict_acc_from_tsb_by_eg", "failed to load fit summary: %s", fit_txt.Data());
        return;
    }

    const int nb = static_cast<int>(rows.size());
    std::vector<long long> n_aw(nb, 0LL);
    std::vector<long long> n_tsb(nb, 0LL);
    std::vector<double> n_acc_pred(nb, 0.0);
    std::vector<double> n_acc_pred_err(nb, 0.0);
    std::vector<double> n_excess(nb, 0.0);

    std::ifstream fin(infile);
    if (!fin) {
        Error("predict_acc_from_tsb_by_eg", "failed to open input file: %s", infile);
        return;
    }

    long long n_lines = 0;
    long long n_parsed = 0;
    long long n_selected = 0;

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
        if (Eg < analysis_window.Eg_min || Eg > analysis_window.Eg_max) continue;
        if (theta_eg < analysis_window.theta_min || theta_eg > analysis_window.theta_max) continue;

        const int ib = FindEgBin(rows, Eg);
        if (ib < 0) continue;

        ++n_selected;
        if (t >= analysis_window.t_min && t <= analysis_window.t_max) {
            ++n_aw[ib];
        } else {
            ++n_tsb[ib];
        }
    }

    double n_acc_pred_total = 0.0;
    double n_acc_pred_err2_total = 0.0;
    double n_aw_total = 0.0;
    for (int i = 0; i < nb; ++i) {
        n_acc_pred[i] = rows[i].r_fit * static_cast<double>(n_tsb[i]);
        const double err_stat = rows[i].r_fit * std::sqrt(static_cast<double>(n_tsb[i]));
        const double err_model = static_cast<double>(n_tsb[i]) * std::fabs(rows[i].r_fit - rows[i].r_data);
        n_acc_pred_err[i] = std::sqrt(err_stat * err_stat + err_model * err_model);
        n_excess[i] = static_cast<double>(n_aw[i]) - n_acc_pred[i];
        n_acc_pred_total += n_acc_pred[i];
        n_acc_pred_err2_total += n_acc_pred_err[i] * n_acc_pred_err[i];
        n_aw_total += static_cast<double>(n_aw[i]);

        Info("predict_acc_from_tsb_by_eg",
             "Eg=[%.0f,%.0f] N_AW=%lld N_TSB=%lld R_fit=%.6g R_data=%.6g => N_ACC_pred=%.6g +/- %.6g excess=%.6g",
             rows[i].eg_low, rows[i].eg_high, n_aw[i], n_tsb[i], rows[i].r_fit,
             rows[i].r_data, n_acc_pred[i], n_acc_pred_err[i], n_excess[i]);
    }
    const double n_acc_pred_err_total = std::sqrt(n_acc_pred_err2_total);

    const TString outpdf = MakeOutputPdfPath(infile);

    TCanvas c0("c0_predict_acc", "meta", 1200, 800);
    c0.cd();
    DrawMetaPage(infile, fit_txt.Data(), outpdf.Data(),
                 n_lines, n_parsed, n_selected,
                 n_acc_pred_total, n_acc_pred_err_total, n_aw_total - n_acc_pred_total);

    TH1D hObs("hObs", "Observed vs predicted counts in AW;Eg bin;Counts", nb, 0.0, static_cast<double>(nb));
    TH1D hPred("hPred", "Observed vs predicted counts in AW;Eg bin;Counts", nb, 0.0, static_cast<double>(nb));
    TH1D hExc("hExc", "AW excess after ACC prediction;Eg bin;Counts", nb, 0.0, static_cast<double>(nb));

    hObs.SetLineColor(kBlack);
    hObs.SetMarkerColor(kBlack);
    hObs.SetMarkerStyle(20);
    hObs.SetLineWidth(2);

    hPred.SetLineColor(kRed + 1);
    hPred.SetMarkerColor(kRed + 1);
    hPred.SetMarkerStyle(21);
    hPred.SetLineWidth(2);

    hExc.SetFillColorAlpha(kBlue + 1, 0.35);
    hExc.SetLineColor(kBlue + 2);
    hExc.SetLineWidth(2);

    for (int i = 0; i < nb; ++i) {
        const int ib = i + 1;
        const TString label = Form("%.0f-%.0f", rows[i].eg_low, rows[i].eg_high);
        hObs.GetXaxis()->SetBinLabel(ib, label.Data());
        hPred.GetXaxis()->SetBinLabel(ib, label.Data());
        hExc.GetXaxis()->SetBinLabel(ib, label.Data());

        hObs.SetBinContent(ib, static_cast<double>(n_aw[i]));
        hPred.SetBinContent(ib, n_acc_pred[i]);
        hPred.SetBinError(ib, n_acc_pred_err[i]);
        hExc.SetBinContent(ib, n_excess[i]);
    }

    TCanvas c1("c1_predict_acc", "counts", 1200, 800);
    c1.Divide(1, 2);

    c1.cd(1);
    gPad->SetGrid();
    hObs.SetMinimum(0.0);
    hObs.SetMaximum(1.2 * std::max(hObs.GetMaximum(), hPred.GetMaximum() + 1.0));
    hObs.Draw("hist p");
    hPred.Draw("E1 same");
    TLegend lg(0.62, 0.72, 0.88, 0.88);
    lg.SetBorderSize(0);
    lg.SetFillStyle(0);
    lg.AddEntry(&hObs, "observed N_AW", "lp");
    lg.AddEntry(&hPred, "predicted ACC from TSB", "lp");
    lg.Draw();

    c1.cd(2);
    gPad->SetGrid();
    hExc.Draw("hist");

    c0.Print(Form("%s[", outpdf.Data()));
    c0.Print(outpdf.Data());
    c1.Print(outpdf.Data());
    c1.Print(Form("%s]", outpdf.Data()));

    Info("predict_acc_from_tsb_by_eg", "total observed N_AW=%.6g", n_aw_total);
    Info("predict_acc_from_tsb_by_eg", "total predicted ACC N_AW=%.6g +/- %.6g", n_acc_pred_total, n_acc_pred_err_total);
    Info("predict_acc_from_tsb_by_eg", "total excess=%.6g", n_aw_total - n_acc_pred_total);
    Info("predict_acc_from_tsb_by_eg", "wrote %s", outpdf.Data());
}
