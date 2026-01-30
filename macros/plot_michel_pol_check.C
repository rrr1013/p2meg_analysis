// macros/plot_michel_pol_check.C
//
// Michel 偏極測定（別実験）チェック用プロット（PDF出力版 / 一発実行対応）
//
// 重要:
//  - root -l -q 'macros/plot_michel_pol_check.C("...")' で未解決シンボルを出さないため、
//    ReadMichelEData / BuildMichelPolKi の実体を「マクロと同一モジュール」に入れます。
//  - 具体的には src の .cc をこのマクロ内で #include します（ROOTマクロとして一般的な手）。
//
// 入力: Ee-only .dat（1行1事象、1列Ee[MeV]。空行/#コメント行は無視）
//   - A_plus  : Run#1 で +側（cosθ>0）
//   - A_minus : Run#2 で -側（cosθ<0）
//   - B_plus  : Run#2 で +側（cosθ>0）
//   - B_minus : Run#1 で -側（cosθ<0）
//
// 出力: doc/michel_pol_check_<basename(A_plus)>.pdf（複数ページ）
//   1) meta（入力/設定/統計/フィット結果）
//   2) A_data(E) と A_th(E)=P_hat*K(E)
//   3) K(E)
//   4) pull(E)（フィット範囲内）
//
// 実行例（リポジトリ直下から）:
//   root -l -q 'macros/plot_michel_pol_check.C("data/mockdata_michel/run1_A_plus.dat","data/mockdata_michel/run2_A_minus.dat","data/mockdata_michel/run2_B_plus.dat","data/mockdata_michel/run1_B_minus.dat")'

R__ADD_INCLUDE_PATH(./include)

#include <vector>
#include <string>
#include <cmath>

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TError.h"
#include "TString.h"
#include "TLatex.h"
#include "TPad.h"

// 設定・型の宣言（ヘッダ）
#include "p2meg/MichelEData.h"
#include "p2meg/MichelPolConfig.h"

// ============================================================
// ここがポイント：未解決を避けるため実装を同一モジュールに取り込む
// （相対パスは macros/ から見た ../src）
// ============================================================
#include "../src/MichelEData.cc"
#include "../src/MichelSpectrum.cc"
#include "../src/MichelPolTemplate.cc"


namespace {

// ---- utility（src 側の static 関数名と被らないように prefix を付ける）----
static inline int PlotMichelPol_BinIndex(double Ee, double Emin, double Emax, int nbins)
{
    if (!(Emax > Emin) || !(nbins > 0)) return -1;
    if (!(Ee >= Emin) || !(Ee < Emax)) return -1;
    const double binw = (Emax - Emin) / (double)nbins;
    const int idx = (int)((Ee - Emin) / binw);
    if (idx < 0 || idx >= nbins) return -1;
    return idx;
}

static inline double PlotMichelPol_BinCenter(int i, double Emin, double Emax, int nbins)
{
    const double binw = (Emax - Emin) / (double)nbins;
    return Emin + (i + 0.5) * binw;
}

static void PlotMichelPol_FillCounts(const std::vector<MichelEEvent>& evs,
                                     const MichelPolConfig& cfg,
                                     std::vector<double>& counts)
{
    for (const auto& ev : evs) {
        const int ib = PlotMichelPol_BinIndex(ev.Ee, cfg.Ee_min, cfg.Ee_max, cfg.nbins_Ee);
        if (ib >= 0) counts[ib] += 1.0;
    }
}

// cross-ratio 非対称と誤差（近似）
// a=nA+(E), b=nB+(E), c=nA-(E), d=nB-(E)
static bool PlotMichelPol_CrossRatioAsymAndErr(double a, double b, double c, double d,
                                              int min_counts_each,
                                              double& A, double& sigmaA)
{
    if (!(a >= 0.0 && b >= 0.0 && c >= 0.0 && d >= 0.0)) return false;
    if (a < min_counts_each || b < min_counts_each || c < min_counts_each || d < min_counts_each) return false;

    const double ab = a * b;
    const double cd = c * d;
    if (!(ab > 0.0 && cd > 0.0) || !std::isfinite(ab) || !std::isfinite(cd)) return false;

    const double r = std::sqrt(ab / cd);
    if (!(r > 0.0) || !std::isfinite(r)) return false;

    A = (r - 1.0) / (r + 1.0);
    if (!std::isfinite(A)) return false;

    // Var(ln r) ≈ (1/4)(1/a+1/b+1/c+1/d)
    // u=(1/2)ln r, A=tanh(u) => dA/du = 1-A^2
    const double var_ln_r = 0.25 * (1.0/a + 1.0/b + 1.0/c + 1.0/d);
    const double var_u    = 0.25 * var_ln_r;
    const double one_minus_A2 = 1.0 - A*A;
    const double var_A = (one_minus_A2 * one_minus_A2) * var_u;

    if (!(var_A > 0.0) || !std::isfinite(var_A)) return false;
    sigmaA = std::sqrt(var_A);
    return (sigmaA > 0.0 && std::isfinite(sigmaA));
}

static inline bool PlotMichelPol_InFitRange(double Ecen, const MichelPolConfig& cfg)
{
    return (Ecen >= cfg.fit_Ee_min && Ecen <= cfg.fit_Ee_max);
}

static TString PlotMichelPol_MakeOutputPdfPath(const char* path_A_plus)
{
    TString base = gSystem->BaseName(path_A_plus); // e.g. xxx.dat
    Ssiz_t dot = base.Last('.');
    if (dot != kNPOS) base.Remove(dot);

    gSystem->mkdir("doc", /*recursive=*/true);
    return Form("doc/michel_pol_check_%s.pdf", base.Data());
}

static void PlotMichelPol_DrawMetaPage(const char* outpdf,
                                       const char* A_plus, const char* A_minus,
                                       const char* B_plus, const char* B_minus,
                                       long long r1,long long s1,
                                       long long r2,long long s2,
                                       long long r3,long long s3,
                                       long long r4,long long s4,
                                       int n_bins_all, int n_bins_fit,
                                       double P_hat, double errP, double chi2, int ndf)
{
    const auto& cfg = michel_pol_config;

    TLatex lat;
    lat.SetNDC(true);
    lat.SetTextAlign(13);

    lat.SetTextSize(0.040);
    lat.DrawLatex(0.05, 0.92, "Michel polarization check (p2MEG: separate experiment)");

    lat.SetTextSize(0.030);
    lat.DrawLatex(0.05, 0.85, Form("output : %s", outpdf));

    lat.SetTextSize(0.028);
    lat.DrawLatex(0.05, 0.78, "inputs:");
    lat.DrawLatex(0.08, 0.74, Form("A_plus  : %s", A_plus));
    lat.DrawLatex(0.08, 0.70, Form("A_minus : %s", A_minus));
    lat.DrawLatex(0.08, 0.66, Form("B_plus  : %s", B_plus));
    lat.DrawLatex(0.08, 0.62, Form("B_minus : %s", B_minus));

    lat.SetTextSize(0.028);
    lat.DrawLatex(0.05, 0.54, "config:");
    lat.DrawLatex(0.08, 0.50, Form("|cos#theta| = %.6g", cfg.cos_theta_abs));
    lat.DrawLatex(0.08, 0.46, Form("Ee range   = [%.6g, %.6g] MeV  nbins=%d",
                                   cfg.Ee_min, cfg.Ee_max, cfg.nbins_Ee));
    lat.DrawLatex(0.08, 0.42, Form("fit range  = [%.6g, %.6g] MeV",
                                   cfg.fit_Ee_min, cfg.fit_Ee_max));
    lat.DrawLatex(0.08, 0.38, Form("min_counts_each (a,b,c,d) = %d", cfg.min_counts_each));

    lat.SetTextSize(0.028);
    lat.DrawLatex(0.05, 0.30, "I/O summary (Ee-only parse):");
    lat.DrawLatex(0.08, 0.26, Form("A_plus  : read=%lld  skipped=%lld",  r1, s1));
    lat.DrawLatex(0.08, 0.22, Form("A_minus : read=%lld  skipped=%lld",  r2, s2));
    lat.DrawLatex(0.08, 0.18, Form("B_plus  : read=%lld  skipped=%lld",  r3, s3));
    lat.DrawLatex(0.08, 0.14, Form("B_minus : read=%lld  skipped=%lld",  r4, s4));

    lat.SetTextSize(0.030);
    lat.DrawLatex(0.55, 0.54, "fit (A_data = P_mu K):");
    lat.SetTextSize(0.028);
    lat.DrawLatex(0.55, 0.50, Form("P_mu_hat = %.6f  +/- %.6f", P_hat, errP));
    if (ndf > 0) {
        lat.DrawLatex(0.55, 0.46, Form("chi2/ndf = %.2f / %d = %.3f", chi2, ndf, chi2 / (double)ndf));
    } else {
        lat.DrawLatex(0.55, 0.46, Form("chi2    = %.2f   ndf=%d", chi2, ndf));
    }
    lat.DrawLatex(0.55, 0.42, Form("bins with A_data (all) = %d", n_bins_all));
    lat.DrawLatex(0.55, 0.38, Form("bins used in fit       = %d", n_bins_fit));

    lat.SetTextSize(0.024);
    lat.DrawLatex(0.05, 0.04, "pages: (1) meta  (2) A_data vs A_th  (3) K(E)  (4) pulls");
}

} // namespace

void plot_michel_pol_check(const char* path_A_plus,
                           const char* path_A_minus,
                           const char* path_B_plus,
                           const char* path_B_minus)
{
    gStyle->SetOptStat(0);

    const TString outpdf = PlotMichelPol_MakeOutputPdfPath(path_A_plus);

    const auto& cfg = michel_pol_config;
    const int nb = cfg.nbins_Ee;
    if (!(nb > 0) || !(cfg.Ee_max > cfg.Ee_min)) {
        Error("plot_michel_pol_check", "invalid michel_pol_config.");
        return;
    }

    // ---- 読み込み ----
    long long r1=0,s1=0, r2=0,s2=0, r3=0,s3=0, r4=0,s4=0;
    const auto evAp = ReadMichelEData(path_A_plus,  &r1, &s1);
    const auto evAm = ReadMichelEData(path_A_minus, &r2, &s2);
    const auto evBp = ReadMichelEData(path_B_plus,  &r3, &s3);
    const auto evBm = ReadMichelEData(path_B_minus, &r4, &s4);

    std::vector<double> nAp(nb,0.0), nAm(nb,0.0), nBp(nb,0.0), nBm(nb,0.0);
    PlotMichelPol_FillCounts(evAp, cfg, nAp);
    PlotMichelPol_FillCounts(evAm, cfg, nAm);
    PlotMichelPol_FillCounts(evBp, cfg, nBp);
    PlotMichelPol_FillCounts(evBm, cfg, nBm);

    // ---- 理論テンプレ K_i（src/MichelPolTemplate.cc を取り込んである）----
    const std::vector<double> K = BuildMichelPolKi(cfg);
    if ((int)K.size() != nb) {
        Error("plot_michel_pol_check", "BuildMichelPolKi returned wrong size.");
        return;
    }

    struct Pt { double E; double A; double s; double K; };
    std::vector<Pt> pts_all;
    std::vector<Pt> pts_fit;

    double S1 = 0.0;
    double S2 = 0.0;

    for (int i = 0; i < nb; ++i) {
        const double Ecen = PlotMichelPol_BinCenter(i, cfg.Ee_min, cfg.Ee_max, cfg.nbins_Ee);

        // cross-ratio: a=nA+(i), b=nB+(i), c=nA-(i), d=nB-(i)
        const double a = nAp[i];
        const double b = nBp[i];
        const double c = nAm[i];
        const double d = nBm[i];

        double A=0.0, sA=0.0;
        if (!PlotMichelPol_CrossRatioAsymAndErr(a,b,c,d,cfg.min_counts_each,A,sA)) continue;

        const double Ki = K[i];
        if (!std::isfinite(Ki)) continue;

        pts_all.push_back(Pt{Ecen, A, sA, Ki});

        if (PlotMichelPol_InFitRange(Ecen, cfg)) {
            const double w = 1.0 / (sA*sA);
            S1 += Ki * A * w;
            S2 += Ki * Ki * w;
            pts_fit.push_back(Pt{Ecen, A, sA, Ki});
        }
    }

    if (!(S2 > 0.0) || (int)pts_fit.size() < 2) {
        Error("plot_michel_pol_check", "not enough bins for fit.");
        return;
    }

    const double P_hat = S1 / S2;
    const double errP  = 1.0 / std::sqrt(S2);

    double chi2 = 0.0;
    for (const auto& p : pts_fit) {
        const double pull = (p.A - P_hat * p.K) / p.s;
        chi2 += pull * pull;
    }
    const int ndf = (int)pts_fit.size() - 1;

    // 1) meta
    TCanvas c0("c_michel_meta", "meta", 1200, 800);
    c0.cd();
    PlotMichelPol_DrawMetaPage(outpdf.Data(),
                               path_A_plus, path_A_minus, path_B_plus, path_B_minus,
                               r1,s1,r2,s2,r3,s3,r4,s4,
                               (int)pts_all.size(), (int)pts_fit.size(),
                               P_hat, errP, chi2, ndf);

    // 2) A_data vs A_th
    TCanvas c1("c_michel_asym", "asymmetry", 1200, 800);
    c1.cd();
    gPad->SetGrid();

    TGraphErrors gA;
    gA.SetTitle("Michel polarization: A_{data}(E) vs A_{th}(E);E_{e} [MeV];A(E)");
    gA.SetMarkerStyle(20);
    gA.SetMarkerSize(0.9);

    for (int ip = 0; ip < (int)pts_all.size(); ++ip) {
        gA.SetPoint(ip, pts_all[ip].E, pts_all[ip].A);
        gA.SetPointError(ip, 0.0, pts_all[ip].s);
    }

    TGraph gTh;
    gTh.SetLineWidth(2);
    for (int ip = 0; ip < (int)pts_all.size(); ++ip) {
        gTh.SetPoint(ip, pts_all[ip].E, P_hat * pts_all[ip].K);
    }

    gA.Draw("AP");
    gTh.Draw("L SAME");

    {
        const double ylo = gPad->GetUymin();
        const double yhi = gPad->GetUymax();
        TLine lfit1(cfg.fit_Ee_min, ylo, cfg.fit_Ee_min, yhi);
        TLine lfit2(cfg.fit_Ee_max, ylo, cfg.fit_Ee_max, yhi);
        lfit1.SetLineStyle(2);
        lfit2.SetLineStyle(2);
        lfit1.Draw("SAME");
        lfit2.Draw("SAME");
    }

    TLegend leg(0.58, 0.74, 0.88, 0.88);
    leg.AddEntry(&gA,  "A_{data} (cross-ratio)", "pe");
    leg.AddEntry(&gTh, "A_{th}=#hat{P}_{#mu} K(E)", "l");
    leg.AddEntry((TObject*)0, Form("#hat{P}_{#mu}=%.4f #pm %.4f", P_hat, errP), "");
    if (ndf > 0) leg.AddEntry((TObject*)0, Form("#chi^{2}/ndf=%.2f/%d", chi2, ndf), "");
    leg.Draw();

    // 3) K(E)
    TCanvas c2("c_michel_K", "K(E)", 1200, 800);
    c2.cd();
    gPad->SetGrid();

    TGraph gK;
    gK.SetTitle("Sensitivity curve K(E);E_{e} [MeV];K(E)");
    gK.SetLineWidth(2);

    int ipk = 0;
    for (int i = 0; i < nb; ++i) {
        const double Ecen = PlotMichelPol_BinCenter(i, cfg.Ee_min, cfg.Ee_max, cfg.nbins_Ee);
        const double Ki = K[i];
        if (!std::isfinite(Ki)) continue;
        gK.SetPoint(ipk++, Ecen, Ki);
    }
    gK.Draw("AL");

    // 4) pull(E) in fit range
    TCanvas c3("c_michel_pull", "pulls", 1200, 800);
    c3.cd();
    gPad->SetGrid();

    TGraph gP;
    gP.SetTitle("Pulls in fit range;E_{e} [MeV];(A_{data}-A_{th})/#sigma_{A}");
    gP.SetMarkerStyle(20);
    gP.SetMarkerSize(0.9);

    int ipp = 0;
    for (const auto& p : pts_fit) {
        const double pull = (p.A - P_hat * p.K) / p.s;
        gP.SetPoint(ipp++, p.E, pull);
    }
    gP.Draw("AP");

    {
        TLine l0(cfg.fit_Ee_min, 0.0, cfg.fit_Ee_max, 0.0);
        l0.SetLineStyle(2);
        l0.Draw("SAME");
    }

    // ---- PDF（複数ページ）----
    c0.Print(Form("%s[", outpdf.Data()));
    c0.Print(outpdf.Data());
    c1.Print(outpdf.Data());
    c2.Print(outpdf.Data());
    c3.Print(outpdf.Data());
    c3.Print(Form("%s]", outpdf.Data()));

    Info("plot_michel_pol_check", "wrote: %s", outpdf.Data());
}
