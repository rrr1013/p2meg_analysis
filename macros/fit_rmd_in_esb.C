// macros/fit_rmd_in_esb.C
//
// ESB(Eg) の t 分布をフィットして RMD と ACC を分離し、
// RMDExtrapolation (MC) で k を推定し、AW 内の RMD 数を予測する。
//
// 実行例:
//   root -l -q 'macros/fit_rmd_in_esb.C("data/mockdata/testdata1.dat")'
//

R__ADD_INCLUDE_PATH(./include)

#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <algorithm>

#include "TCanvas.h"
#include "TH1D.h"
#include "TF1.h"
#include "TFitResultPtr.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"
#include "TLatex.h"
#include "TMath.h"

#include "p2meg/AnalysisWindow.h"
#include "p2meg/AnalysisWindowUtils.h"
#include "p2meg/DetectorResolution.h"
#include "p2meg/MathUtils.h"
#include "p2meg/RMDExtrapolation.h"

// 単独実行時に未定義シンボルにならないよう、実装を取り込む
#include "../src/RMDSpectrum.cc"
#include "../src/RMDExtrapolation.cc"

// 5列 double のみ event として読む（それ以外はメタ行として捨てる）
static bool ParseEventLine5Doubles(const std::string& line,
                                  double& Ee, double& Eg, double& t,
                                  double& phi_e, double& phi_g)
{
  if (line.empty()) return false;
  std::istringstream iss(line);
  if (!(iss >> Ee >> Eg >> t >> phi_e >> phi_g)) return false;
  if (!std::isfinite(Ee) || !std::isfinite(Eg) || !std::isfinite(t) ||
      !std::isfinite(phi_e) || !std::isfinite(phi_g)) return false;
  return true;
}

static TString MakeOutputPdfPath(const char* infile)
{
  TString base = gSystem->BaseName(infile);
  Ssiz_t dot = base.Last('.');
  if (dot != kNPOS) base.Remove(dot);
  gSystem->mkdir("doc", /*recursive=*/true);
  return Form("doc/fit_rmd_in_esb_%s.pdf", base.Data());
}

// 正規化ガウシアンの区間積分
static double GaussFrac(double t1, double t2, double mu, double sigma)
{
  if (!(sigma > 0.0) || !std::isfinite(mu) || !std::isfinite(sigma)) return 0.0;
  if (t2 < t1) std::swap(t1, t2);
  const double a = (t1 - mu) / (std::sqrt(2.0) * sigma);
  const double b = (t2 - mu) / (std::sqrt(2.0) * sigma);
  return 0.5 * (TMath::Erf(b) - TMath::Erf(a));
}

void fit_rmd_in_esb(const char* infile = "data/mockdata/testdata1.dat",
                    double t_all_min = -10.0,
                    double t_all_max =  10.0,
                    int nbins_t = 240,
                    double Eg_all_min = 10.0,
                    double Eg_all_max = 70.0,
                    long long n_mc_k = 300000,
                    unsigned long long seed_k = 1)
{
  gStyle->SetOptStat(0);

  if (!(t_all_max > t_all_min)) {
    ::Error("fit_rmd_in_esb", "invalid t_all range");
    return;
  }
  if (nbins_t <= 10) nbins_t = 10;
  if (!std::isfinite(Eg_all_min) || !std::isfinite(Eg_all_max) || !(Eg_all_max > Eg_all_min)) {
    ::Error("fit_rmd_in_esb", "invalid Eg_all range");
    return;
  }

  // t window（AW）
  const double t_win_min = analysis_window.t_min;
  const double t_win_max = analysis_window.t_max;

  // 入力読み込み
  std::ifstream fin(infile);
  if (!fin) {
    ::Error("fit_rmd_in_esb", "failed to open input file: %s", infile);
    return;
  }

  TH1D* hT = new TH1D("hT_esb", "t in ESB;t [ns];Entries", nbins_t, t_all_min, t_all_max);
  hT->Sumw2();

  long long n_lines = 0;
  long long n_parsed = 0;
  long long n_esb = 0;

  const int nphi_e = Math_GetNPhiE(detres);
  const int nphi_g = Math_GetNPhiG(detres);

  std::string line;
  while (std::getline(fin, line)) {
    ++n_lines;

    double Ee = 0.0, Eg = 0.0, t = 0.0, phi_e = 0.0, phi_g = 0.0;
    if (!ParseEventLine5Doubles(line, Ee, Eg, t, phi_e, phi_g)) continue;
    ++n_parsed;

    // phi を検出器グリッドへ丸め
    const int idx_e = Detector_PhiIndexFromValue(phi_e, detres.phi_e_min, detres.phi_e_max, nphi_e);
    const int idx_g = Detector_PhiIndexFromValue(phi_g, detres.phi_g_min, detres.phi_g_max, nphi_g);
    if (idx_e < 0 || idx_g < 0) continue;
    if (!Detector_IsAllowedPhiPairIndex(idx_e, idx_g, detres)) continue;

    const double pe = Detector_PhiGridPoint(idx_e, detres.phi_e_min, detres.phi_e_max, nphi_e);
    const double pg = Detector_PhiGridPoint(idx_g, detres.phi_g_min, detres.phi_g_max, nphi_g);

    // MakeRMDGridPdf と同じ：theta_eg = |phi_e - phi_g|
    const double theta = std::fabs(pe - pg);

    // ESB（Eg sideband）判定は AnalysisWindowUtils の定義を使用
    if (!AnalysisWindow_InESB_Eg(analysis_window,
                                 Ee, Eg, theta, t,
                                 Eg_all_min, Eg_all_max,
                                 t_all_min, t_all_max)) {
      continue;
    }

    ++n_esb;
    hT->Fill(t);
  }

  if (hT->GetEntries() <= 0) {
    ::Warning("fit_rmd_in_esb", "no ESB entries in the specified t range. Nothing to fit.");
    return;
  }

  // ==========================================================
  // ヒスト bin 内容（カウント）に対して、
  // “bin あたり期待カウント” を直接フィットする。
  //   f_bin(t) = Δt * ( N_rmd * GausN + N_acc / T )
  // ==========================================================
  const double binw = hT->GetBinWidth(1);
  if (!(binw > 0.0) || !std::isfinite(binw)) {
    ::Error("fit_rmd_in_esb", "invalid bin width");
    return;
  }

  TF1* f = new TF1("f_rmdacc_bin",
                   "[6]*([0]*TMath::Gaus(x,[1],[2],true) + [3]/([5]-[4]))",
                   t_all_min, t_all_max);

  f->SetParameter(4, t_all_min);
  f->SetParameter(5, t_all_max);
  f->SetParameter(6, binw);
  f->FixParameter(4, t_all_min);
  f->FixParameter(5, t_all_max);
  f->FixParameter(6, binw);

  const double Ntot = hT->Integral(1, hT->GetNbinsX());
  const double mu0 = detres.t_mean;
  const double sig0 = (detres.sigma_t > 0.0 ? detres.sigma_t : 0.2);

  f->SetParName(0, "N_RMD");
  f->SetParName(1, "mu");
  f->SetParName(2, "sigma");
  f->SetParName(3, "N_ACC");

  f->SetParameter(0, 0.2 * Ntot);
  f->SetParameter(1, mu0);
  f->SetParameter(2, sig0);
  f->SetParameter(3, 0.8 * Ntot);

  f->SetParLimits(0, 0.0, 10.0 * Ntot);
  f->SetParLimits(2, 1e-6, 10.0);
  f->SetParLimits(3, 0.0, 10.0 * Ntot);

  (void)hT->Fit(f, "SRQ0");

  // ここでの N_rmd_esb, N_acc_esb は「ESB の事象数」
  const double N_rmd_esb = f->GetParameter(0);
  const double mu_fit    = f->GetParameter(1);
  const double sig_fit   = f->GetParameter(2);
  const double N_acc_esb = f->GetParameter(3);

  const double eN_rmd_esb = f->GetParError(0);
  const double eN_acc_esb = f->GetParError(3);

  const double N_sum_esb = N_rmd_esb + N_acc_esb;

  // time window 内に入る割合（fit の mu, sigma を使う）
  const double frac_time = GaussFrac(t_win_min, t_win_max, mu_fit, sig_fit);
  const double N_rmd_esb_in_tw  = N_rmd_esb * frac_time;
  const double eN_rmd_esb_in_tw = eN_rmd_esb * frac_time; // mu/sigma 誤差はここでは無視

  // k の推定（RMDExtrapolation.cc を使用）
  const KFactorResult kres = RMD_EstimateK_AW_over_ESB_MC(analysis_window,
                                                          detres,
                                                          Eg_all_min, Eg_all_max,
                                                          t_all_min, t_all_max,
                                                          n_mc_k,
                                                          (std::uint64_t)seed_k);

  const double k  = kres.k;
  const double ek = kres.sigma_k;

  // AW 内の RMD 予測数
  const double N_rmd_aw_pred  = k * N_rmd_esb_in_tw;
  const double eN_rmd_aw_pred = std::sqrt(std::pow(k * eN_rmd_esb_in_tw, 2) +
                                          std::pow(N_rmd_esb_in_tw * ek, 2));

  // 端末出力
  ::Info("fit_rmd_in_esb", "input: %s", infile);
  ::Info("fit_rmd_in_esb", "lines=%lld parsed=%lld ESB(entries)=%lld", n_lines, n_parsed, n_esb);
  ::Info("fit_rmd_in_esb", "t_all=[%.6g, %.6g], t_win=[%.6g, %.6g]", t_all_min, t_all_max, t_win_min, t_win_max);
  ::Info("fit_rmd_in_esb", "Eg_all=[%.6g, %.6g]", Eg_all_min, Eg_all_max);
  ::Info("fit_rmd_in_esb", "binw=%.6g", binw);

  ::Info("fit_rmd_in_esb", "Fit (t in ESB, count-per-bin model):");
  ::Info("fit_rmd_in_esb", "  N_RMD_ESB = %.6g +/- %.6g", N_rmd_esb, eN_rmd_esb);
  ::Info("fit_rmd_in_esb", "  N_ACC_ESB = %.6g +/- %.6g", N_acc_esb, eN_acc_esb);
  ::Info("fit_rmd_in_esb", "  N_SUM_ESB = %.6g (check vs entries %.6g)", N_sum_esb, (double)hT->GetEntries());
  ::Info("fit_rmd_in_esb", "  mu = %.6g  sigma = %.6g", mu_fit, sig_fit);
  ::Info("fit_rmd_in_esb", "  frac_time(t_win) = %.6g", frac_time);
  ::Info("fit_rmd_in_esb", "  N_RMD_ESB(in t_win) = %.6g +/- %.6g", N_rmd_esb_in_tw, eN_rmd_esb_in_tw);

  ::Info("fit_rmd_in_esb", "k estimation (MC, RMDExtrapolation):");
  ::Info("fit_rmd_in_esb", "  k = %.6g +/- %.6g", k, ek);
  ::Info("fit_rmd_in_esb", "  (n_gen=%lld, n_aw=%lld, n_esb=%lld)", kres.n_gen, kres.n_aw, kres.n_esb);

  ::Info("fit_rmd_in_esb", "RMD prediction in AW:");
  ::Info("fit_rmd_in_esb", "  N_RMD_AW_pred = %.6g +/- %.6g", N_rmd_aw_pred, eN_rmd_aw_pred);

  // PDF 出力
  const TString outpdf = MakeOutputPdfPath(infile);

  TCanvas c("c_fit", "fit_rmd_in_esb", 1100, 800);
  c.cd();
  hT->SetLineWidth(2);
  hT->Draw("E");

  f->SetLineWidth(2);
  f->Draw("same");

  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextAlign(13);
  lat.SetTextSize(0.030);
  lat.DrawLatex(0.12, 0.92, Form("input: %s", gSystem->BaseName(infile)));
  lat.DrawLatex(0.12, 0.88, Form("N_{RMD}^{ESB}=%.3g #pm %.3g", N_rmd_esb, eN_rmd_esb));
  lat.DrawLatex(0.12, 0.84, Form("N_{ACC}^{ESB}=%.3g #pm %.3g", N_acc_esb, eN_acc_esb));
  lat.DrawLatex(0.12, 0.80, Form("frac_{twin}=%.3g  ->  N_{RMD}^{ESB}(twin)=%.3g", frac_time, N_rmd_esb_in_tw));
  lat.DrawLatex(0.12, 0.76, Form("k=%.3g #pm %.3g", k, ek));
  lat.DrawLatex(0.12, 0.72, Form("N_{RMD}^{AW,pred}=%.3g #pm %.3g", N_rmd_aw_pred, eN_rmd_aw_pred));

  c.Print(outpdf);
  ::Info("fit_rmd_in_esb", "wrote: %s", outpdf.Data());
}
