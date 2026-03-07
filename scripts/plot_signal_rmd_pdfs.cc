// scripts/plot_signal_rmd_pdfs.cc
//
// 現在の解析窓・分解能設定で Signal PDF と RMD PDF のヒストを作る。
//
// 出力:
//   - doc/finalanalysis/signal_rmd_pdf_hist_current.pdf   : 4ページ
//   - doc/finalanalysis/signal_rmd_pdf_hist_current.root  : ヒスト保存
//   - doc/finalanalysis/signal_rmd_pdf_hist_current_*.png : 各ページPNG
//
// 注意:
//   - RMD は「現在の設定」で格子 PDF を /tmp に再生成してから読む。
//   - ヒストは PDF の数値積分から作る。1D/2D とも密度表示に直している。

#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"

#include "p2meg/AnalysisWindow.h"
#include "p2meg/Constants.h"
#include "p2meg/DetectorResolution.h"
#include "p2meg/MakeRMDGridPdf.h"
#include "p2meg/MathUtils.h"
#include "p2meg/RMDGridPdf.h"
#include "p2meg/SignalPdf.h"

static constexpr int kNBinsE = 60;
static constexpr int kNBinsT = 120;
static constexpr int kNBinsTheta = 18;

static constexpr int kColorSignal = kRed + 1;
static constexpr int kColorRmd = kBlue + 1;

static const char* kOutPdf = "doc/finalanalysis/signal_rmd_pdf_hist_current.pdf";
static const char* kOutRoot = "doc/finalanalysis/signal_rmd_pdf_hist_current.root";
static const char* kOutPngMeta = "doc/finalanalysis/signal_rmd_pdf_hist_current_meta.png";
static const char* kOutPng1D = "doc/finalanalysis/signal_rmd_pdf_hist_current_1d.png";
static const char* kOutPngSig2D = "doc/finalanalysis/signal_rmd_pdf_hist_current_signal2d.png";
static const char* kOutPngRmd2D = "doc/finalanalysis/signal_rmd_pdf_hist_current_rmd2d.png";
static const char* kTmpRmdGrid = "data/pdf_cache/rmd_grid.root";
static const char* kTmpRmdKey = "rmd_grid";

static bool IsFinite(double x)
{
  return std::isfinite(x);
}

static void ConvertHist1ToDensity(TH1D* h)
{
  if (!h) return;
  for (int ib = 1; ib <= h->GetNbinsX(); ++ib) {
    const double w = h->GetXaxis()->GetBinWidth(ib);
    if (!(w > 0.0)) continue;
    h->SetBinContent(ib, h->GetBinContent(ib) / w);
    h->SetBinError(ib, h->GetBinError(ib) / w);
  }
}

static void ConvertHist2ToDensity(TH2D* h)
{
  if (!h) return;
  for (int ix = 1; ix <= h->GetNbinsX(); ++ix) {
    const double wx = h->GetXaxis()->GetBinWidth(ix);
    if (!(wx > 0.0)) continue;
    for (int iy = 1; iy <= h->GetNbinsY(); ++iy) {
      const double wy = h->GetYaxis()->GetBinWidth(iy);
      if (!(wy > 0.0)) continue;
      const double area = wx * wy;
      h->SetBinContent(ix, iy, h->GetBinContent(ix, iy) / area);
      h->SetBinError(ix, iy, h->GetBinError(ix, iy) / area);
    }
  }
}

static double ThetaDisplayMin()
{
  double th_det_min = 0.0;
  double th_det_max = 0.0;
  if (Detector_ThetaRangeFromAllowedPhi(detres, th_det_min, th_det_max)) {
    const double th_min = (analysis_window.theta_min > th_det_min)
                            ? analysis_window.theta_min
                            : th_det_min;
    if (th_min < analysis_window.theta_max) return th_min;
  }
  return analysis_window.theta_min;
}

static double ThetaDisplayMax()
{
  double th_det_min = 0.0;
  double th_det_max = 0.0;
  if (Detector_ThetaRangeFromAllowedPhi(detres, th_det_min, th_det_max)) {
    const double th_max = (analysis_window.theta_max < th_det_max)
                            ? analysis_window.theta_max
                            : th_det_max;
    if (analysis_window.theta_min < th_max) return th_max;
  }
  return analysis_window.theta_max;
}

static void Style1D(TH1D* h, int color)
{
  if (!h) return;
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetLineWidth(2);
  h->SetStats(0);
}

static void WriteHistPair(TFile& fout,
                          TH1D* h_sig,
                          TH1D* h_rmd,
                          TH2D* h2_sig,
                          TH2D* h2_rmd)
{
  fout.cd();
  if (h_sig) h_sig->Write();
  if (h_rmd) h_rmd->Write();
  if (h2_sig) h2_sig->Write();
  if (h2_rmd) h2_rmd->Write();
}

static TH1D* MakePhiHist1D(const char* name, const char* title,
                           const std::vector<double>& edges)
{
  return new TH1D(name, title, static_cast<int>(edges.size()) - 1, edges.data());
}

static TH2D* MakePhiHist2D(const char* name, const char* title,
                           const std::vector<double>& edges_e,
                           const std::vector<double>& edges_g)
{
  return new TH2D(name, title,
                  static_cast<int>(edges_e.size()) - 1, edges_e.data(),
                  static_cast<int>(edges_g.size()) - 1, edges_g.data());
}

static void DrawMetaPage(double sig_norm,
                         double rmd_norm,
                         long long n_cells,
                         long long n_sig_pos,
                         long long n_rmd_pos)
{
  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextAlign(13);

  lat.SetTextSize(0.040);
  lat.DrawLatex(0.05, 0.93, "p2MEG Signal / RMD PDF histograms (current settings)");

  lat.SetTextSize(0.030);
  lat.DrawLatex(0.05, 0.85, Form("output pdf : %s", kOutPdf));
  lat.DrawLatex(0.05, 0.81, Form("output root: %s", kOutRoot));
  lat.DrawLatex(0.05, 0.77, Form("RMD grid    : %s", kTmpRmdGrid));

  lat.DrawLatex(0.05, 0.70, Form("analysis window: Ee=[%.6g, %.6g] MeV, Eg=[%.6g, %.6g] MeV",
                                 analysis_window.Ee_min, analysis_window.Ee_max,
                                 analysis_window.Eg_min, analysis_window.Eg_max));
  lat.DrawLatex(0.05, 0.66, Form("analysis window: t=[%.6g, %.6g] ns, theta=[%.6g, %.6g] rad",
                                 analysis_window.t_min, analysis_window.t_max,
                                 analysis_window.theta_min, analysis_window.theta_max));
  lat.DrawLatex(0.05, 0.62, Form("detres: sigma_t=%.6g ns, t_mean=%.6g ns, P_mu=%.6g, N_theta=%d",
                                 detres.sigma_t, detres.t_mean, detres.P_mu, detres.N_theta));
  lat.DrawLatex(0.05, 0.58, Form("phi_e=[%.6g, %.6g] rad (N=%d), phi_g=[%.6g, %.6g] rad (N=%d)",
                                 detres.phi_e_min, detres.phi_e_max, detres.N_phi_e,
                                 detres.phi_g_min, detres.phi_g_max, detres.N_phi_g));

  lat.DrawLatex(0.05, 0.50, Form("integration cells: %lld", n_cells));
  lat.DrawLatex(0.05, 0.46, Form("signal positive cells: %lld", n_sig_pos));
  lat.DrawLatex(0.05, 0.42, Form("RMD positive cells   : %lld", n_rmd_pos));

  lat.DrawLatex(0.05, 0.34, Form("signal norm check = %.10f", sig_norm));
  lat.DrawLatex(0.05, 0.30, Form("RMD norm check    = %.10f", rmd_norm));

  lat.DrawLatex(0.05, 0.20, "1D/2D histograms are converted to density.");
  lat.DrawLatex(0.05, 0.16, "Theta is derived from |phi_{e}-phi_{g}| after detector-grid snapping.");
  lat.DrawLatex(0.05, 0.12, "Current theta window keeps only 120 deg and 180 deg detector pairs.");
}

int main()
{
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat(".3g");

  gSystem->mkdir("doc/finalanalysis", true);

  const int N_phi_e = Math_GetNPhiE(detres);
  const int N_phi_g = Math_GetNPhiG(detres);
  const std::vector<double> phi_edges_e =
      Detector_PhiEdgesFromGrid(detres.phi_e_min, detres.phi_e_max, N_phi_e);
  const std::vector<double> phi_edges_g =
      Detector_PhiEdgesFromGrid(detres.phi_g_min, detres.phi_g_max, N_phi_g);

  if (MakeRMDGridPdfWithTruthWindow(kTmpRmdGrid, kTmpRmdKey, analysis_window) != 0) {
    std::cerr << "[plot_signal_rmd_pdfs] failed to generate current RMD grid\n";
    return 1;
  }
  if (!RMDGridPdf_Load(kTmpRmdGrid, kTmpRmdKey)) {
    std::cerr << "[plot_signal_rmd_pdfs] failed to load current RMD grid\n";
    return 1;
  }

  const double theta_min = ThetaDisplayMin();
  const double theta_max = ThetaDisplayMax();
  const double theta_axis_max = Math_AxisMaxInclusive(theta_max);

  TH1D* h_sig_Ee = new TH1D("h_sig_Ee", "Signal;Ee [MeV];Density [1/MeV]",
                            kNBinsE, analysis_window.Ee_min, analysis_window.Ee_max);
  TH1D* h_rmd_Ee = new TH1D("h_rmd_Ee", "RMD;Ee [MeV];Density [1/MeV]",
                            kNBinsE, analysis_window.Ee_min, analysis_window.Ee_max);

  TH1D* h_sig_Eg = new TH1D("h_sig_Eg", "Signal;Eg [MeV];Density [1/MeV]",
                            kNBinsE, analysis_window.Eg_min, analysis_window.Eg_max);
  TH1D* h_rmd_Eg = new TH1D("h_rmd_Eg", "RMD;Eg [MeV];Density [1/MeV]",
                            kNBinsE, analysis_window.Eg_min, analysis_window.Eg_max);

  TH1D* h_sig_t = new TH1D("h_sig_t", "Signal;#Deltat [ns];Density [1/ns]",
                           kNBinsT, analysis_window.t_min, Math_AxisMaxInclusive(analysis_window.t_max));
  TH1D* h_rmd_t = new TH1D("h_rmd_t", "RMD;#Deltat [ns];Density [1/ns]",
                           kNBinsT, analysis_window.t_min, Math_AxisMaxInclusive(analysis_window.t_max));

  TH1D* h_sig_phi_e = MakePhiHist1D("h_sig_phi_e",
                                    "Signal;#phi_{detector,e} [rad];Density [1/rad]",
                                    phi_edges_e);
  TH1D* h_rmd_phi_e = MakePhiHist1D("h_rmd_phi_e",
                                    "RMD;#phi_{detector,e} [rad];Density [1/rad]",
                                    phi_edges_e);

  TH1D* h_sig_phi_g = MakePhiHist1D("h_sig_phi_g",
                                    "Signal;#phi_{detector,#gamma} [rad];Density [1/rad]",
                                    phi_edges_g);
  TH1D* h_rmd_phi_g = MakePhiHist1D("h_rmd_phi_g",
                                    "RMD;#phi_{detector,#gamma} [rad];Density [1/rad]",
                                    phi_edges_g);

  TH1D* h_sig_theta = new TH1D("h_sig_theta",
                               "Signal;#theta_{eg} [rad];Density [1/rad]",
                               kNBinsTheta, theta_min, theta_axis_max);
  TH1D* h_rmd_theta = new TH1D("h_rmd_theta",
                               "RMD;#theta_{eg} [rad];Density [1/rad]",
                               kNBinsTheta, theta_min, theta_axis_max);

  TH2D* h2_sig_EeEg = new TH2D("h2_sig_EeEg",
                               "Signal;Ee [MeV];Eg [MeV]",
                               kNBinsE, analysis_window.Ee_min, analysis_window.Ee_max,
                               kNBinsE, analysis_window.Eg_min, analysis_window.Eg_max);
  TH2D* h2_rmd_EeEg = new TH2D("h2_rmd_EeEg",
                               "RMD;Ee [MeV];Eg [MeV]",
                               kNBinsE, analysis_window.Ee_min, analysis_window.Ee_max,
                               kNBinsE, analysis_window.Eg_min, analysis_window.Eg_max);

  TH2D* h2_sig_phi = MakePhiHist2D("h2_sig_phi",
                                   "Signal;#phi_{detector,e} [rad];#phi_{detector,#gamma} [rad]",
                                   phi_edges_e, phi_edges_g);
  TH2D* h2_rmd_phi = MakePhiHist2D("h2_rmd_phi",
                                   "RMD;#phi_{detector,e} [rad];#phi_{detector,#gamma} [rad]",
                                   phi_edges_e, phi_edges_g);

  Style1D(h_sig_Ee, kColorSignal);
  Style1D(h_rmd_Ee, kColorRmd);
  Style1D(h_sig_Eg, kColorSignal);
  Style1D(h_rmd_Eg, kColorRmd);
  Style1D(h_sig_t, kColorSignal);
  Style1D(h_rmd_t, kColorRmd);
  Style1D(h_sig_phi_e, kColorSignal);
  Style1D(h_rmd_phi_e, kColorRmd);
  Style1D(h_sig_phi_g, kColorSignal);
  Style1D(h_rmd_phi_g, kColorRmd);
  Style1D(h_sig_theta, kColorSignal);
  Style1D(h_rmd_theta, kColorRmd);

  double sig_norm = 0.0;
  double rmd_norm = 0.0;
  long long n_cells = 0;
  long long n_sig_pos = 0;
  long long n_rmd_pos = 0;

  for (int ib_e = 1; ib_e <= h_sig_Ee->GetNbinsX(); ++ib_e) {
    const double Ee = h_sig_Ee->GetXaxis()->GetBinCenter(ib_e);
    const double dEe = h_sig_Ee->GetXaxis()->GetBinWidth(ib_e);
    for (int ib_g = 1; ib_g <= h_sig_Eg->GetNbinsX(); ++ib_g) {
      const double Eg = h_sig_Eg->GetXaxis()->GetBinCenter(ib_g);
      const double dEg = h_sig_Eg->GetXaxis()->GetBinWidth(ib_g);
      for (int ib_t = 1; ib_t <= h_sig_t->GetNbinsX(); ++ib_t) {
        const double t = h_sig_t->GetXaxis()->GetBinCenter(ib_t);
        const double dt = h_sig_t->GetXaxis()->GetBinWidth(ib_t);
        for (int ib_pe = 1; ib_pe <= h_sig_phi_e->GetNbinsX(); ++ib_pe) {
          const double phi_e = h_sig_phi_e->GetXaxis()->GetBinCenter(ib_pe);
          const double dphi_e = h_sig_phi_e->GetXaxis()->GetBinWidth(ib_pe);
          for (int ib_pg = 1; ib_pg <= h_sig_phi_g->GetNbinsX(); ++ib_pg) {
            const double phi_g = h_sig_phi_g->GetXaxis()->GetBinCenter(ib_pg);
            const double dphi_g = h_sig_phi_g->GetXaxis()->GetBinWidth(ib_pg);
            const double dV = dEe * dEg * dt * dphi_e * dphi_g;
            const double theta = std::fabs(
                Detector_PhiSnapToGrid(phi_e, detres.phi_e_min, detres.phi_e_max, N_phi_e) -
                Detector_PhiSnapToGrid(phi_g, detres.phi_g_min, detres.phi_g_max, N_phi_g));

            const double p_sig = SignalPdf(Ee, Eg, t, phi_e, phi_g,
                                           analysis_window, detres, kMassesPDG);
            if (p_sig > 0.0 && IsFinite(p_sig)) {
              const double mass = p_sig * dV;
              sig_norm += mass;
              ++n_sig_pos;
              h_sig_Ee->Fill(Ee, mass);
              h_sig_Eg->Fill(Eg, mass);
              h_sig_t->Fill(t, mass);
              h_sig_phi_e->Fill(phi_e, mass);
              h_sig_phi_g->Fill(phi_g, mass);
              h_sig_theta->Fill(theta, mass);
              h2_sig_EeEg->Fill(Ee, Eg, mass);
              h2_sig_phi->Fill(phi_e, phi_g, mass);
            }

            const double p_rmd = RMDGridPdf(Ee, Eg, t, phi_e, phi_g);
            if (p_rmd > 0.0 && IsFinite(p_rmd)) {
              const double mass = p_rmd * dV;
              rmd_norm += mass;
              ++n_rmd_pos;
              h_rmd_Ee->Fill(Ee, mass);
              h_rmd_Eg->Fill(Eg, mass);
              h_rmd_t->Fill(t, mass);
              h_rmd_phi_e->Fill(phi_e, mass);
              h_rmd_phi_g->Fill(phi_g, mass);
              h_rmd_theta->Fill(theta, mass);
              h2_rmd_EeEg->Fill(Ee, Eg, mass);
              h2_rmd_phi->Fill(phi_e, phi_g, mass);
            }

            ++n_cells;
          }
        }
      }
    }
  }

  ConvertHist1ToDensity(h_sig_Ee);
  ConvertHist1ToDensity(h_rmd_Ee);
  ConvertHist1ToDensity(h_sig_Eg);
  ConvertHist1ToDensity(h_rmd_Eg);
  ConvertHist1ToDensity(h_sig_t);
  ConvertHist1ToDensity(h_rmd_t);
  ConvertHist1ToDensity(h_sig_phi_e);
  ConvertHist1ToDensity(h_rmd_phi_e);
  ConvertHist1ToDensity(h_sig_phi_g);
  ConvertHist1ToDensity(h_rmd_phi_g);
  ConvertHist1ToDensity(h_sig_theta);
  ConvertHist1ToDensity(h_rmd_theta);

  ConvertHist2ToDensity(h2_sig_EeEg);
  ConvertHist2ToDensity(h2_rmd_EeEg);
  ConvertHist2ToDensity(h2_sig_phi);
  ConvertHist2ToDensity(h2_rmd_phi);

  TFile fout(kOutRoot, "RECREATE");
  WriteHistPair(fout, h_sig_Ee, h_rmd_Ee, nullptr, nullptr);
  WriteHistPair(fout, h_sig_Eg, h_rmd_Eg, nullptr, nullptr);
  WriteHistPair(fout, h_sig_t, h_rmd_t, nullptr, nullptr);
  WriteHistPair(fout, h_sig_phi_e, h_rmd_phi_e, nullptr, nullptr);
  WriteHistPair(fout, h_sig_phi_g, h_rmd_phi_g, nullptr, nullptr);
  WriteHistPair(fout, h_sig_theta, h_rmd_theta, h2_sig_EeEg, h2_rmd_EeEg);
  WriteHistPair(fout, nullptr, nullptr, h2_sig_phi, h2_rmd_phi);
  fout.Close();

  TLegend leg(0.62, 0.76, 0.88, 0.88);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.AddEntry(h_sig_Ee, "Signal PDF", "l");
  leg.AddEntry(h_rmd_Ee, "RMD PDF", "l");

  TCanvas c0("c0", "meta", 1200, 800);
  c0.cd();
  DrawMetaPage(sig_norm, rmd_norm, n_cells, n_sig_pos, n_rmd_pos);

  TCanvas c1("c1", "1d", 1400, 900);
  c1.Divide(3, 2);

  auto draw_overlay = [&c1, &leg](int ipad, TH1D* hs, TH1D* hr) {
    if (!hs || !hr) return;
    c1.cd(ipad);
    gPad->SetGrid();
    const double ymax = (hs->GetMaximum() > hr->GetMaximum()) ? hs->GetMaximum() : hr->GetMaximum();
    hs->SetMaximum((ymax > 0.0) ? 1.20 * ymax : 1.0);
    hs->Draw("hist");
    hr->Draw("hist same");
    leg.Draw();
  };

  draw_overlay(1, h_sig_Ee, h_rmd_Ee);
  draw_overlay(2, h_sig_Eg, h_rmd_Eg);
  draw_overlay(3, h_sig_t, h_rmd_t);
  draw_overlay(4, h_sig_phi_e, h_rmd_phi_e);
  draw_overlay(5, h_sig_phi_g, h_rmd_phi_g);
  draw_overlay(6, h_sig_theta, h_rmd_theta);

  auto draw_2d = [](TH2D* h, bool draw_text) {
    if (!h) return;
    gPad->SetGrid();
    gPad->SetRightMargin(0.14);
    if (draw_text) h->Draw("colz text");
    else h->Draw("colz");
  };

  TCanvas c2("c2", "signal2d", 1200, 800);
  c2.Divide(2, 1);
  c2.cd(1); draw_2d(h2_sig_EeEg, false);
  c2.cd(2); draw_2d(h2_sig_phi, true);

  TCanvas c3("c3", "rmd2d", 1200, 800);
  c3.Divide(2, 1);
  c3.cd(1); draw_2d(h2_rmd_EeEg, false);
  c3.cd(2); draw_2d(h2_rmd_phi, true);

  c0.Print(Form("%s[", kOutPdf));
  c0.Print(kOutPdf);
  c1.Print(kOutPdf);
  c2.Print(kOutPdf);
  c3.Print(kOutPdf);
  c3.Print(Form("%s]", kOutPdf));

  c0.SaveAs(kOutPngMeta);
  c1.SaveAs(kOutPng1D);
  c2.SaveAs(kOutPngSig2D);
  c3.SaveAs(kOutPngRmd2D);

  std::cout << "[plot_signal_rmd_pdfs] wrote: " << kOutPdf << "\n";
  std::cout << "[plot_signal_rmd_pdfs] wrote: " << kOutRoot << "\n";
  std::cout << "[plot_signal_rmd_pdfs] signal norm = " << sig_norm << "\n";
  std::cout << "[plot_signal_rmd_pdfs] RMD norm    = " << rmd_norm << "\n";

  return 0;
}
