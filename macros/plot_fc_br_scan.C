#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "TBox.h"
#include "TCanvas.h"
#include "TEfficiency.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TPad.h"
#include "TStyle.h"
#include "TH1.h"

struct FCBRPoint {
  double br;
  double n_sig_nom;
  double q_obs;
  double p_value;
  int accepted;
  int valid_toys;
  double n_rmd_prof;
  double n_acc_prof;
};

static bool ParseFCBRPoint(const std::string& line, FCBRPoint& out) {
  if (line.empty()) return false;
  if (line[0] == '#') return false;

  std::istringstream iss(line);
  if (!(iss >> out.br
            >> out.n_sig_nom
            >> out.q_obs
            >> out.p_value
            >> out.accepted
            >> out.valid_toys
            >> out.n_rmd_prof
            >> out.n_acc_prof)) {
    return false;
  }
  return true;
}

static bool LoadFCBRPoints(const char* input_txt, std::vector<FCBRPoint>& points) {
  points.clear();
  if (!input_txt || !input_txt[0]) return false;

  std::ifstream fin(input_txt);
  if (!fin) return false;

  std::string line;
  bool in_scan_block = false;
  while (std::getline(fin, line)) {
    if (line.find("# BR") != std::string::npos) {
      in_scan_block = true;
      continue;
    }
    if (!in_scan_block) continue;
    if (line.find("====") != std::string::npos) break;

    FCBRPoint p{};
    if (ParseFCBRPoint(line, p)) points.push_back(p);
  }
  return !points.empty();
}

static bool ClopperPearsonInterval(const FCBRPoint& p, double& ylo, double& yhi) {
  ylo = p.p_value;
  yhi = p.p_value;
  if (p.valid_toys <= 0) return false;
  if (!(p.p_value >= 0.0) || !(p.p_value <= 1.0) || !std::isfinite(p.p_value)) return false;

  const int n = p.valid_toys;
  int k = static_cast<int>(std::llround(p.p_value * static_cast<double>(n)));
  if (k < 0) k = 0;
  if (k > n) k = n;

  const double alpha = 1.0 - 0.682689492; // 約 1 sigma
  ylo = TEfficiency::ClopperPearson(n, k, alpha * 0.5, false);
  yhi = TEfficiency::ClopperPearson(n, k, alpha * 0.5, true);
  return true;
}

static int RepresentativeToyCount(const std::vector<FCBRPoint>& points) {
  if (points.empty()) return 0;

  std::vector<int> vals;
  vals.reserve(points.size());
  for (const auto& p : points) {
    if (p.valid_toys > 0) vals.push_back(p.valid_toys);
  }
  if (vals.empty()) return 0;

  std::sort(vals.begin(), vals.end());
  return vals[vals.size() / 2];
}

static TGraph* BuildQGraph(const std::vector<FCBRPoint>& points,
                          double x_shift = 0.0) {
  const int n = static_cast<int>(points.size());
  TGraph* gr = new TGraph(n);
  for (int i = 0; i < n; ++i) {
    gr->SetPoint(i, points[i].br + x_shift, points[i].q_obs);
  }
  return gr;
}

static TGraphAsymmErrors* BuildPGraphErrors(const std::vector<FCBRPoint>& points,
                                            double x_shift = 0.0,
                                            double xmin = -1.0,
                                            double xmax = -1.0) {
  std::vector<FCBRPoint> kept;
  kept.reserve(points.size());
  for (const auto& p : points) {
    if (xmin < xmax && (p.br < xmin || p.br > xmax)) continue;
    kept.push_back(p);
  }

  const int n = static_cast<int>(kept.size());
  TGraphAsymmErrors* gr = new TGraphAsymmErrors(n);
  for (int i = 0; i < n; ++i) {
    double ylo = kept[i].p_value;
    double yhi = kept[i].p_value;
    ClopperPearsonInterval(kept[i], ylo, yhi);
    gr->SetPoint(i, kept[i].br + x_shift, kept[i].p_value);
    gr->SetPointError(i, 0.0, 0.0,
                      kept[i].p_value - ylo,
                      yhi - kept[i].p_value);
  }
  return gr;
}

static void StylePGraph(TGraphAsymmErrors* gr, int color, int marker_style) {
  if (!gr) return;
  gr->SetLineWidth(2);
  gr->SetLineColor(color);
  gr->SetMarkerColor(color);
  gr->SetMarkerStyle(marker_style);
  gr->SetMarkerSize(0.95);
}

void plot_fc_br_scan(const char* input_wide_txt,
                     const char* output_pdf,
                     double br90 = -1.0,
                     const char* input_focus_txt = "") {
  std::vector<FCBRPoint> wide_points;
  if (!LoadFCBRPoints(input_wide_txt, wide_points)) {
    Error("plot_fc_br_scan", "cannot load wide scan from %s", input_wide_txt);
    return;
  }

  std::vector<FCBRPoint> focus_points;
  const bool has_focus = LoadFCBRPoints(input_focus_txt, focus_points);

  gStyle->SetOptStat(0);

  double br_min = wide_points.front().br;
  double br_max = wide_points.back().br;
  double q_max = 0.0;
  for (const auto& p : wide_points) q_max = std::max(q_max, p.q_obs);
  for (const auto& p : focus_points) q_max = std::max(q_max, p.q_obs);
  if (!(q_max > 0.0)) q_max = 1.0;

  double zoom_min = br_min;
  double zoom_max = br_max;
  double zoom_p_min = 1.0;
  double zoom_p_max = 0.0;
  if (has_focus) {
    zoom_min = focus_points.front().br;
    zoom_max = focus_points.back().br;
    for (const auto& p : focus_points) {
      double ylo = p.p_value;
      double yhi = p.p_value;
      ClopperPearsonInterval(p, ylo, yhi);
      zoom_p_min = std::min(zoom_p_min, ylo);
      zoom_p_max = std::max(zoom_p_max, yhi);
    }
  }
  if (!(zoom_p_max > zoom_p_min)) {
    zoom_p_min = 0.08;
    zoom_p_max = 0.12;
  }
  zoom_p_min = std::max(0.0, zoom_p_min - 0.01);
  zoom_p_max = std::min(1.0, zoom_p_max + 0.01);

  const double x_shift_wide = 0.0;
  const double x_shift_focus = 0.0;

  TGraph* gr_q_wide = BuildQGraph(wide_points, x_shift_wide);
  TGraphAsymmErrors* gr_p_wide = BuildPGraphErrors(wide_points, x_shift_wide);
  TGraph* gr_q_focus = has_focus ? BuildQGraph(focus_points, x_shift_focus) : nullptr;
  TGraphAsymmErrors* gr_p_focus =
      has_focus ? BuildPGraphErrors(focus_points, x_shift_focus) : nullptr;
  TGraphAsymmErrors* gr_p_focus_zoom =
      has_focus ? BuildPGraphErrors(focus_points, x_shift_focus, zoom_min, zoom_max) : nullptr;
  TGraphAsymmErrors* gr_p_wide_zoom =
      has_focus ? BuildPGraphErrors(wide_points, x_shift_wide, zoom_min, zoom_max) : nullptr;
  const int wide_toys = RepresentativeToyCount(wide_points);
  const int focus_toys = RepresentativeToyCount(focus_points);

  TCanvas* c = new TCanvas("c_fc_br", "FC BR scan", 900, 820);
  c->Divide(1, 2);

  c->cd(1);
  TH1F* hframe_q = new TH1F("hframe_q", ";BR test;q_{obs}(BR)", 100, br_min, br_max);
  hframe_q->SetMinimum(0.0);
  hframe_q->SetMaximum(1.15 * q_max);
  hframe_q->Draw();

  gr_q_wide->SetLineWidth(2);
  gr_q_wide->SetLineColor(kGray + 2);
  gr_q_wide->SetMarkerColor(kGray + 2);
  gr_q_wide->SetMarkerStyle(20);
  gr_q_wide->SetMarkerSize(1.0);
  gr_q_wide->Draw("LP SAME");

  if (has_focus) {
    gr_q_focus->SetLineWidth(2);
    gr_q_focus->SetLineColor(kRed + 1);
    gr_q_focus->SetMarkerColor(kRed + 1);
    gr_q_focus->SetMarkerStyle(21);
    gr_q_focus->SetMarkerSize(1.0);
    gr_q_focus->Draw("LP SAME");

    TBox* zoom_box_q = new TBox(zoom_min, 0.0, zoom_max, 1.15 * q_max);
    zoom_box_q->SetFillStyle(0);
    zoom_box_q->SetLineStyle(3);
    zoom_box_q->SetLineColor(kRed + 1);
    zoom_box_q->Draw();
  }

  if (br90 > 0.0) {
    TLine* l90q = new TLine(br90, 0.0, br90, 1.15 * q_max);
    l90q->SetLineStyle(2);
    l90q->SetLineColor(kBlue + 2);
    l90q->Draw();
  }

  c->cd(2);
  TH1F* hframe_p = new TH1F("hframe_p", ";BR test;p-value", 100, br_min, br_max);
  hframe_p->SetMinimum(0.0);
  hframe_p->SetMaximum(1.05);
  hframe_p->Draw();

  StylePGraph(gr_p_wide, kGray + 2, 24);
  gr_p_wide->Draw("E SAME");
  gr_p_wide->Draw("P SAME");

  if (has_focus) {
    StylePGraph(gr_p_focus, kRed + 1, 25);
    gr_p_focus->Draw("E SAME");
    gr_p_focus->Draw("P SAME");

    TBox* zoom_box_p = new TBox(zoom_min, 0.0, zoom_max, 1.05);
    zoom_box_p->SetFillStyle(0);
    zoom_box_p->SetLineStyle(3);
    zoom_box_p->SetLineColor(kRed + 1);
    zoom_box_p->Draw();
  }

  TLine* lcl = new TLine(br_min, 0.1, br_max, 0.1);
  lcl->SetLineStyle(2);
  lcl->SetLineColor(kGreen + 2);
  lcl->Draw();

  if (br90 > 0.0) {
    TLine* l90p = new TLine(br90, 0.0, br90, 1.05);
    l90p->SetLineStyle(2);
    l90p->SetLineColor(kBlue + 2);
    l90p->Draw();
  }

  TLegend* leg = new TLegend(0.14, 0.73, 0.48, 0.89);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(gr_p_wide,
                Form("wide scan (%d toy/point)", wide_toys > 0 ? wide_toys : 0),
                "lep");
  if (has_focus) {
    leg->AddEntry(gr_p_focus,
                  Form("focused high-stat (%d toy/point)", focus_toys > 0 ? focus_toys : 0),
                  "lep");
  }
  leg->AddEntry(lcl, "90% C.L. threshold", "l");
  leg->Draw();

  if (has_focus) {
    TPad* inset = new TPad("p_zoom", "p_zoom", 0.56, 0.32, 0.95, 0.70);
    inset->SetFillStyle(1001);
    inset->SetFillColor(kWhite);
    inset->SetMargin(0.16, 0.06, 0.16, 0.08);
    inset->Draw();
    inset->cd();

    TH1F* hzoom = new TH1F("hframe_p_zoom", ";BR test;p-value", 100, zoom_min, zoom_max);
    hzoom->SetMinimum(zoom_p_min);
    hzoom->SetMaximum(zoom_p_max);
    hzoom->GetXaxis()->SetLabelSize(0.055);
    hzoom->GetYaxis()->SetLabelSize(0.055);
    hzoom->GetXaxis()->SetTitleSize(0.055);
    hzoom->GetYaxis()->SetTitleSize(0.055);
    hzoom->Draw();

    if (gr_p_wide_zoom) {
      StylePGraph(gr_p_wide_zoom, kGray + 2, 24);
      gr_p_wide_zoom->Draw("E SAME");
      gr_p_wide_zoom->Draw("P SAME");
    }
    if (gr_p_focus_zoom) {
      StylePGraph(gr_p_focus_zoom, kRed + 1, 25);
      gr_p_focus_zoom->Draw("E SAME");
      gr_p_focus_zoom->Draw("P SAME");
    }

    TLine* lcl_zoom = new TLine(zoom_min, 0.1, zoom_max, 0.1);
    lcl_zoom->SetLineStyle(2);
    lcl_zoom->SetLineColor(kGreen + 2);
    lcl_zoom->Draw();

    if (br90 > 0.0) {
      TLine* l90p_zoom = new TLine(br90, zoom_p_min, br90, zoom_p_max);
      l90p_zoom->SetLineStyle(2);
      l90p_zoom->SetLineColor(kBlue + 2);
      l90p_zoom->Draw();
    }
  }

  c->cd(2);
  if (br90 > 0.0) {
    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.040);
    lat.DrawLatex(0.52, 0.78, Form("BR_{90} = %.6g", br90));
  }

  c->Print(output_pdf);
}
