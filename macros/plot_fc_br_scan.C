#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TGraph.h"
#include "TLine.h"
#include "TLatex.h"
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

void plot_fc_br_scan(const char* input_txt,
                     const char* output_pdf,
                     double br90 = -1.0) {
  std::ifstream fin(input_txt);
  if (!fin) {
    Error("plot_fc_br_scan", "cannot open %s", input_txt);
    return;
  }

  std::vector<FCBRPoint> points;
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

  if (points.empty()) {
    Error("plot_fc_br_scan", "no scan points found in %s", input_txt);
    return;
  }

  gStyle->SetOptStat(0);

  const int n = static_cast<int>(points.size());
  TGraph* gr_q = new TGraph(n);
  TGraph* gr_p = new TGraph(n);

  double br_min = points.front().br;
  double br_max = points.back().br;
  double q_max = 0.0;
  for (int i = 0; i < n; ++i) {
    gr_q->SetPoint(i, points[i].br, points[i].q_obs);
    gr_p->SetPoint(i, points[i].br, points[i].p_value);
    if (points[i].q_obs > q_max) q_max = points[i].q_obs;
  }
  if (!(q_max > 0.0)) q_max = 1.0;

  TCanvas* c = new TCanvas("c_fc_br", "FC BR scan", 900, 800);
  c->Divide(1, 2);

  c->cd(1);
  TH1F* hframe_q = new TH1F("hframe_q", ";BR test;q_{obs}(BR)", 100, br_min, br_max);
  hframe_q->SetMinimum(0.0);
  hframe_q->SetMaximum(1.15 * q_max);
  hframe_q->Draw();
  gr_q->SetLineWidth(2);
  gr_q->SetMarkerStyle(20);
  gr_q->SetMarkerSize(1.0);
  gr_q->Draw("LP SAME");
  if (br90 > 0.0) {
    TLine* l90q = new TLine(br90, 0.0, br90, 1.15 * q_max);
    l90q->SetLineStyle(2);
    l90q->SetLineColor(kRed + 1);
    l90q->Draw();
  }

  c->cd(2);
  TH1F* hframe_p = new TH1F("hframe_p", ";BR test;p-value", 100, br_min, br_max);
  hframe_p->SetMinimum(0.0);
  hframe_p->SetMaximum(1.05);
  hframe_p->Draw();
  gr_p->SetLineWidth(2);
  gr_p->SetMarkerStyle(20);
  gr_p->SetMarkerSize(1.0);
  gr_p->SetLineColor(kBlue + 1);
  gr_p->SetMarkerColor(kBlue + 1);
  gr_p->Draw("LP SAME");

  TLine* lcl = new TLine(br_min, 0.1, br_max, 0.1);
  lcl->SetLineStyle(2);
  lcl->SetLineColor(kGreen + 2);
  lcl->Draw();

  if (br90 > 0.0) {
    TLine* l90p = new TLine(br90, 0.0, br90, 1.05);
    l90p->SetLineStyle(2);
    l90p->SetLineColor(kRed + 1);
    l90p->Draw();

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.04);
    lat.DrawLatex(0.16, 0.86, Form("BR_{90} = %.6g", br90));
  }

  c->Print(output_pdf);
}
