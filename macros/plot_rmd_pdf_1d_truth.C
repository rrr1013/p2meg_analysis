#include <cmath>
#include <iostream>

#include "TCanvas.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TROOT.h"

// RMDSpectrum.cc で定義されている関数（ヘッダは使わず宣言だけする）
extern double RMD_d3B_dEe_dEg_dcos(double Ee, double Eg, double cosTheta, double d_min = 1e-6);


// ============================================================
// 理論RMD（スメアなし）1D周辺分布の可視化
//
// ・RMD_d3B_dEe_dEg_dcos を直接使用
// ・解析窓内で数値積分して 1D pdf を作る
// ・正規化は最後に width 付き積分で行う
// ============================================================

// ---- 解析窓（grid生成時と一致させる） ----
static const double Ee_min = 40.0;
static const double Ee_max = 55.0;
static const double Eg_min = 40.0;
static const double Eg_max = 55.0;
static const double th_min = 2.2;
static const double th_max = 3.141592653589793;

// ---- 数値積分の刻み（粗すぎなければOK） ----
static const int NEe_int = 120;
static const int NEg_int = 120;
static const int Nth_int = 120;

// ------------------------------------------------------------
// p(Ee) = ∫ dEg dθ  d3B/dEe dEg dcosθ
// ------------------------------------------------------------
TH1D* MakeTruth_pEe() {
  TH1D* h = new TH1D(
    "h_truth_Ee",
    "RMD truth marginal p(E_{e});E_{e} [MeV];p(E_{e}) [MeV^{-1}]",
    60, Ee_min, Ee_max
  );
  h->SetDirectory(nullptr);

  const double dEg = (Eg_max - Eg_min) / NEg_int;
  const double dth = (th_max - th_min) / Nth_int;

  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    const double Ee = h->GetBinCenter(i);
    double sum = 0.0;

    for (int j = 0; j < NEg_int; ++j) {
      const double Eg = Eg_min + (j + 0.5) * dEg;
      for (int k = 0; k < Nth_int; ++k) {
        const double th = th_min + (k + 0.5) * dth;
        const double c  = std::cos(th);

        const double w = RMD_d3B_dEe_dEg_dcos(Ee, Eg, c);
        if (w > 0.0) {
          sum += w * dEg * dth;
        }
      }
    }
    h->SetBinContent(i, sum);
  }

  h->Scale(1.0 / h->Integral("width"));
  return h;
}

// ------------------------------------------------------------
// p(Eg)
// ------------------------------------------------------------
TH1D* MakeTruth_pEg() {
  TH1D* h = new TH1D(
    "h_truth_Eg",
    "RMD truth marginal p(E_{#gamma});E_{#gamma} [MeV];p(E_{#gamma}) [MeV^{-1}]",
    60, Eg_min, Eg_max
  );
  h->SetDirectory(nullptr);

  const double dEe = (Ee_max - Ee_min) / NEe_int;
  const double dth = (th_max - th_min) / Nth_int;

  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    const double Eg = h->GetBinCenter(i);
    double sum = 0.0;

    for (int j = 0; j < NEe_int; ++j) {
      const double Ee = Ee_min + (j + 0.5) * dEe;
      for (int k = 0; k < Nth_int; ++k) {
        const double th = th_min + (k + 0.5) * dth;
        const double c  = std::cos(th);

        const double w = RMD_d3B_dEe_dEg_dcos(Ee, Eg, c);
        if (w > 0.0) {
          sum += w * dEe * dth;
        }
      }
    }
    h->SetBinContent(i, sum);
  }

  h->Scale(1.0 / h->Integral("width"));
  return h;
}

// ------------------------------------------------------------
// p(theta)
// ------------------------------------------------------------
TH1D* MakeTruth_pTheta() {
  TH1D* h = new TH1D(
    "h_truth_theta",
    "RMD truth marginal p(#theta);#theta [rad];p(#theta) [rad^{-1}]",
    60, th_min, th_max
  );
  h->SetDirectory(nullptr);

  const double dEe = (Ee_max - Ee_min) / NEe_int;
  const double dEg = (Eg_max - Eg_min) / NEg_int;

  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    const double th = h->GetBinCenter(i);
    const double c  = std::cos(th);
    double sum = 0.0;

    for (int j = 0; j < NEe_int; ++j) {
      const double Ee = Ee_min + (j + 0.5) * dEe;
      for (int k = 0; k < NEg_int; ++k) {
        const double Eg = Eg_min + (k + 0.5) * dEg;

        const double w = RMD_d3B_dEe_dEg_dcos(Ee, Eg, c);
        if (w > 0.0) {
          sum += w * dEe * dEg;
        }
      }
    }
    h->SetBinContent(i, sum);
  }

  h->Scale(1.0 / h->Integral("width"));
  return h;
}

// ------------------------------------------------------------
// メイン描画
// ------------------------------------------------------------
void plot_rmd_pdf_1d_truth() {
  gROOT->cd();
  gStyle->SetOptStat(0);

  TH1D* hEe = MakeTruth_pEe();
  TH1D* hEg = MakeTruth_pEg();
  TH1D* hTh = MakeTruth_pTheta();

  std::cout << "Integral p(Ee) = " << hEe->Integral("width") << std::endl;
  std::cout << "Integral p(Eg) = " << hEg->Integral("width") << std::endl;
  std::cout << "Integral p(th) = " << hTh->Integral("width") << std::endl;

  TCanvas* c = new TCanvas("c_truth_1d", "RMD truth 1D marginals", 1200, 900);
  c->Divide(2,2);

  c->cd(1); hEe->SetMinimum(0); hEe->Draw("hist");
  c->cd(2); hEg->SetMinimum(0); hEg->Draw("hist");
  c->cd(3); hTh->SetMinimum(0); hTh->Draw("hist");

  c->Update();
}
