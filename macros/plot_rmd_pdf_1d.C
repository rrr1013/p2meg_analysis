#include <cmath>
#include <iostream>
#include <vector>

#include "TFile.h"
#include "TObject.h"
#include "TAxis.h"
#include "TH1D.h"
#include "THn.h"
#include "THnBase.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"

// ------------------------------------------------------------
// t の解析的 PDF（窓内正規化）
// DetectorResolution.h / AnalysisWindow.h と一致させること
// ------------------------------------------------------------
static const double t_min   = -2.0;  // [ns]
static const double t_max   =  2.0;  // [ns]
static const double t_mean  =  0.0;  // [ns]
static const double sigma_t =  0.5;  // [ns]

static double PtNormFactor_At(double tmin, double tmax, double mean, double sigma) {
  const double s2 = std::sqrt(2.0) * sigma;
  const double u1 = (tmax - mean) / s2;
  const double u0 = (tmin - mean) / s2;
  return 0.5 * (std::erf(u1) - std::erf(u0));
}

static double PtWindowNormalized(double t, double tmin, double tmax, double mean, double sigma) {
  if (!(sigma > 0.0)) return 0.0;
  const double pi = 3.14159265358979323846;
  const double z  = (t - mean) / sigma;
  const double g  = (1.0 / (std::sqrt(2.0 * pi) * sigma)) * std::exp(-0.5 * z * z);
  const double At = PtNormFactor_At(tmin, tmax, mean, sigma);
  if (!(At > 0.0)) return 0.0;
  return g / At;
}

// ------------------------------------------------------------
// THnBase の「ビン内容取得」を確実にする
// - THnBase::GetBin(const Int_t*) は Long64_t を返す
// - THnBase::GetBinContent(Long64_t) で読む
// ------------------------------------------------------------
static double GetHnBinContent(const THnBase& h, const std::vector<int>& idx) {
  const Long64_t b = h.GetBin(idx.data());
  return h.GetBinContent(b);
}

// ------------------------------------------------------------
// 3D格子 p3(Ee,Eg,theta) から 1D 周辺分布を作る
// （p3 は密度 [MeV^-2 rad^-1] を想定）
// ------------------------------------------------------------
static TH1D* MarginalTheta(const THnBase& h3) {
  const TAxis* axEe = h3.GetAxis(0);
  const TAxis* axEg = h3.GetAxis(1);
  const TAxis* axTh = h3.GetAxis(2);

  TH1D* hout = new TH1D("h_theta",
                        "RMD marginal p(#theta);#theta [rad];p(#theta) [rad^{-1}]",
                        axTh->GetNbins(), axTh->GetXmin(), axTh->GetXmax());
  hout->SetDirectory(nullptr);

  std::vector<int> idx(3, 1);
  for (int iEe = 1; iEe <= axEe->GetNbins(); ++iEe) {
    idx[0] = iEe;
    const double dEe = axEe->GetBinWidth(iEe);
    for (int iEg = 1; iEg <= axEg->GetNbins(); ++iEg) {
      idx[1] = iEg;
      const double dEg = axEg->GetBinWidth(iEg);
      for (int iTh = 1; iTh <= axTh->GetNbins(); ++iTh) {
        idx[2] = iTh;
        const double p3 = GetHnBinContent(h3, idx);
        if (p3 > 0.0 && std::isfinite(p3)) {
          hout->AddBinContent(iTh, p3 * dEe * dEg);
        }
      }
    }
  }
  return hout;
}

static TH1D* MarginalEe(const THnBase& h3) {
  const TAxis* axEe = h3.GetAxis(0);
  const TAxis* axEg = h3.GetAxis(1);
  const TAxis* axTh = h3.GetAxis(2);

  TH1D* hout = new TH1D("h_Ee",
                        "RMD marginal p(E_{e});E_{e} [MeV];p(E_{e}) [MeV^{-1}]",
                        axEe->GetNbins(), axEe->GetXmin(), axEe->GetXmax());
  hout->SetDirectory(nullptr);

  std::vector<int> idx(3, 1);
  for (int iEe = 1; iEe <= axEe->GetNbins(); ++iEe) {
    idx[0] = iEe;
    for (int iEg = 1; iEg <= axEg->GetNbins(); ++iEg) {
      idx[1] = iEg;
      const double dEg = axEg->GetBinWidth(iEg);
      for (int iTh = 1; iTh <= axTh->GetNbins(); ++iTh) {
        idx[2] = iTh;
        const double dTh = axTh->GetBinWidth(iTh);
        const double p3 = GetHnBinContent(h3, idx);
        if (p3 > 0.0 && std::isfinite(p3)) {
          hout->AddBinContent(iEe, p3 * dEg * dTh);
        }
      }
    }
  }
  return hout;
}

static TH1D* MarginalEg(const THnBase& h3) {
  const TAxis* axEe = h3.GetAxis(0);
  const TAxis* axEg = h3.GetAxis(1);
  const TAxis* axTh = h3.GetAxis(2);

  TH1D* hout = new TH1D("h_Eg",
                        "RMD marginal p(E_{#gamma});E_{#gamma} [MeV];p(E_{#gamma}) [MeV^{-1}]",
                        axEg->GetNbins(), axEg->GetXmin(), axEg->GetXmax());
  hout->SetDirectory(nullptr);

  std::vector<int> idx(3, 1);
  for (int iEe = 1; iEe <= axEe->GetNbins(); ++iEe) {
    idx[0] = iEe;
    const double dEe = axEe->GetBinWidth(iEe);
    for (int iEg = 1; iEg <= axEg->GetNbins(); ++iEg) {
      idx[1] = iEg;
      for (int iTh = 1; iTh <= axTh->GetNbins(); ++iTh) {
        idx[2] = iTh;
        const double dTh = axTh->GetBinWidth(iTh);
        const double p3 = GetHnBinContent(h3, idx);
        if (p3 > 0.0 && std::isfinite(p3)) {
          hout->AddBinContent(iEg, p3 * dEe * dTh);
        }
      }
    }
  }
  return hout;
}

static TH1D* MarginalT(int nbins = 200) {
  TH1D* ht = new TH1D("h_t",
                      "RMD marginal p(t);t [ns];p(t) [ns^{-1}]",
                      nbins, t_min, t_max);
  ht->SetDirectory(nullptr);

  for (int i = 1; i <= ht->GetNbinsX(); ++i) {
    const double t = ht->GetBinCenter(i);
    ht->SetBinContent(i, PtWindowNormalized(t, t_min, t_max, t_mean, sigma_t));
  }
  return ht;
}

// ------------------------------------------------------------
// メイン：4つの 1D 周辺分布を描画
// ------------------------------------------------------------
void plot_rmd_pdf_1d(const char* filepath="data/pdf_cache/rmd_grid.root",
                     const char* key="rmd_grid") {
  gROOT->cd();
  gStyle->SetOptStat(0);

  TFile f(filepath, "READ");
  if (f.IsZombie()) {
    std::cerr << "[plot] cannot open: " << filepath << "\n";
    return;
  }

  TObject* obj = f.Get(key);
  if (!obj) {
    std::cerr << "[plot] key not found: " << key << "\n";
    return;
  }

  THnBase* h3 = dynamic_cast<THnBase*>(obj);
  if (!h3) {
    std::cerr << "[plot] object is not THnBase: " << key << "\n";
    return;
  }
  if (h3->GetNdimensions() != 3) {
    std::cerr << "[plot] grid dimension is not 3: ndims=" << h3->GetNdimensions() << "\n";
    return;
  }

  TH1D* hEe = MarginalEe(*h3);
  TH1D* hEg = MarginalEg(*h3);
  TH1D* hTh = MarginalTheta(*h3);
  TH1D* ht  = MarginalT(200);

  std::cout << "Integral p(Ee) dEe      = " << hEe->Integral("width") << "\n";
  std::cout << "Integral p(Eg) dEg      = " << hEg->Integral("width") << "\n";
  std::cout << "Integral p(theta) dth   = " << hTh->Integral("width") << "\n";
  std::cout << "Integral p(t) dt        = " << ht ->Integral("width") << "\n";

  std::cout << "Max p(Ee)    = " << hEe->GetMaximum() << "\n";
  std::cout << "Max p(Eg)    = " << hEg->GetMaximum() << "\n";
  std::cout << "Max p(theta) = " << hTh->GetMaximum() << "\n";
  std::cout << "Max p(t)     = " << ht ->GetMaximum() << "\n";

  // 最大値がゼロなら「投影が空」
  if (hEe->GetMaximum() <= 0 && hEg->GetMaximum() <= 0 && hTh->GetMaximum() <= 0) {
    std::cerr << "[plot] projections look empty (max==0). Something is wrong in bin reading.\n";
    return;
  }

  TCanvas* c = new TCanvas("c_rmd_1d", "RMD 1D marginals", 1200, 900);
  c->Divide(2,2);

  c->cd(1);
  hEe->SetMinimum(0.0);
  hEe->SetMaximum(hEe->GetMaximum() * 1.2);
  hEe->Draw("hist");
  gPad->Modified(); gPad->Update();

  c->cd(2);
  hEg->SetMinimum(0.0);
  hEg->SetMaximum(hEg->GetMaximum() * 1.2);
  hEg->Draw("hist");
  gPad->Modified(); gPad->Update();

  c->cd(3);
  hTh->SetMinimum(0.0);
  hTh->SetMaximum(hTh->GetMaximum() * 1.2);
  hTh->Draw("hist");
  gPad->Modified(); gPad->Update();

  c->cd(4);
  ht->SetMinimum(0.0);
  ht->SetMaximum(ht->GetMaximum() * 1.2);
  ht->Draw("hist");
  gPad->Modified(); gPad->Update();

  c->Modified();
  c->Update();
  f.Close();
}
