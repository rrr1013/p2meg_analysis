#include <cmath>
#include <iostream>
#include <vector>

#include "TFile.h"
#include "TObject.h"
#include "TAxis.h"
#include "TH1D.h"
#include "THnBase.h"
#include "TCanvas.h"

//------------------------------------------------------------
// 解析窓・分解能（いまの設定に合わせて）
//  ※AnalysisWindow.h / DetectorResolution.h と一致しているか注意
//------------------------------------------------------------
static const double Ee_min = 40.0;     // [MeV]
static const double Ee_max = 55.0;     // [MeV]
static const double Eg_min = 40.0;     // [MeV]
static const double Eg_max = 55.0;     // [MeV]
static const double th_min = 2.2;      // [rad]
static const double th_max = 3.14159;  // [rad]

static const double t_min = -2.0;      // [ns]
static const double t_max =  2.0;      // [ns]
static const double sigma_t = 0.5;     // [ns]
static const double t_mean  = 0.0;     // [ns]

//------------------------------------------------------------
// pt の窓内正規化（連続式）
//------------------------------------------------------------
static double PtNormFactor_At(double tmin, double tmax, double mean, double sigma) {
  const double s2 = std::sqrt(2.0) * sigma;
  const double u1 = (tmax - mean) / s2;
  const double u0 = (tmin - mean) / s2;
  return 0.5 * (std::erf(u1) - std::erf(u0));
}

static double PtWindowNormalized(double t, double tmin, double tmax, double mean, double sigma) {
  if (!(sigma > 0.0)) return 0.0;
  const double pi = 3.14159265358979323846;
  const double z = (t - mean) / sigma;
  const double g = (1.0 / (std::sqrt(2.0 * pi) * sigma)) * std::exp(-0.5 * z * z);
  const double At = PtNormFactor_At(tmin, tmax, mean, sigma);
  if (!(At > 0.0)) return 0.0;
  return g / At;
}

//------------------------------------------------------------
// 3D格子の離散積分： Σ p3(bin) * ΔEe*ΔEg*Δθ
//------------------------------------------------------------
static double Integral3D(const THnBase& h) {
  const int ndim = h.GetNdimensions();
  if (ndim != 3) return 0.0;

  const TAxis* ax0 = h.GetAxis(0);
  const TAxis* ax1 = h.GetAxis(1);
  const TAxis* ax2 = h.GetAxis(2);

  const int n0 = ax0->GetNbins();
  const int n1 = ax1->GetNbins();
  const int n2 = ax2->GetNbins();

  std::vector<int> idx(3, 1);
  double sum = 0.0;

  for (int i0 = 1; i0 <= n0; ++i0) {
    idx[0] = i0;
    const double w0 = ax0->GetBinWidth(i0);
    for (int i1 = 1; i1 <= n1; ++i1) {
      idx[1] = i1;
      const double w1 = ax1->GetBinWidth(i1);
      for (int i2 = 1; i2 <= n2; ++i2) {
        idx[2] = i2;
        const double w2 = ax2->GetBinWidth(i2);

        const Long64_t bin = h.GetBin(idx.data());
        const double p3 = h.GetBinContent(bin); // 密度 [MeV^-2 rad^-1] のはず

        if (p3 > 0.0 && std::isfinite(p3)) {
          sum += p3 * (w0 * w1 * w2);
        }
      }
    }
  }
  return sum;
}

//------------------------------------------------------------
// 周辺分布（投影）を作る：Ee / Eg / theta
// ここでは「積分で周辺密度」になるように Δ(他軸) を掛けて足し込む
//------------------------------------------------------------
static TH1D* ProjectTheta(const THnBase& h) {
  const TAxis* axEe = h.GetAxis(0);
  const TAxis* axEg = h.GetAxis(1);
  const TAxis* axTh = h.GetAxis(2);

  auto* hout = new TH1D("h_theta_marginal", "marginal in #theta;#theta [rad];p(#theta) [rad^{-1}]",
                        axTh->GetNbins(), axTh->GetXmin(), axTh->GetXmax());

  std::vector<int> idx(3, 1);
  for (int iEe = 1; iEe <= axEe->GetNbins(); ++iEe) {
    idx[0] = iEe;
    const double dEe = axEe->GetBinWidth(iEe);
    for (int iEg = 1; iEg <= axEg->GetNbins(); ++iEg) {
      idx[1] = iEg;
      const double dEg = axEg->GetBinWidth(iEg);
      for (int iTh = 1; iTh <= axTh->GetNbins(); ++iTh) {
        idx[2] = iTh;
        const Long64_t bin = h.GetBin(idx.data());
        const double p3 = h.GetBinContent(bin);

        if (p3 > 0.0 && std::isfinite(p3)) {
          // p(theta) = ∫ p3 dEe dEg なので、ΔEeΔEg を掛けて theta bin に加算
          hout->AddBinContent(iTh, p3 * dEe * dEg);
        }
      }
    }
  }
  return hout;
}

static TH1D* ProjectEe(const THnBase& h) {
  const TAxis* axEe = h.GetAxis(0);
  const TAxis* axEg = h.GetAxis(1);
  const TAxis* axTh = h.GetAxis(2);

  auto* hout = new TH1D("h_Ee_marginal", "marginal in E_{e};E_{e} [MeV];p(E_{e}) [MeV^{-1}]",
                        axEe->GetNbins(), axEe->GetXmin(), axEe->GetXmax());

  std::vector<int> idx(3, 1);
  for (int iEe = 1; iEe <= axEe->GetNbins(); ++iEe) {
    idx[0] = iEe;
    for (int iEg = 1; iEg <= axEg->GetNbins(); ++iEg) {
      idx[1] = iEg;
      const double dEg = axEg->GetBinWidth(iEg);
      for (int iTh = 1; iTh <= axTh->GetNbins(); ++iTh) {
        idx[2] = iTh;
        const double dTh = axTh->GetBinWidth(iTh);

        const Long64_t bin = h.GetBin(idx.data());
        const double p3 = h.GetBinContent(bin);

        if (p3 > 0.0 && std::isfinite(p3)) {
          // p(Ee) = ∫ p3 dEg dθ なので、ΔEgΔθ を掛けて Ee bin に加算
          hout->AddBinContent(iEe, p3 * dEg * dTh);
        }
      }
    }
  }
  return hout;
}

static TH1D* ProjectEg(const THnBase& h) {
  const TAxis* axEe = h.GetAxis(0);
  const TAxis* axEg = h.GetAxis(1);
  const TAxis* axTh = h.GetAxis(2);

  auto* hout = new TH1D("h_Eg_marginal", "marginal in E_{#gamma};E_{#gamma} [MeV];p(E_{#gamma}) [MeV^{-1}]",
                        axEg->GetNbins(), axEg->GetXmin(), axEg->GetXmax());

  std::vector<int> idx(3, 1);
  for (int iEe = 1; iEe <= axEe->GetNbins(); ++iEe) {
    idx[0] = iEe;
    const double dEe = axEe->GetBinWidth(iEe);
    for (int iEg = 1; iEg <= axEg->GetNbins(); ++iEg) {
      idx[1] = iEg;
      for (int iTh = 1; iTh <= axTh->GetNbins(); ++iTh) {
        idx[2] = iTh;
        const double dTh = axTh->GetBinWidth(iTh);

        const Long64_t bin = h.GetBin(idx.data());
        const double p3 = h.GetBinContent(bin);

        if (p3 > 0.0 && std::isfinite(p3)) {
          // p(Eg) = ∫ p3 dEe dθ
          hout->AddBinContent(iEg, p3 * dEe * dTh);
        }
      }
    }
  }
  return hout;
}

//------------------------------------------------------------
// メイン
//------------------------------------------------------------
void check_rmd_grid_3d(const char* filepath="data/pdf_cache/rmd_grid.root",
                      const char* key="rmd_grid",
                      bool draw_projections=true) {
  TFile f(filepath, "READ");
  if (f.IsZombie()) {
    std::cerr << "[check] cannot open: " << filepath << "\n";
    return;
  }

  TObject* obj = f.Get(key);
  if (!obj) {
    std::cerr << "[check] key not found: " << key << "\n";
    return;
  }

  THnBase* h = dynamic_cast<THnBase*>(obj);
  if (!h) {
    std::cerr << "[check] object is not THnBase/THnD: " << key << "\n";
    return;
  }

  std::cout << "---- [1] grid info ----\n";
  std::cout << "ndim = " << h->GetNdimensions() << "\n";
  for (int d = 0; d < h->GetNdimensions(); ++d) {
    const TAxis* ax = h->GetAxis(d);
    std::cout << "axis " << d
              << " nbins=" << ax->GetNbins()
              << " range=[" << ax->GetXmin() << "," << ax->GetXmax() << "]\n";
  }

  std::cout << "\n---- [2] normalization of p3 ----\n";
  const double I3 = Integral3D(*h);
  std::cout << "Sum_bins p3*ΔEeΔEgΔθ = " << I3 << "\n";
  std::cout << "Expected ~ 1 (discretization error should be small)\n";

  std::cout << "\n---- [3] normalization of pt (analytic) ----\n";
  const double At = PtNormFactor_At(t_min, t_max, t_mean, sigma_t);
  std::cout << "At = ∫ N(t_mean,sigma_t) dt over window = " << At << "\n";
  // ざっくり数値積分でも確認（台形）
  const int N = 20000;
  double s = 0.0;
  for (int i = 0; i <= N; ++i) {
    const double t = t_min + (t_max - t_min) * (double(i) / N);
    const double w = (i==0 || i==N) ? 0.5 : 1.0;
    s += w * PtWindowNormalized(t, t_min, t_max, t_mean, sigma_t);
  }
  s *= (t_max - t_min) / N;
  std::cout << "Numeric ∫ pt dt over window = " << s << "\n";
  std::cout << "Expected ~ 1\n";

  std::cout << "\n---- [4] 4D normalization by separability ----\n";
  std::cout << "Since pdf4 = p3 * pt, integral over window should be I3 * 1 ~ " << I3 << "\n";

  if (draw_projections) {
    std::cout << "\n---- [5] projections ----\n";
    TH1D* hth = ProjectTheta(*h);
    TH1D* hEe = ProjectEe(*h);
    TH1D* hEg = ProjectEg(*h);

    // それぞれ積分が1になっているか（離散）
    std::cout << "Integral p(theta) dtheta ~ " << hth->Integral("width") << "\n";
    std::cout << "Integral p(Ee) dEe ~ " << hEe->Integral("width") << "\n";
    std::cout << "Integral p(Eg) dEg ~ " << hEg->Integral("width") << "\n";

    auto* c1 = new TCanvas("c1", "RMD grid projections", 1200, 800);
    c1->Divide(2,2);
    c1->cd(1); hth->Draw("hist");
    c1->cd(2); hEe->Draw("hist");
    c1->cd(3); hEg->Draw("hist");
    c1->cd(4);
    // pt の形も描く（参考）
    TH1D* ht = new TH1D("h_pt", "p_{t}(t) (window-normalized);t [ns];p_{t}(t) [ns^{-1}]",
                        200, t_min, t_max);
    for (int i = 1; i <= ht->GetNbinsX(); ++i) {
      const double t = ht->GetBinCenter(i);
      ht->SetBinContent(i, PtWindowNormalized(t, t_min, t_max, t_mean, sigma_t));
    }
    ht->Draw("hist");
  }
}
