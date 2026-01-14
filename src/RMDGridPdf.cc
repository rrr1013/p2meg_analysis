// src/RMDGridPdf.cc
#include "p2meg/RMDGridPdf.h"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "THn.h"
#include "TAxis.h"

#include "p2meg/AnalysisWindow.h"
#include "p2meg/DetectorResolution.h"
#include "p2meg/Constants.h"

//============================================================
// 内部状態
//============================================================

// 4D 格子
static THnD* gHist = nullptr;

static bool IsFinite(double x) { return std::isfinite(x); }

static double Clamp(double x, double lo, double hi) {
  if (x < lo) return lo;
  if (x > hi) return hi;
  return x;
}

static bool IsInsideWindow4D(double Ee, double Eg, double t, double theta_eg) {
  if (Ee < analysis_window.Ee_min || Ee > analysis_window.Ee_max) return false;
  if (Eg < analysis_window.Eg_min || Eg > analysis_window.Eg_max) return false;
  if (t  < analysis_window.t_min  || t  > analysis_window.t_max ) return false;
  if (theta_eg < analysis_window.theta_min || theta_eg > analysis_window.theta_max) return false;
  return true;
}

// θ を 0..π に収めた上で |Δφ| を作る（0..π を想定）
static double ThetaFromPhi(double phi_e, double phi_g) {
  if (!IsFinite(phi_e) || !IsFinite(phi_g)) return 0.0;
  const double pe = Clamp(phi_e, 0.0, pi);
  const double pg = Clamp(phi_g, 0.0, pi);
  return std::fabs(pe - pg);
}

// 軸が等間隔ビンである前提で、座標 x を隣り合う2ビンに落とす（Ee,Eg 用）。
// 値は「ビン中心に定義された格子値」と見做して補間するための (i0, i1, f) を返す。
//  - i0, i1 は [1..nbins] のビン番号
//  - f は 0..1 の補間係数（i0 側が (1-f)、i1 側が f）
static int AxisBracketUniform(const TAxis& ax, double x, int& i0, int& i1, double& f) {
  const int n = ax.GetNbins();
  if (n < 2) return 1;

  const double xmin = ax.GetXmin();
  const double xmax = ax.GetXmax();
  const double dx = (xmax - xmin) / n;
  if (!(dx > 0.0) || !IsFinite(dx)) return 2;

  // ビン中心基準の連続座標（1番ビン中心が 1.0）
  double i_float = (x - xmin) / dx + 0.5;

  // i0 と i1=i0+1 が必要なので i0 は 1..n-1
  if (i_float < 1.0) i_float = 1.0;
  const double eps = 1e-12;
  const double max_ifloat = static_cast<double>(n) - eps; // < n
  if (i_float > max_ifloat) i_float = max_ifloat;

  i0 = static_cast<int>(std::floor(i_float));
  if (i0 < 1) i0 = 1;
  if (i0 > n - 1) i0 = n - 1;

  i1 = i0 + 1;
  f = i_float - static_cast<double>(i0);
  if (f < 0.0) f = 0.0;
  if (f > 1.0) f = 1.0;

  return 0;
}

// Ee,Eg を 2D（4点）補間して、phi軸は固定ビン（最近傍）で評価する。
static double InterpEeEg_4(const THnD& h, double Ee, double Eg, int bin_phi_e, int bin_phi_g) {
  const TAxis* axE = h.GetAxis(0);
  const TAxis* axG = h.GetAxis(1);

  int i0, i1, j0, j1;
  double fE, fG;

  if (AxisBracketUniform(*axE, Ee, i0, i1, fE) != 0) return 0.0;
  if (AxisBracketUniform(*axG, Eg, j0, j1, fG) != 0) return 0.0;

  std::vector<int> idx(4, 1);

  auto get = [&](int ie, int ig) -> double {
    idx[0] = ie;
    idx[1] = ig;
    idx[2] = bin_phi_e;
    idx[3] = bin_phi_g;
    const Long64_t bin = h.GetBin(idx.data());
    const double v = h.GetBinContent(bin);
    if (v > 0.0 && IsFinite(v)) return v;
    return 0.0;
  };

  const double v00 = get(i0, j0);
  const double v10 = get(i1, j0);
  const double v01 = get(i0, j1);
  const double v11 = get(i1, j1);

  const double a0 = (1.0 - fE) * v00 + fE * v10;
  const double a1 = (1.0 - fE) * v01 + fE * v11;
  const double v  = (1.0 - fG) * a0  + fG * a1;

  return (v > 0.0 && IsFinite(v)) ? v : 0.0;
}

// 窓内で正規化された時間ガウシアン（密度）
// pt(t) = N(t_mean, sigma_t) / A_t, ただし A_t = ∫_{tmin}^{tmax} N dt
static double PtWindowNormalized(double t) {
  const double sigma_t = detres.sigma_t;
  const double t_mean  = detres.t_mean;

  const double tmin = analysis_window.t_min;
  const double tmax = analysis_window.t_max;

  if (!(sigma_t > 0.0) || !IsFinite(sigma_t)) return 0.0;

  const double inv_s = 1.0 / sigma_t;
  const double z = (t - t_mean) * inv_s;

  const double norm = 1.0 / (std::sqrt(2.0 * pi) * sigma_t);
  const double g = norm * std::exp(-0.5 * z * z);

  const double s2 = std::sqrt(2.0) * sigma_t;
  const double u1 = (tmax - t_mean) / s2;
  const double u0 = (tmin - t_mean) / s2;
  const double At = 0.5 * (std::erf(u1) - std::erf(u0));

  if (!(At > 0.0) || !IsFinite(At)) return 0.0;
  return g / At;
}

//============================================================
// 公開関数
//============================================================

bool RMDGridPdf_Load(const char* filepath, const char* key) {
  if (!filepath || !key) return false;

  if (gHist) { delete gHist; gHist = nullptr; }

  TFile f(filepath, "READ");
  if (f.IsZombie()) {
    std::cerr << "[RMDGridPdf_Load] cannot open: " << filepath << "\n";
    return false;
  }

  TObject* obj = f.Get(key);

  if (!obj) {
    std::cerr << "[RMDGridPdf_Load] key not found: " << key << "\n";
    f.Close();
    return false;
  }

  THnD* h = dynamic_cast<THnD*>(obj);

  if (!h) {
    std::cerr << "[RMDGridPdf_Load] object is not THnD: " << key << "\n";
    f.Close();
    return false;
  }

  if (h->GetNdimensions() != 4) {
    std::cerr << "[RMDGridPdf_Load] THnD dimension is not 4: "
              << "N=" << h->GetNdimensions() << "\n";
    f.Close();
    return false;
  }

  gHist = dynamic_cast<THnD*>(h->Clone());
  f.Close();

  if (!gHist) {
    std::cerr << "[RMDGridPdf_Load] clone failed\n";
    if (gHist) { delete gHist; gHist = nullptr; }
    return false;
  }

  return true;
}

bool RMDGridPdf_IsLoaded() {
  return (gHist != nullptr);
}

double RMDGridPdf(double Ee, double Eg, double t,
                  double phi_detector_e, double phi_detector_g) {
  if (!gHist) return 0.0;

  if (!IsFinite(Ee) || !IsFinite(Eg) || !IsFinite(t) ||
      !IsFinite(phi_detector_e) || !IsFinite(phi_detector_g)) {
    return 0.0;
  }

  // theta_eg = |phi_e - phi_g|
  const double theta_eg = ThetaFromPhi(phi_detector_e, phi_detector_g);

  // 解析窓（4D: Ee,Eg,t,theta）でカット
  if (!IsInsideWindow4D(Ee, Eg, t, theta_eg)) return 0.0;

  // phi は物理範囲にクリップ
  const double pe_in = Clamp(phi_detector_e, 0.0, pi);
  const double pg_in = Clamp(phi_detector_g, 0.0, pi);

  // phi 軸は「最近傍θ格子」になる可変ビンのはずなので FindBin で最近傍に落ちる
  const TAxis* axPe = gHist->GetAxis(2);
  const TAxis* axPg = gHist->GetAxis(3);
  const int bin_pe = axPe->FindBin(pe_in);
  const int bin_pg = axPg->FindBin(pg_in);

  // 4D格子から p4(Ee,Eg,phi_e,phi_g) を評価（Ee,Eg は補間、phi は固定ビン）
  const double p4 = InterpEeEg_4(*gHist, Ee, Eg, bin_pe, bin_pg);
  if (!(p4 > 0.0) || !IsFinite(p4)) return 0.0;

  // 時間因子を解析的に掛ける（窓内正規化）
  const double pt = PtWindowNormalized(t);
  if (!(pt > 0.0) || !IsFinite(pt)) return 0.0;

  const double pdf = p4 * pt;
  if (!(pdf > 0.0) || !IsFinite(pdf)) return 0.0;
  return pdf;
}
