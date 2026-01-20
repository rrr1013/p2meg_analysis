// src/RMDGridPdf.cc
#include "p2meg/RMDGridPdf.h"

#include <cmath>
#include <iostream>
#include <string>
#include <limits>

#include "TFile.h"
#include "THn.h"
#include "TAxis.h"

#include "p2meg/AnalysisWindow.h"
#include "p2meg/AnalysisWindowUtils.h"
#include "p2meg/AngleUtils.h"
#include "p2meg/DetectorResolution.h"
#include "p2meg/HistUtils.h"
#include "p2meg/MathUtils.h"

//============================================================
// 内部状態
//============================================================

// 4D 格子
static THnD* gHist = nullptr;

//============================================================
// 窓内で正規化された時間ガウシアン（密度）
// pt(t) = N(t_mean, sigma_t) / A_t, ただし A_t = ∫_{tmin}^{tmax} N dt
//============================================================
static double PtWindowNormalized(double t) {
  const double sigma_t = detres.sigma_t;
  const double t_mean  = detres.t_mean;

  const double tmin = analysis_window.t_min;
  const double tmax = analysis_window.t_max;

  if (!(sigma_t > 0.0) || !Math_IsFinite(sigma_t)) return 0.0;

  const double inv_s = 1.0 / sigma_t;
  const double z = (t - t_mean) * inv_s;

  const double norm = 1.0 / (std::sqrt(2.0 * pi) * sigma_t);
  const double g = norm * std::exp(-0.5 * z * z);

  const double s2 = std::sqrt(2.0) * sigma_t;
  const double u1 = (tmax - t_mean) / s2;
  const double u0 = (tmin - t_mean) / s2;
  const double At = 0.5 * (std::erf(u1) - std::erf(u0));

  if (!(At > 0.0) || !Math_IsFinite(At)) return 0.0;
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

  if (!Math_IsFinite(Ee) || !Math_IsFinite(Eg) || !Math_IsFinite(t) ||
      !Math_IsFinite(phi_detector_e) || !Math_IsFinite(phi_detector_g)) {
    return 0.0;
  }

  // theta_eg = |phi_e - phi_g|
  const double theta_eg = Angle_ThetaFromPhiClipped(phi_detector_e, phi_detector_g);

  // 解析窓（4D: Ee,Eg,t,theta）でカット
  if (!AnalysisWindow_In4D(analysis_window, Ee, Eg, t, theta_eg)) return 0.0;

  // phi は軸範囲にクリップ（x==xmax は overflow になり得るので、必ず xmax の内側へ落とす）
  const TAxis* axPe = gHist->GetAxis(2);
  const TAxis* axPg = gHist->GetAxis(3);
  const double pe_hi = std::nextafter(axPe->GetXmax(), 0.0);
  const double pg_hi = std::nextafter(axPg->GetXmax(), 0.0);
  const double pe_in = Math_Clamp(phi_detector_e, 0.0, pe_hi);
  const double pg_in = Math_Clamp(phi_detector_g, 0.0, pg_hi);

  // phi 軸は可変ビンなので FindBin で最近傍に落ちる
  const int bin_pe = axPe->FindBin(pe_in);
  const int bin_pg = axPg->FindBin(pg_in);

  // 4D格子から p4(Ee,Eg,phi_e,phi_g) を評価（Ee,Eg は補間、phi は固定ビン）
  const double p4 = Hist_InterpEeEg4(*gHist, Ee, Eg, bin_pe, bin_pg);
  if (!(p4 > 0.0) || !Math_IsFinite(p4)) return 0.0;

  // 時間因子を解析的に掛ける（窓内正規化）
  const double pt = PtWindowNormalized(t);
  if (!(pt > 0.0) || !Math_IsFinite(pt)) return 0.0;

  const double pdf = p4 * pt;
  if (!(pdf > 0.0) || !Math_IsFinite(pdf)) return 0.0;
  return pdf;
}
