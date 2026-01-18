// src/ACCGridPdf.cc
#include "p2meg/ACCGridPdf.h"

#include <cmath>
#include <iostream>
#include <string>

#include "TFile.h"
#include "TAxis.h"
#include "THn.h"
#include "TParameter.h"

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
// 公開関数
//============================================================

bool ACCGridPdf_Load(const char* filepath, const char* key) {
  if (!filepath || !key) return false;

  if (gHist) { delete gHist; gHist = nullptr; }

  TFile f(filepath, "READ");
  if (f.IsZombie()) {
    std::cerr << "[ACCGridPdf_Load] cannot open: " << filepath << "\n";
    return false;
  }

  TObject* obj = f.Get(key);

  if (!obj) {
    std::cerr << "[ACCGridPdf_Load] key not found: " << key << "\n";
    f.Close();
    return false;
  }

  THnD* h = dynamic_cast<THnD*>(obj);

  if (!h) {
    std::cerr << "[ACCGridPdf_Load] object is not THnD: " << key << "\n";
    f.Close();
    return false;
  }

  // ---- メタ情報の整合性チェック（警告のみ）----
  {
    const std::string kNtheta = std::string(key) + "_N_theta";
    TObject* par = f.Get(kNtheta.c_str());
    TParameter<int>* pN = dynamic_cast<TParameter<int>*>(par);
    if (pN) {
      const int fileN = pN->GetVal();
      const int nowN  = Math_GetNTheta(detres);
      if (fileN != nowN) {
        std::cerr << "[ACCGridPdf_Load] WARNING: N_theta mismatch (file="
                  << fileN << ", current detres=" << nowN
                  << "). ACC pdf may become sparse/zero.\n";
      }
    }
  }

  if (h->GetNdimensions() != 4) {
    std::cerr << "[ACCGridPdf_Load] THnD dimension is not 4: "
              << "N=" << h->GetNdimensions() << "\n";
    f.Close();
    return false;
  }

  gHist = dynamic_cast<THnD*>(h->Clone());
  f.Close();

  if (!gHist) {
    std::cerr << "[ACCGridPdf_Load] clone failed\n";
    if (gHist) { delete gHist; gHist = nullptr; }
    return false;
  }

  return true;
}

bool ACCGridPdf_IsLoaded() {
  return (gHist != nullptr);
}

double ACCGridPdf(double Ee, double Eg, double t,
                  double phi_detector_e, double phi_detector_g) {
  if (!gHist) return 0.0;

  if (!Math_IsFinite(Ee) || !Math_IsFinite(Eg) || !Math_IsFinite(t) ||
      !Math_IsFinite(phi_detector_e) || !Math_IsFinite(phi_detector_g)) {
    return 0.0;
  }

  const int N_theta = Math_GetNTheta(detres);
  const double phi_e_disc = Angle_DiscretizePhi(phi_detector_e, N_theta);
  const double phi_g_disc = Angle_DiscretizePhi(phi_detector_g, N_theta);

  // theta_eg = |phi_e - phi_g|
  const double theta_eg = Angle_ThetaFromPhiClipped(phi_e_disc, phi_g_disc);

  // 解析窓（Ee,Eg,t,theta）でカット
  if (!AnalysisWindow_In3D(analysis_window, Ee, Eg, theta_eg)) return 0.0;
  if (!AnalysisWindow_InTime(analysis_window, t)) return 0.0;

  // phi 軸は離散点のため FindBin で対応ビンに落とす
  const TAxis* axPe = gHist->GetAxis(2);
  const TAxis* axPg = gHist->GetAxis(3);
  const int bin_pe = axPe->FindBin(phi_e_disc);
  const int bin_pg = axPg->FindBin(phi_g_disc);

  // 4D格子から p4(Ee,Eg,phi_e,phi_g) を評価（Ee,Eg は補間）
  const double p4 = Hist_InterpEeEg4(*gHist, Ee, Eg, bin_pe, bin_pg);
  if (!(p4 > 0.0) || !Math_IsFinite(p4)) return 0.0;

  // 時間因子（解析窓内一様）
  const double dt = analysis_window.t_max - analysis_window.t_min;
  if (!(dt > 0.0) || !Math_IsFinite(dt)) return 0.0;
  const double pt = 1.0 / dt;

  const double pdf = p4 * pt;
  if (!(pdf > 0.0) || !Math_IsFinite(pdf)) return 0.0;
  return pdf;
}
