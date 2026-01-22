// src/RMDGridPdf.cc
#include "p2meg/RMDGridPdf.h"

#include <cmath>
#include <iostream>
#include <string>
#include <limits>

#include "TFile.h"
#include "TH2.h"
#include "THn.h"
#include "TAxis.h"
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

  // ---- メタ情報の整合性チェック（警告のみ）----
  {
    const std::string kNphiE = std::string(key) + "_N_phi_e";
    const std::string kNphiG = std::string(key) + "_N_phi_g";
    TParameter<int>* pNphiE = dynamic_cast<TParameter<int>*>(f.Get(kNphiE.c_str()));
    TParameter<int>* pNphiG = dynamic_cast<TParameter<int>*>(f.Get(kNphiG.c_str()));
    if (pNphiE) {
      const int fileN = pNphiE->GetVal();
      const int nowN  = Math_GetNPhiE(detres);
      if (fileN != nowN) {
        std::cerr << "[RMDGridPdf_Load] WARNING: N_phi_e mismatch (file="
                  << fileN << ", current detres=" << nowN
                  << ").\n";
      }
    }
    if (pNphiG) {
      const int fileN = pNphiG->GetVal();
      const int nowN  = Math_GetNPhiG(detres);
      if (fileN != nowN) {
        std::cerr << "[RMDGridPdf_Load] WARNING: N_phi_g mismatch (file="
                  << fileN << ", current detres=" << nowN
                  << ").\n";
      }
    }

    const std::string kPhiEMin = std::string(key) + "_phi_e_min";
    const std::string kPhiEMax = std::string(key) + "_phi_e_max";
    const std::string kPhiGMin = std::string(key) + "_phi_g_min";
    const std::string kPhiGMax = std::string(key) + "_phi_g_max";
    TParameter<double>* pPhiEMin = dynamic_cast<TParameter<double>*>(f.Get(kPhiEMin.c_str()));
    TParameter<double>* pPhiEMax = dynamic_cast<TParameter<double>*>(f.Get(kPhiEMax.c_str()));
    TParameter<double>* pPhiGMin = dynamic_cast<TParameter<double>*>(f.Get(kPhiGMin.c_str()));
    TParameter<double>* pPhiGMax = dynamic_cast<TParameter<double>*>(f.Get(kPhiGMax.c_str()));
    if (pPhiEMin && pPhiEMax) {
      if (pPhiEMin->GetVal() != detres.phi_e_min ||
          pPhiEMax->GetVal() != detres.phi_e_max) {
        std::cerr << "[RMDGridPdf_Load] WARNING: phi_e range mismatch.\n";
      }
    }
    if (pPhiGMin && pPhiGMax) {
      if (pPhiGMin->GetVal() != detres.phi_g_min ||
          pPhiGMax->GetVal() != detres.phi_g_max) {
        std::cerr << "[RMDGridPdf_Load] WARNING: phi_g range mismatch.\n";
      }
    }

    const std::string kMask = std::string(key) + "_phi_mask";
    TH2I* hmask = dynamic_cast<TH2I*>(f.Get(kMask.c_str()));
    if (hmask) {
      bool mismatch = false;
      const int nxe = hmask->GetXaxis()->GetNbins();
      const int nxg = hmask->GetYaxis()->GetNbins();
      const int nowNe = Math_GetNPhiE(detres);
      const int nowNg = Math_GetNPhiG(detres);
      if (nxe != nowNe + 1 || nxg != nowNg + 1) {
        mismatch = true;
      } else {
        for (int ie = 0; ie <= nowNe && !mismatch; ++ie) {
          for (int ig = 0; ig <= nowNg; ++ig) {
            const int fileAllowed = static_cast<int>(hmask->GetBinContent(ie + 1, ig + 1));
            const int nowAllowed = Detector_IsAllowedPhiPairIndex(ie, ig, detres) ? 1 : 0;
            if (fileAllowed != nowAllowed) {
              mismatch = true;
              break;
            }
          }
        }
      }
      if (mismatch) {
        std::cerr << "[RMDGridPdf_Load] WARNING: phi mask mismatch with current detres.\n";
      }
    }
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

  int idx_e = -1;
  int idx_g = -1;
  if (!Detector_IsAllowedPhiPairValue(phi_detector_e, phi_detector_g, detres, idx_e, idx_g)) {
    return 0.0;
  }

  const int N_phi_e = Math_GetNPhiE(detres);
  const int N_phi_g = Math_GetNPhiG(detres);
  const double phi_e_disc =
      Detector_PhiGridPoint(idx_e, detres.phi_e_min, detres.phi_e_max, N_phi_e);
  const double phi_g_disc =
      Detector_PhiGridPoint(idx_g, detres.phi_g_min, detres.phi_g_max, N_phi_g);

  // theta_eg = |phi_e - phi_g|
  const double theta_eg = std::fabs(phi_e_disc - phi_g_disc);

  // 解析窓（4D: Ee,Eg,t,theta）でカット
  if (!AnalysisWindow_In4D(analysis_window, Ee, Eg, t, theta_eg)) return 0.0;

  // phi 軸は可変ビンなので FindBin で最近傍に落ちる
  const TAxis* axPe = gHist->GetAxis(2);
  const TAxis* axPg = gHist->GetAxis(3);
  const int bin_pe = axPe->FindBin(phi_e_disc);
  const int bin_pg = axPg->FindBin(phi_g_disc);

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
