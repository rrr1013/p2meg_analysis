// src/ACCGridPdf.cc
#include "p2meg/ACCGridPdf.h"

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

// phi を [0,pi] にクリップする
static double ClipPhi0Pi(double phi) {
  if (!IsFinite(phi)) return 0.0;
  if (phi < 0.0) return 0.0;
  if (phi > pi) return pi;
  return phi;
}

// detres.N_theta を int として使う（型が double でもここで丸める）
static int GetNTheta() {
  const double x = static_cast<double>(detres.N_theta);
  int N = static_cast<int>(std::lround(x));
  if (N < 1) N = 1;
  return N;
}

// phi を [0,pi] にクリップして N_theta の離散点に丸める
static double DiscretizePhi(double phi, int N_theta) {
  const double p = ClipPhi0Pi(phi);
  const double step = pi / static_cast<double>(N_theta);
  if (!(step > 0.0) || !IsFinite(step)) return 0.0;

  long long i = std::llround(p / step);
  if (i < 0LL) i = 0LL;
  if (i > static_cast<long long>(N_theta)) i = static_cast<long long>(N_theta);
  return step * static_cast<double>(i);
}

// 離散化後の phi から theta_eg = |phi_e - phi_g| を作る
static double ThetaEgFromPhi(double phi_e_disc, double phi_g_disc) {
  if (!IsFinite(phi_e_disc) || !IsFinite(phi_g_disc)) return 0.0;
  return std::fabs(phi_e_disc - phi_g_disc);
}

// 解析窓の Ee/Eg/theta のみ判定
static bool InEnergyAngleWindow3D(double Ee, double Eg, double theta) {
  if (Ee < analysis_window.Ee_min || Ee > analysis_window.Ee_max) return false;
  if (Eg < analysis_window.Eg_min || Eg > analysis_window.Eg_max) return false;
  if (theta < analysis_window.theta_min || theta > analysis_window.theta_max) return false;
  return true;
}

static bool InTimeWindow(double t) {
  if (t < analysis_window.t_min || t > analysis_window.t_max) return false;
  return true;
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

  if (!IsFinite(Ee) || !IsFinite(Eg) || !IsFinite(t) ||
      !IsFinite(phi_detector_e) || !IsFinite(phi_detector_g)) {
    return 0.0;
  }

  const int N_theta = GetNTheta();
  const double phi_e_disc = DiscretizePhi(phi_detector_e, N_theta);
  const double phi_g_disc = DiscretizePhi(phi_detector_g, N_theta);

  // theta_eg = |phi_e - phi_g|
  const double theta_eg = ThetaEgFromPhi(phi_e_disc, phi_g_disc);

  // 解析窓（Ee,Eg,t,theta）でカット
  if (!InEnergyAngleWindow3D(Ee, Eg, theta_eg)) return 0.0;
  if (!InTimeWindow(t)) return 0.0;

  // phi 軸は離散点のため FindBin で対応ビンに落とす
  const TAxis* axPe = gHist->GetAxis(2);
  const TAxis* axPg = gHist->GetAxis(3);
  const int bin_pe = axPe->FindBin(phi_e_disc);
  const int bin_pg = axPg->FindBin(phi_g_disc);

  // 4D格子から p4(Ee,Eg,phi_e,phi_g) を評価（Ee,Eg は補間）
  const double p4 = InterpEeEg_4(*gHist, Ee, Eg, bin_pe, bin_pg);
  if (!(p4 > 0.0) || !IsFinite(p4)) return 0.0;

  // 時間因子（解析窓内一様）
  const double dt = analysis_window.t_max - analysis_window.t_min;
  if (!(dt > 0.0) || !IsFinite(dt)) return 0.0;
  const double pt = 1.0 / dt;

  const double pdf = p4 * pt;
  if (!(pdf > 0.0) || !IsFinite(pdf)) return 0.0;
  return pdf;
}
