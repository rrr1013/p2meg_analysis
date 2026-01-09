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

// +枝（cosΔφ=+1）と -枝（cosΔφ=-1）の 4D 格子
static THnD* gHistP = nullptr;
static THnD* gHistM = nullptr;

// ---- 正規化モード（MakeRMDGridPdf.cc 側と合わせること） ----
// false: plus-only（+枝のみ正規化で生成した格子を使う。評価も常に +枝。）
// true : both（(+枝 + -枝) で正規化して生成した格子を使う。評価は生thetaで枝をハード判定。）
static constexpr bool kBothBranchesNormalized = true;

static bool IsFinite(double x) { return std::isfinite(x); }

static double Clamp(double x, double lo, double hi) {
  if (x < lo) return lo;
  if (x > hi) return hi;
  return x;
}

static bool IsInsideWindow4D(double Ee, double Eg, double t, double theta) {
  if (Ee < analysis_window.Ee_min || Ee > analysis_window.Ee_max) return false;
  if (Eg < analysis_window.Eg_min || Eg > analysis_window.Eg_max) return false;
  if (t  < analysis_window.t_min  || t  > analysis_window.t_max ) return false;
  if (theta < analysis_window.theta_min || theta > analysis_window.theta_max) return false;
  return true;
}

static int GetNTheta() {
  const double x = static_cast<double>(detres.N_theta);
  int N = static_cast<int>(std::lround(x));
  if (N < 1) N = 1;
  return N;
}

// θ を 0..π の (π/N) 格子へ最近傍丸め
static double QuantizeTheta(double theta_raw, int N_theta) {
  if (!IsFinite(theta_raw)) return 0.0;
  const double th = Clamp(theta_raw, 0.0, pi);
  const double step = pi / static_cast<double>(N_theta);
  if (!(step > 0.0) || !IsFinite(step)) return 0.0;

  const double kf = th / step;
  int k = static_cast<int>(std::lround(kf));
  if (k < 0) k = 0;
  if (k > N_theta) k = N_theta;
  return static_cast<double>(k) * step;
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

// Ee,Eg を 2D（4点）補間して、cos軸は固定ビン（最近傍）で評価する。
static double InterpEeEg_4(const THnD& h, double Ee, double Eg, int bin_cos_e, int bin_cos_g) {
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
    idx[2] = bin_cos_e;
    idx[3] = bin_cos_g;
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

  if (gHistP) { delete gHistP; gHistP = nullptr; }
  if (gHistM) { delete gHistM; gHistM = nullptr; }

  TFile f(filepath, "READ");
  if (f.IsZombie()) {
    std::cerr << "[RMDGridPdf_Load] cannot open: " << filepath << "\n";
    return false;
  }

  const std::string key_p = std::string(key) + "_p";
  const std::string key_m = std::string(key) + "_m";

  TObject* objP = f.Get(key_p.c_str());
  TObject* objM = f.Get(key_m.c_str());

  if (!objP) {
    std::cerr << "[RMDGridPdf_Load] key not found: " << key_p << "\n";
    f.Close();
    return false;
  }
  if (!objM) {
    std::cerr << "[RMDGridPdf_Load] key not found: " << key_m << "\n";
    f.Close();
    return false;
  }

  THnD* hP = dynamic_cast<THnD*>(objP);
  THnD* hM = dynamic_cast<THnD*>(objM);

  if (!hP || !hM) {
    std::cerr << "[RMDGridPdf_Load] object is not THnD: " << key_p << " or " << key_m << "\n";
    f.Close();
    return false;
  }

  if (hP->GetNdimensions() != 4 || hM->GetNdimensions() != 4) {
    std::cerr << "[RMDGridPdf_Load] THnD dimension is not 4: "
              << "P=" << hP->GetNdimensions() << " M=" << hM->GetNdimensions() << "\n";
    f.Close();
    return false;
  }

  gHistP = dynamic_cast<THnD*>(hP->Clone());
  gHistM = dynamic_cast<THnD*>(hM->Clone());
  f.Close();

  if (!gHistP || !gHistM) {
    std::cerr << "[RMDGridPdf_Load] clone failed\n";
    if (gHistP) { delete gHistP; gHistP = nullptr; }
    if (gHistM) { delete gHistM; gHistM = nullptr; }
    return false;
  }

  return true;
}

bool RMDGridPdf_IsLoaded() {
  return (gHistP != nullptr && gHistM != nullptr);
}

double RMDGridPdf(double Ee, double Eg, double t, double theta,
                  double cos_detector_e, double cos_detector_g) {
  if (!gHistP || !gHistM) return 0.0;

  if (!IsFinite(Ee) || !IsFinite(Eg) || !IsFinite(t) || !IsFinite(theta) ||
      !IsFinite(cos_detector_e) || !IsFinite(cos_detector_g)) {
    return 0.0;
  }

  // theta は 0..pi に丸めてから最近傍格子点へ（run_nll_fit 側で丸めていても二重で害は小さい）
  const int N_theta = GetNTheta();
  const double theta_q = QuantizeTheta(theta, N_theta);

  // 解析窓（4D: Ee,Eg,t,theta）でカット
  if (!IsInsideWindow4D(Ee, Eg, t, theta_q)) return 0.0;

  // cos は物理範囲にクリップ
  const double ce_in = Clamp(cos_detector_e, -1.0, 1.0);
  const double cg_in = Clamp(cos_detector_g, -1.0, 1.0);

  // cos 軸は「最近傍θ格子」になる可変ビンのはずなので FindBin で最近傍に落ちる
  const TAxis* axCe = gHistP->GetAxis(2);
  const TAxis* axCg = gHistP->GetAxis(3);
  const int bin_ce = axCe->FindBin(ce_in);
  const int bin_cg = axCg->FindBin(cg_in);

  // ビン中心の cos を採用（生成側と整合）
  const double cosE = axCe->GetBinCenter(bin_ce);
  const double cosG = axCg->GetBinCenter(bin_cg);

  // θe, θg（0..π）
  const double th_e = std::acos(Clamp(cosE, -1.0, 1.0));
  const double th_g = std::acos(Clamp(cosG, -1.0, 1.0));

  const double sinE = std::sin(th_e);
  const double sinG = std::sin(th_g);

  // 2枝の θeγ 候補（いずれも 0..π）
  const double cosEG_p = Clamp(cosE * cosG + sinE * sinG, -1.0, 1.0); // cosΔφ=+1
  const double cosEG_m = Clamp(cosE * cosG - sinE * sinG, -1.0, 1.0); // cosΔφ=-1
  const double thEG_p  = std::acos(cosEG_p);
  const double thEG_m  = std::acos(cosEG_m);

  // 候補も同じ θ 格子へ最近傍丸めしてから比較（ハード判定）
  const double thEG_p_q = QuantizeTheta(thEG_p, N_theta);
  const double thEG_m_q = QuantizeTheta(thEG_m, N_theta);

  const double dp = std::fabs(theta_q - thEG_p_q);
  const double dm = std::fabs(theta_q - thEG_m_q);

  // 使用する枝（plus-only のときは常に +枝）
  const THnD* hUse = gHistP;
  if (kBothBranchesNormalized) {
    hUse = (dp <= dm) ? gHistP : gHistM;
  } else {
    hUse = gHistP;
  }

  // 4D格子から p4(Ee,Eg,cosE,cosG) を評価（Ee,Eg は補間、cos は固定ビン）
  const double p4 = InterpEeEg_4(*hUse, Ee, Eg, bin_ce, bin_cg);
  if (!(p4 > 0.0) || !IsFinite(p4)) return 0.0;

  // 時間因子を解析的に掛ける（窓内正規化）
  const double pt = PtWindowNormalized(t);
  if (!(pt > 0.0) || !IsFinite(pt)) return 0.0;

  const double pdf = p4 * pt;
  if (!(pdf > 0.0) || !IsFinite(pdf)) return 0.0;
  return pdf;
}
