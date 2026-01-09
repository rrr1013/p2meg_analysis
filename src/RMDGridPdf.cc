// src/RMDGridPdf.cc
#include "p2meg/RMDGridPdf.h"

#include <cmath>
#include <iostream>
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

static THnD* gHist4 = nullptr;

static bool IsFinite(double x) {
  return std::isfinite(x);
}

static double Clamp(double x, double lo, double hi) {
  if (x < lo) return lo;
  if (x > hi) return hi;
  return x;
}

//============================================================
// 角度離散化（cos -> theta -> 最近傍格子）
//============================================================

// theta_i = i*pi/N_theta (i=0..N_theta)
// 入力 theta を最近傍の格子インデックスに丸める（端はクリップ）
// 不正なら -1
static int ThetaNearestIndex(double theta, int N_theta) {
  if (!(N_theta >= 1)) return -1;
  if (!IsFinite(theta)) return -1;
  if (theta < 0.0 || theta > pi) return -1;

  const double step = pi / static_cast<double>(N_theta);
  if (!(step > 0.0) || !IsFinite(step)) return -1;

  long long i_ll = std::llround(theta / step);

  // 数値丸めの端のはみ出しを安全にクリップ（物理カットではない）
  if (i_ll < 0LL) i_ll = 0LL;
  if (i_ll > static_cast<long long>(N_theta)) i_ll = static_cast<long long>(N_theta);

  return static_cast<int>(i_ll);
}

// cos -> acos -> 最近傍格子で ie を決める
static int CosToIndex(double cosv, int N_theta) {
  if (!(N_theta >= 1)) return -1;
  if (!IsFinite(cosv)) return -1;

  const double c = Clamp(cosv, -1.0, 1.0);
  const double theta = std::acos(c);
  if (!IsFinite(theta)) return -1;

  return ThetaNearestIndex(theta, N_theta);
}

// (ie, ig) から相対角 thetaEG を作る（cosΔφ=+1 固定）
static bool ThetaEG_FromIndices(int ie, int ig, int N_theta, double& thetaEG_out) {
  if (!(N_theta >= 1)) return false;
  if (ie < 0 || ie > N_theta) return false;
  if (ig < 0 || ig > N_theta) return false;

  const double thetaE = pi * static_cast<double>(ie) / static_cast<double>(N_theta);
  const double thetaG = pi * static_cast<double>(ig) / static_cast<double>(N_theta);

  // cos(thetaE - thetaG)（cosΔφ=+1）
  const double cEG = std::cos(thetaE - thetaG);
  const double c = Clamp(cEG, -1.0, 1.0);
  const double th = std::acos(c);

  if (!IsFinite(th)) return false;
  thetaEG_out = th;
  return true;
}

//============================================================
// 軸補間（Ee,Eg だけ補間、ie/ig は固定）
//============================================================

// 軸が等間隔ビンである前提で、座標 x を隣り合う2ビンに落とす。
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

// 4D THnD から、(Ee,Eg) のみ双一次補間（4点）し、(ie,ig) は固定
static double Interp2D_4_FixedI(const THnD& h, double Ee, double Eg, int ie, int ig) {
  const TAxis* ax0 = h.GetAxis(0); // Ee
  const TAxis* ax1 = h.GetAxis(1); // Eg
  const TAxis* ax2 = h.GetAxis(2); // ie
  const TAxis* ax3 = h.GetAxis(3); // ig

  if (!ax0 || !ax1 || !ax2 || !ax3) return 0.0;

  // ie, ig は整数インデックスなので、対応するビン番号を確定
  // MakeRMDGridPdf 側では軸範囲 [-0.5, N+0.5] かつ nbins=N+1 を想定
  // よって ie=0 -> bin1, ie=1 -> bin2, ... になる
  const int nb_ie = ax2->GetNbins();
  const int nb_ig = ax3->GetNbins();

  const int ie_bin = ie + 1;
  const int ig_bin = ig + 1;

  if (ie_bin < 1 || ie_bin > nb_ie) return 0.0;
  if (ig_bin < 1 || ig_bin > nb_ig) return 0.0;

  int i0e = 0, i1e = 0;
  int i0g = 0, i1g = 0;
  double fe = 0.0, fg = 0.0;

  if (AxisBracketUniform(*ax0, Ee, i0e, i1e, fe) != 0) return 0.0;
  if (AxisBracketUniform(*ax1, Eg, i0g, i1g, fg) != 0) return 0.0;

  // 4点
  // (i0e,i0g), (i1e,i0g), (i0e,i1g), (i1e,i1g)
  double sum = 0.0;

  // idx: [EeBin, EgBin, ieBin, igBin]
  std::vector<int> idx(4, 1);
  idx[2] = ie_bin;
  idx[3] = ig_bin;

  for (int be = 0; be < 2; ++be) {
    const int EeBin = (be == 0) ? i0e : i1e;
    const double we = (be == 0) ? (1.0 - fe) : fe;

    idx[0] = EeBin;

    for (int bg = 0; bg < 2; ++bg) {
      const int EgBin = (bg == 0) ? i0g : i1g;
      const double wg = (bg == 0) ? (1.0 - fg) : fg;

      idx[1] = EgBin;

      const Long64_t bin = h.GetBin(idx.data());
      const double v = h.GetBinContent(bin);
      if (v > 0.0 && IsFinite(v)) {
        sum += (we * wg) * v;
      }
    }
  }

  return sum;
}

//============================================================
// 時間因子（窓内正規化ガウシアン）
//============================================================

// pt(t) = N(t_mean, sigma_t) / A_t, ただし A_t = ∫_{tmin}^{tmax} N dt
static double PtWindowNormalized(double t) {
  const double sigma_t = detres.sigma_t;
  const double t_mean  = detres.t_mean;

  const double tmin = analysis_window.t_min;
  const double tmax = analysis_window.t_max;

  if (!(sigma_t > 0.0) || !IsFinite(sigma_t)) return 0.0;
  if (!IsFinite(t) || !IsFinite(t_mean)) return 0.0;

  const double z = (t - t_mean) / sigma_t;

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

  if (gHist4) {
    delete gHist4;
    gHist4 = nullptr;
  }

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

  // 4D格子のみを受け付ける
  if (h->GetNdimensions() != 4) {
    std::cerr << "[RMDGridPdf_Load] THnD dimension is not 4: ndims=" << h->GetNdimensions() << "\n";
    f.Close();
    return false;
  }

  // 軸サイズの整合（N_theta+1 を期待）
  const int N_theta = detres.N_theta;
  if (!(N_theta >= 1)) {
    std::cerr << "[RMDGridPdf_Load] detres.N_theta must be >= 1\n";
    f.Close();
    return false;
  }

  const int nb_ie = h->GetAxis(2)->GetNbins();
  const int nb_ig = h->GetAxis(3)->GetNbins();
  if (nb_ie != (N_theta + 1) || nb_ig != (N_theta + 1)) {
    std::cerr << "[RMDGridPdf_Load] axis bins mismatch: "
              << "ie=" << nb_ie << ", ig=" << nb_ig
              << " (expected " << (N_theta + 1) << ")\n";
    f.Close();
    return false;
  }

  gHist4 = dynamic_cast<THnD*>(h->Clone());
  f.Close();

  if (!gHist4) {
    std::cerr << "[RMDGridPdf_Load] clone failed\n";
    return false;
  }

  return true;
}

bool RMDGridPdf_IsLoaded() {
  return (gHist4 != nullptr);
}

double RMDGridPdf(double Ee, double Eg, double t, double /*theta*/,
                  double cos_detector_e, double cos_detector_g) {
  if (!gHist4) return 0.0;

  // 基本チェック
  if (!IsFinite(Ee) || !IsFinite(Eg) || !IsFinite(t)) return 0.0;
  if (!IsFinite(cos_detector_e) || !IsFinite(cos_detector_g)) return 0.0;

  // 解析窓（Ee,Eg,t）
  if (Ee < analysis_window.Ee_min || Ee > analysis_window.Ee_max) return 0.0;
  if (Eg < analysis_window.Eg_min || Eg > analysis_window.Eg_max) return 0.0;
  if (t  < analysis_window.t_min  || t  > analysis_window.t_max ) return 0.0;

  const int N_theta = detres.N_theta;
  if (!(N_theta >= 1)) return 0.0;

  // cos -> (ie,ig) に丸め（角度離散化）
  const int ie = CosToIndex(cos_detector_e, N_theta);
  const int ig = CosToIndex(cos_detector_g, N_theta);
  if (ie < 0 || ig < 0) return 0.0;

  // (ie,ig) から thetaEG を作り、解析窓の theta カットを適用
  double thetaEG = 0.0;
  if (!ThetaEG_FromIndices(ie, ig, N_theta, thetaEG)) return 0.0;

  if (thetaEG < analysis_window.theta_min || thetaEG > analysis_window.theta_max) return 0.0;

  // 4D格子から p2(Ee,Eg|ie,ig) を補間（Ee,Eg のみ）
  const double p2 = Interp2D_4_FixedI(*gHist4, Ee, Eg, ie, ig);
  if (!(p2 > 0.0) || !IsFinite(p2)) return 0.0;

  // 時間因子を解析的に掛ける（窓内正規化）
  const double pt = PtWindowNormalized(t);
  if (!(pt > 0.0) || !IsFinite(pt)) return 0.0;

  const double pdf = p2 * pt;
  if (!(pdf > 0.0) || !IsFinite(pdf)) return 0.0;
  return pdf;
}
