#include "p2meg/RMDGridPdf.h"

#include <cmath>
#include <iostream>
#include <vector>

#include "TFile.h"
#include "THn.h"
#include "TAxis.h"

#include "p2meg/AnalysisWindow.h"
#include "p2meg/DetectorResolution.h"

//============================================================
// 内部状態
//============================================================

static THnD* gHist3 = nullptr;

static bool IsFinite(double x) {
    return std::isfinite(x);
}

static bool IsInsideWindow4D(double Ee, double Eg, double t, double theta) {
    if (Ee < analysis_window_rmd.Ee_min || Ee > analysis_window_rmd.Ee_max) return false;
    if (Eg < analysis_window_rmd.Eg_min || Eg > analysis_window_rmd.Eg_max) return false;
    if (t  < analysis_window_rmd.t_min  || t  > analysis_window_rmd.t_max ) return false;
    if (theta < analysis_window_rmd.theta_min || theta > analysis_window_rmd.theta_max) return false;
    return true;
}

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

// 3D multilinear interpolation（8点）
static double Interp3D_8(const THnD& h, double Ee, double Eg, double theta) {
    const TAxis* ax0 = h.GetAxis(0);
    const TAxis* ax1 = h.GetAxis(1);
    const TAxis* ax2 = h.GetAxis(2);

    int i0[3], i1[3];
    double f[3];

    if (AxisBracketUniform(*ax0, Ee,    i0[0], i1[0], f[0]) != 0) return 0.0;
    if (AxisBracketUniform(*ax1, Eg,    i0[1], i1[1], f[1]) != 0) return 0.0;
    if (AxisBracketUniform(*ax2, theta, i0[2], i1[2], f[2]) != 0) return 0.0;

    double sum = 0.0;
    std::vector<int> idx(3, 1);

    for (int mask = 0; mask < 8; ++mask) {
        double w = 1.0;
        for (int d = 0; d < 3; ++d) {
            const bool upper = ((mask >> d) & 1) != 0;
            if (upper) {
                idx[d] = i1[d];
                w *= f[d];
            } else {
                idx[d] = i0[d];
                w *= (1.0 - f[d]);
            }
        }
        const Long64_t bin = h.GetBin(idx.data());
        const double v = h.GetBinContent(bin);
        if (v > 0.0 && IsFinite(v)) sum += w * v;
    }

    return sum;
}

// 窓内で正規化された時間ガウシアン（密度）
// pt(t) = N(t_mean, sigma_t) / A_t, ただし A_t = ∫_{tmin}^{tmax} N dt
static double PtWindowNormalized(double t) {
    const double sigma_t = detres_rmd.sigma_t;
    const double t_mean  = detres_rmd.t_mean;

    const double tmin = analysis_window_rmd.t_min;
    const double tmax = analysis_window_rmd.t_max;

    if (!(sigma_t > 0.0) || !IsFinite(sigma_t)) return 0.0;

    const double inv_s = 1.0 / sigma_t;
    const double z = (t - t_mean) * inv_s;

    const double pi = 3.14159265358979323846;
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

    if (gHist3) {
        delete gHist3;
        gHist3 = nullptr;
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

    // 3D格子のみを受け付ける
    if (h->GetNdimensions() != 3) {
        std::cerr << "[RMDGridPdf_Load] THnD dimension is not 3: ndims=" << h->GetNdimensions() << "\n";
        f.Close();
        return false;
    }

    gHist3 = dynamic_cast<THnD*>(h->Clone());
    f.Close();

    if (!gHist3) {
        std::cerr << "[RMDGridPdf_Load] clone failed\n";
        return false;
    }

    return true;
}

bool RMDGridPdf_IsLoaded() {
    return (gHist3 != nullptr);
}

double RMDGridPdf(double Ee, double Eg, double t, double theta) {
    if (!gHist3) return 0.0;

    // 解析窓（4D）でカット
    if (!IsInsideWindow4D(Ee, Eg, t, theta)) return 0.0;

    // 3D格子から p3(Ee,Eg,theta) を補間
    const double p3 = Interp3D_8(*gHist3, Ee, Eg, theta);
    if (!(p3 > 0.0) || !IsFinite(p3)) return 0.0;

    // 時間因子を解析的に掛ける（窓内正規化）
    const double pt = PtWindowNormalized(t);
    if (!(pt > 0.0) || !IsFinite(pt)) return 0.0;

    const double pdf = p3 * pt;
    if (!(pdf > 0.0) || !IsFinite(pdf)) return 0.0;
    return pdf;
}
