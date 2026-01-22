#ifndef P2MEG_DETECTOR_RESOLUTION_H
#define P2MEG_DETECTOR_RESOLUTION_H

#include <cmath>
#include <limits>
#include <vector>
#include "TRandom3.h"

#include "p2meg/Constants.h"

// ============================================================
// p2MEG 分解能モデル
//
// 単位:
//  - Ee, Eg: MeV
//  - t     : ns
//
// 角度 theta の扱い:
//  - 測定設定に合わせて、角度は離散化する
//      theta_i = i * pi / N_theta   (i = 0..N_theta)
//    すなわち、0 から pi までを N_theta 分割した (N_theta+1) 点のみを許す
//  - 崩壊後の e, γ の散乱は無視し、離散化後の角度スメアはかけない
//
// 角度 phi の扱い:
//  - phi_detector_e/g は DetectorResolutionConst の phi_*_min/max, N_phi_* で離散化する
//  - phi==phi_max は overflow に落ちず最後のビンに入るよう丸める
//  - 許可領域（マスク）は Detector_IsAllowedPhiPairIndex で定義する
//
// エネルギー分解能:
//  - energy_response_shape_e/g(E_res, E_true) は
//      「真値 E_true に対する再構成 E_res の分布 shape（未正規化）」を返す
//  - e+ と γ で別の shape を持てるように分離する
//  - 正規化やスメア生成は内部で自動化する
//  - ユーザが変更するのは shape 部分だけでよい
//
// t_mean は Δt の平均値（現状は 0 以外にもなり得るので分離）
// ============================================================

struct DetectorResolutionConst {
    double sigma_t;   // [ns]
    int    N_theta;   // 角度分割数（theta_i = i*pi/N_theta, i=0..N_theta）
    double t_mean;    // [ns]  Δt の平均値
    double P_mu;      // muon polarization (signed, [-1,1])
    double phi_e_min; // [rad] phi_detector_e の最小値
    double phi_e_max; // [rad] phi_detector_e の最大値
    int    N_phi_e;   // phi_detector_e の分割数（i=0..N_phi_e）
    double phi_g_min; // [rad] phi_detector_g の最小値
    double phi_g_max; // [rad] phi_detector_g の最大値
    int    N_phi_g;   // phi_detector_g の分割数（i=0..N_phi_g）
};

inline constexpr DetectorResolutionConst detres{
    0.1561,  // sigma_t  [ns]
    18,      // N_theta  （例：0..pi を 18 分割 → 19 点）
    -0.1479, // t_mean [ns]
    -0.8,    // P_mu

    pi * (10.0/180.0),     // phi_e_min [rad]
    pi * (190.0/180.0),      // phi_e_max [rad]
    18,      // N_phi_e

    pi * (10.0/180.0),     // phi_g_min [rad]
    pi * (190.0/180.0),      // phi_g_max [rad]
    18       // N_phi_g
};

// ------------------------------------------------------------
// phi の離散化・ビン境界の共通ヘルパ
// ------------------------------------------------------------
static inline bool Detector_IsPhiRangeValid(double phi_min, double phi_max, int N_phi) {
    if (!(N_phi >= 1)) return false;
    if (!std::isfinite(phi_min) || !std::isfinite(phi_max)) return false;
    if (!(phi_min < phi_max)) return false;
    return true;
}

static inline double Detector_PhiStep(double phi_min, double phi_max, int N_phi) {
    if (!Detector_IsPhiRangeValid(phi_min, phi_max, N_phi)) return 0.0;
    const double step = (phi_max - phi_min) / static_cast<double>(N_phi);
    return (step > 0.0 && std::isfinite(step)) ? step : 0.0;
}

// phi を [phi_min, phi_max] にクリップ（不正入力は phi_min に落とす）
static inline double Detector_PhiClamp(double phi, double phi_min, double phi_max) {
    if (!std::isfinite(phi)) return phi_min;
    if (phi < phi_min) return phi_min;
    if (phi > phi_max) return phi_max;
    return phi;
}

// phi の離散点（格子点）を返す
//  - phi_i = phi_min + i * step (i=0..N_phi)
static inline double Detector_PhiGridPoint(int i, double phi_min, double phi_max, int N_phi) {
    const double step = Detector_PhiStep(phi_min, phi_max, N_phi);
    if (!(step > 0.0)) return phi_min;
    long long idx = static_cast<long long>(i);
    if (idx < 0LL) idx = 0LL;
    if (idx > static_cast<long long>(N_phi)) idx = static_cast<long long>(N_phi);
    return phi_min + step * static_cast<double>(idx);
}

// phi を離散格子に丸め、インデックス i を返す
//  - phi==phi_max は overflow に落とさず i=N_phi へ入れる
static inline int Detector_PhiIndexFromValue(double phi,
                                             double phi_min, double phi_max,
                                             int N_phi) {
    const double step = Detector_PhiStep(phi_min, phi_max, N_phi);
    if (!(step > 0.0)) return -1;
    const double p = Detector_PhiClamp(phi, phi_min, phi_max);
    long long idx = std::llround((p - phi_min) / step);
    if (idx < 0LL) idx = 0LL;
    if (idx > static_cast<long long>(N_phi)) idx = static_cast<long long>(N_phi);
    return static_cast<int>(idx);
}

// phi を離散格子点に丸める（物理カットではない）
static inline double Detector_PhiSnapToGrid(double phi,
                                            double phi_min, double phi_max,
                                            int N_phi) {
    const int idx = Detector_PhiIndexFromValue(phi, phi_min, phi_max, N_phi);
    if (idx < 0) return phi_min;
    return Detector_PhiGridPoint(idx, phi_min, phi_max, N_phi);
}

// 「最近傍格子」になるように phi 軸の可変ビン境界を作る
//  - edges サイズ: N_phi+2
//  - ROOT のビン [low, high) で phi_max を in-range に入れるため上端を僅かに広げる
static inline std::vector<double> Detector_PhiEdgesFromGrid(double phi_min,
                                                            double phi_max,
                                                            int N_phi) {
    std::vector<double> edges;
    edges.resize(static_cast<size_t>(N_phi + 2));

    if (!Detector_IsPhiRangeValid(phi_min, phi_max, N_phi)) {
        for (auto& v : edges) v = phi_min;
        return edges;
    }

    const double step = Detector_PhiStep(phi_min, phi_max, N_phi);
    const double phi_max_plus =
        std::nextafter(phi_max, std::numeric_limits<double>::infinity());

    edges[0] = phi_min;
    for (int i = 1; i <= N_phi; ++i) {
        edges[i] = phi_min + (static_cast<double>(i) - 0.5) * step;
    }
    edges[N_phi + 1] = phi_max_plus;

    for (int k = 1; k < static_cast<int>(edges.size()); ++k) {
        if (edges[k] < edges[k - 1]) edges[k] = edges[k - 1];
        if (edges[k] < phi_min) edges[k] = phi_min;
        if (edges[k] > phi_max_plus) edges[k] = phi_max_plus;
    }
    return edges;
}

// ROOT の [low, high) ビンで phi_max を overflow に落とさない上端補正
// （物理カットではない）
static inline double Detector_PhiAxisMaxInclusive(double phi_max) {
    if (!std::isfinite(phi_max)) return phi_max;
    return std::nextafter(phi_max, std::numeric_limits<double>::infinity());
}

// 端点補正を入れた phi のビン幅（離散格子の「面積」計算用）
static inline double Detector_PhiBinWidth(int i, double phi_min, double phi_max, int N_phi) {
    const double step = Detector_PhiStep(phi_min, phi_max, N_phi);
    if (!(step > 0.0)) return 0.0;
    if (i <= 0 || i >= N_phi) return 0.5 * step;
    return step;
}

// ------------------------------------------------------------
// phi_e, phi_g の許可マスク（必要ならユーザがここを書き換える）
// ------------------------------------------------------------
static inline bool Detector_IsAllowedPhiPairIndex(int i_e, int i_g,
                                                  const DetectorResolutionConst& res) {
    // まず範囲外は拒否
    if (i_e < 0 || i_g < 0) return false;
    if (i_e > res.N_phi_e || i_g > res.N_phi_g) return false;

    // マスク条件: どちらかが phi_max に対応する端点インデックス
    return (i_e == res.N_phi_e) || (i_g == res.N_phi_g);
}

static inline bool Detector_IsAllowedPhiPairValue(double phi_e, double phi_g,
                                                  const DetectorResolutionConst& res,
                                                  int& idx_e_out, int& idx_g_out) {
    idx_e_out = Detector_PhiIndexFromValue(phi_e, res.phi_e_min, res.phi_e_max, res.N_phi_e);
    idx_g_out = Detector_PhiIndexFromValue(phi_g, res.phi_g_min, res.phi_g_max, res.N_phi_g);
    if (idx_e_out < 0 || idx_g_out < 0) return false;
    return Detector_IsAllowedPhiPairIndex(idx_e_out, idx_g_out, res);
}

// ------------------------------------------------------------
// phi 設定から theta_eg の許可範囲を見積もる
// ------------------------------------------------------------
static inline bool Detector_ThetaRangeFromAllowedPhi(const DetectorResolutionConst& res,
                                                     double& theta_min_out,
                                                     double& theta_max_out) {
    // phi の離散格子と許可マスクから theta_eg=|phi_e-phi_g| の範囲を求める
    // 単位: rad
    // 不正入力や許可ペアなしは false（物理カットではない）
    theta_min_out = 0.0;
    theta_max_out = 0.0;

    const int N_phi_e = res.N_phi_e;
    const int N_phi_g = res.N_phi_g;
    if (!Detector_IsPhiRangeValid(res.phi_e_min, res.phi_e_max, N_phi_e)) return false;
    if (!Detector_IsPhiRangeValid(res.phi_g_min, res.phi_g_max, N_phi_g)) return false;

    double tmin = std::numeric_limits<double>::infinity();
    double tmax = -std::numeric_limits<double>::infinity();
    bool any = false;

    for (int ie = 0; ie <= N_phi_e; ++ie) {
        for (int ig = 0; ig <= N_phi_g; ++ig) {
            if (!Detector_IsAllowedPhiPairIndex(ie, ig, res)) continue;
            const double phi_e = Detector_PhiGridPoint(ie, res.phi_e_min, res.phi_e_max, N_phi_e);
            const double phi_g = Detector_PhiGridPoint(ig, res.phi_g_min, res.phi_g_max, N_phi_g);
            const double theta = std::fabs(phi_e - phi_g);
            if (!std::isfinite(theta)) continue;
            if (theta < tmin) tmin = theta;
            if (theta > tmax) tmax = theta;
            any = true;
        }
    }
    if (!any) return false;
    if (!std::isfinite(tmin) || !std::isfinite(tmax) || !(tmin <= tmax)) return false;

    theta_min_out = tmin;
    theta_max_out = tmax;
    return true;
}


// ------------------------------------------------------------
// エネルギー応答 shape（ユーザが自由に変更する部分）
//  - e+ と γ で別々に定義する
// ------------------------------------------------------------
inline double energy_response_shape_e(double E_res, double E_true) {

    // 例: E_true を中心とする幅 0.1*E_true のガウシアン（未正規化）
    //  - 単位: E_res, E_true ともに MeV
    //  - 不正入力や非物理は 0 を返す
    if (!(E_true > 0.0)) return 0.0;
    if (!std::isfinite(E_res) || !std::isfinite(E_true)) return 0.0;

    const double sigma = 0.1 * E_true;
    if (!(sigma > 0.0) || !std::isfinite(sigma)) return 0.0;

    const double z = (E_res - E_true) / sigma;
    const double p = std::exp(-0.5 * z * z);
    return std::isfinite(p) ? p : 0.0;
}

inline double energy_response_shape_g(double E_res, double E_true) {

    // 例: E_true を中心とする幅 0.1*E_true のガウシアン（未正規化）
    //  - 単位: E_res, E_true ともに MeV
    //  - 不正入力や非物理は 0 を返す
    if (!(E_true > 0.0)) return 0.0;
    if (!std::isfinite(E_res) || !std::isfinite(E_true)) return 0.0;

    const double sigma = 0.1 * E_true;
    if (!(sigma > 0.0) || !std::isfinite(sigma)) return 0.0;

    const double z = (E_res - E_true) / sigma;
    const double p = std::exp(-0.5 * z * z);
    return std::isfinite(p) ? p : 0.0;
}

// ------------------------------------------------------------
// 内部: energy_response_shape の共通処理（e/g で共用）
// ------------------------------------------------------------
static inline double energy_response_integral_impl(double E_min, double E_max, double E_true,
                                                    double (*shape)(double, double)) {
    if (shape == nullptr) return 0.0;
    if (!(E_true > 0.0)) return 0.0;
    if (!(E_min < E_max)) return 0.0;
    if (!std::isfinite(E_min) || !std::isfinite(E_max) || !std::isfinite(E_true)) return 0.0;

    const int n = 400;
    const double h = (E_max - E_min) / static_cast<double>(n);
    if (!(h > 0.0) || !std::isfinite(h)) return 0.0;

    double sum = 0.0;
    for (int i = 0; i <= n; ++i) {
        const double x = E_min + h * static_cast<double>(i);
        const double w = (i == 0 || i == n) ? 0.5 : 1.0;
        const double fx = shape(x, E_true);
        if (fx > 0.0 && std::isfinite(fx)) sum += w * fx;
    }
    const double A = h * sum;
    return (A > 0.0 && std::isfinite(A)) ? A : 0.0;
}

// ------------------------------------------------------------
// energy_response_shape の正規化レンジを自動決定
//  - 応答が単峰で E_true 周りに減衰する前提
// ------------------------------------------------------------
static inline bool energy_response_autorange_impl(double E_true, double& E_min, double& E_max,
                                                   double (*shape)(double, double)) {
    if (shape == nullptr) return false;
    if (!(E_true > 0.0)) return false;
    const double p0 = shape(E_true, E_true);
    if (!(p0 > 0.0) || !std::isfinite(p0)) return false;

    const double tail = 1e-6 * p0;
    double range = 0.2 * E_true;
    if (!(range > 0.0)) range = 1.0;

    for (int i = 0; i < 40; ++i) {
        const double lo = E_true - range;
        const double hi = E_true + range;
        const double plo = shape(lo, E_true);
        const double phi = shape(hi, E_true);
        if ((plo < tail || !std::isfinite(plo)) && (phi < tail || !std::isfinite(phi))) {
            E_min = lo;
            E_max = hi;
            return true;
        }
        range *= 2.0;
    }

    E_min = E_true - range;
    E_max = E_true + range;
    return true;
}

static inline double energy_response_pdf_impl(double E_res, double E_true,
                                              double (*shape)(double, double)) {
    if (shape == nullptr) return 0.0;
    double E_min = 0.0;
    double E_max = 0.0;
    if (!energy_response_autorange_impl(E_true, E_min, E_max, shape)) return 0.0;
    if (E_res < E_min || E_res > E_max) return 0.0;
    const double A = energy_response_integral_impl(E_min, E_max, E_true, shape);
    if (!(A > 0.0) || !std::isfinite(A)) return 0.0;
    const double p = shape(E_res, E_true) / A;
    return std::isfinite(p) ? p : 0.0;
}

static inline double energy_response_pdf_window_impl(double E_res, double E_true,
                                                     double E_min, double E_max,
                                                     double (*shape)(double, double)) {
    if (shape == nullptr) return 0.0;
    if (E_res < E_min || E_res > E_max) return 0.0;
    const double A = energy_response_integral_impl(E_min, E_max, E_true, shape);
    if (!(A > 0.0) || !std::isfinite(A)) return 0.0;
    const double p = shape(E_res, E_true) / A;
    return std::isfinite(p) ? p : 0.0;
}

// ------------------------------------------------------------
// energy_response の 0.1 倍点の距離を数値的に求める
//  - 応答が単峰で E_true 周りに対称に減衰する前提
//  - 数値ガードは「物理カットではない」ことに注意
// ------------------------------------------------------------
static inline double energy_response_offset_high_impl(double E_true,
                                                      double (*shape)(double, double)) {
    if (shape == nullptr) return 0.0;
    if (!(E_true > 0.0)) return 0.0;
    const double p0 = shape(E_true, E_true);
    if (!(p0 > 0.0) || !std::isfinite(p0)) return 0.0;
    const double target = 0.1 * p0;

    double lo = 0.0;
    double hi = 0.2 * E_true;
    if (!(hi > 0.0)) hi = 1.0;

    int expand = 0;
    while (expand < 50 && shape(E_true + hi, E_true) > target) {
        hi *= 2.0;
        ++expand;
    }

    for (int i = 0; i < 80; ++i) {
        const double mid = 0.5 * (lo + hi);
        const double p = shape(E_true + mid, E_true);
        if (p > target) lo = mid;
        else hi = mid;
    }
    return hi;
}

static inline double energy_response_offset_low_impl(double E_true,
                                                     double (*shape)(double, double)) {
    if (shape == nullptr) return 0.0;
    if (!(E_true > 0.0)) return 0.0;
    const double p0 = shape(E_true, E_true);
    if (!(p0 > 0.0) || !std::isfinite(p0)) return 0.0;
    const double target = 0.1 * p0;

    double lo = 0.0;
    double hi = 0.2 * E_true;
    if (!(hi > 0.0)) hi = 1.0;

    int expand = 0;
    while (expand < 50 && shape(E_true - hi, E_true) > target) {
        hi *= 2.0;
        ++expand;
    }

    for (int i = 0; i < 80; ++i) {
        const double mid = 0.5 * (lo + hi);
        const double p = shape(E_true - mid, E_true);
        if (p > target) lo = mid;
        else hi = mid;
    }
    return hi;
}

// ------------------------------------------------------------
// energy_response_shape を正規化した分布から乱数サンプル（スメア）
//  - 正規化は自動化され、ユーザは shape を変更するだけでよい
// ------------------------------------------------------------
static inline double energy_response_max_on_range_impl(double E_true, double E_min, double E_max,
                                                       double (*shape)(double, double)) {
    if (shape == nullptr) return 0.0;
    if (!(E_min < E_max)) return 0.0;
    const int n = 200;
    double pmax = 0.0;
    for (int i = 0; i <= n; ++i) {
        const double x = E_min + (E_max - E_min) * (static_cast<double>(i) / static_cast<double>(n));
        const double p = shape(x, E_true);
        if (p > pmax && std::isfinite(p)) pmax = p;
    }
    return pmax;
}

static inline double smear_energy_trandom3_impl(TRandom3& rng, double E_true,
                                                double (*shape)(double, double)) {
    if (shape == nullptr) return E_true;
    double E_min = 0.0;
    double E_max = 0.0;
    if (!energy_response_autorange_impl(E_true, E_min, E_max, shape)) return E_true;

    const double pmax = energy_response_max_on_range_impl(E_true, E_min, E_max, shape);
    if (!(pmax > 0.0) || !std::isfinite(pmax)) return E_true;

    for (int it = 0; it < 100000; ++it) {
        const double x = rng.Uniform(E_min, E_max);
        const double u = rng.Uniform(0.0, pmax);
        if (u < shape(x, E_true)) return x;
    }
    return E_true;
}

// ------------------------------------------------------------
// e+ 用の energy_response ラッパー
// ------------------------------------------------------------
inline double energy_response_integral_e(double E_min, double E_max, double E_true) {
    return energy_response_integral_impl(E_min, E_max, E_true, energy_response_shape_e);
}

inline bool energy_response_autorange_e(double E_true, double& E_min, double& E_max) {
    return energy_response_autorange_impl(E_true, E_min, E_max, energy_response_shape_e);
}

inline double energy_response_pdf_e(double E_res, double E_true) {
    return energy_response_pdf_impl(E_res, E_true, energy_response_shape_e);
}

inline double energy_response_pdf_window_e(double E_res, double E_true,
                                           double E_min, double E_max) {
    return energy_response_pdf_window_impl(E_res, E_true, E_min, E_max, energy_response_shape_e);
}

inline double energy_response_offset_high_e(double E_true) {
    return energy_response_offset_high_impl(E_true, energy_response_shape_e);
}

inline double energy_response_offset_low_e(double E_true) {
    return energy_response_offset_low_impl(E_true, energy_response_shape_e);
}

inline double energy_response_max_on_range_e(double E_true, double E_min, double E_max) {
    return energy_response_max_on_range_impl(E_true, E_min, E_max, energy_response_shape_e);
}

inline double smear_energy_trandom3_e(TRandom3& rng, double E_true) {
    return smear_energy_trandom3_impl(rng, E_true, energy_response_shape_e);
}

// ------------------------------------------------------------
// γ 用の energy_response ラッパー
// ------------------------------------------------------------
inline double energy_response_integral_g(double E_min, double E_max, double E_true) {
    return energy_response_integral_impl(E_min, E_max, E_true, energy_response_shape_g);
}

inline bool energy_response_autorange_g(double E_true, double& E_min, double& E_max) {
    return energy_response_autorange_impl(E_true, E_min, E_max, energy_response_shape_g);
}

inline double energy_response_pdf_g(double E_res, double E_true) {
    return energy_response_pdf_impl(E_res, E_true, energy_response_shape_g);
}

inline double energy_response_pdf_window_g(double E_res, double E_true,
                                           double E_min, double E_max) {
    return energy_response_pdf_window_impl(E_res, E_true, E_min, E_max, energy_response_shape_g);
}

inline double energy_response_offset_high_g(double E_true) {
    return energy_response_offset_high_impl(E_true, energy_response_shape_g);
}

inline double energy_response_offset_low_g(double E_true) {
    return energy_response_offset_low_impl(E_true, energy_response_shape_g);
}

inline double energy_response_max_on_range_g(double E_true, double E_min, double E_max) {
    return energy_response_max_on_range_impl(E_true, E_min, E_max, energy_response_shape_g);
}

inline double smear_energy_trandom3_g(TRandom3& rng, double E_true) {
    return smear_energy_trandom3_impl(rng, E_true, energy_response_shape_g);
}

#endif // P2MEG_DETECTOR_RESOLUTION_H
