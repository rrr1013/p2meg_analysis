#ifndef P2MEG_DETECTOR_RESOLUTION_H
#define P2MEG_DETECTOR_RESOLUTION_H

#include <cmath>
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
};

inline constexpr DetectorResolutionConst detres{
    0.1561,  // sigma_t  [ns]
    18,      // N_theta  （例：0..pi を 18 分割 → 19 点）
    -0.1479, // t_mean [ns]
    -0.8     // P_mu
};


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
