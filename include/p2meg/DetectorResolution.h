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
//  - energy_response_shape(E_res, E_true) は
//      「真値 E_true に対する再構成 E_res の分布 shape（未正規化）」を返す
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
// ------------------------------------------------------------
inline double energy_response_shape(double E_res, double E_true) {

    // 例: E_true を中心とする幅 0.1*E_true のガウシアン（未正規化）
    //  - 単位: E_res, E_true ともに MeV
    //  - 不正入力や非物理は 0 を返す
    if (!(E_true > 0.0)) return 0.0;
    if (!std::isfinite(E_res) || !std::isfinite(E_true)) return 0.0;

    //const double sigma = 0.1 * E_true;
    const double sigma = 1;
    if (!(sigma > 0.0) || !std::isfinite(sigma)) return 0.0;

    const double z = (E_res - E_true) / sigma;
    const double p = std::exp(-0.5 * z * z);
    return std::isfinite(p) ? p : 0.0;
}







// ------------------------------------------------------------
// energy_response_shape の数値積分（未正規化）
//  - 数値積分は分解能モデルの詳細に依存するため、最小限の台形則で扱う
//  - 物理カットではなく数値ガードとして使う
// ------------------------------------------------------------
inline double energy_response_integral(double E_min, double E_max, double E_true) {
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
        const double fx = energy_response_shape(x, E_true);
        if (fx > 0.0 && std::isfinite(fx)) sum += w * fx;
    }
    const double A = h * sum;
    return (A > 0.0 && std::isfinite(A)) ? A : 0.0;
}

// ------------------------------------------------------------
// energy_response_shape の正規化レンジを自動決定
//  - 応答が単峰で E_true 周りに減衰する前提
// ------------------------------------------------------------
inline bool energy_response_autorange(double E_true, double& E_min, double& E_max) {
    if (!(E_true > 0.0)) return false;
    const double p0 = energy_response_shape(E_true, E_true);
    if (!(p0 > 0.0) || !std::isfinite(p0)) return false;

    const double tail = 1e-6 * p0;
    double range = 0.2 * E_true;
    if (!(range > 0.0)) range = 1.0;

    for (int i = 0; i < 40; ++i) {
        const double lo = E_true - range;
        const double hi = E_true + range;
        const double plo = energy_response_shape(lo, E_true);
        const double phi = energy_response_shape(hi, E_true);
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

inline double energy_response_pdf(double E_res, double E_true) {
    double E_min = 0.0;
    double E_max = 0.0;
    if (!energy_response_autorange(E_true, E_min, E_max)) return 0.0;
    if (E_res < E_min || E_res > E_max) return 0.0;
    const double A = energy_response_integral(E_min, E_max, E_true);
    if (!(A > 0.0) || !std::isfinite(A)) return 0.0;
    const double p = energy_response_shape(E_res, E_true) / A;
    return std::isfinite(p) ? p : 0.0;
}

inline double energy_response_pdf_window(double E_res, double E_true,
                                         double E_min, double E_max) {
    if (E_res < E_min || E_res > E_max) return 0.0;
    const double A = energy_response_integral(E_min, E_max, E_true);
    if (!(A > 0.0) || !std::isfinite(A)) return 0.0;
    const double p = energy_response_shape(E_res, E_true) / A;
    return std::isfinite(p) ? p : 0.0;
}

// ------------------------------------------------------------
// energy_response の 0.1 倍点の距離を数値的に求める
//  - 応答が単峰で E_true 周りに対称に減衰する前提
//  - 数値ガードは「物理カットではない」ことに注意
// ------------------------------------------------------------
inline double energy_response_offset_high(double E_true) {
    if (!(E_true > 0.0)) return 0.0;
    const double p0 = energy_response_shape(E_true, E_true);
    if (!(p0 > 0.0) || !std::isfinite(p0)) return 0.0;
    const double target = 0.1 * p0;

    double lo = 0.0;
    double hi = 0.2 * E_true;
    if (!(hi > 0.0)) hi = 1.0;

    int expand = 0;
    while (expand < 50 && energy_response_shape(E_true + hi, E_true) > target) {
        hi *= 2.0;
        ++expand;
    }

    for (int i = 0; i < 80; ++i) {
        const double mid = 0.5 * (lo + hi);
        const double p = energy_response_shape(E_true + mid, E_true);
        if (p > target) lo = mid;
        else hi = mid;
    }
    return hi;
}

inline double energy_response_offset_low(double E_true) {
    if (!(E_true > 0.0)) return 0.0;
    const double p0 = energy_response_shape(E_true, E_true);
    if (!(p0 > 0.0) || !std::isfinite(p0)) return 0.0;
    const double target = 0.1 * p0;

    double lo = 0.0;
    double hi = 0.2 * E_true;
    if (!(hi > 0.0)) hi = 1.0;

    int expand = 0;
    while (expand < 50 && energy_response_shape(E_true - hi, E_true) > target) {
        hi *= 2.0;
        ++expand;
    }

    for (int i = 0; i < 80; ++i) {
        const double mid = 0.5 * (lo + hi);
        const double p = energy_response_shape(E_true - mid, E_true);
        if (p > target) lo = mid;
        else hi = mid;
    }
    return hi;
}

// ------------------------------------------------------------
// energy_response_shape を正規化した分布から乱数サンプル（スメア）
//  - 正規化は自動化され、ユーザは shape を変更するだけでよい
// ------------------------------------------------------------
inline double energy_response_max_on_range(double E_true, double E_min, double E_max) {
    if (!(E_min < E_max)) return 0.0;
    const int n = 200;
    double pmax = 0.0;
    for (int i = 0; i <= n; ++i) {
        const double x = E_min + (E_max - E_min) * (static_cast<double>(i) / static_cast<double>(n));
        const double p = energy_response_shape(x, E_true);
        if (p > pmax && std::isfinite(p)) pmax = p;
    }
    return pmax;
}

inline double smear_energy_trandom3(TRandom3& rng, double E_true) {
    double E_min = 0.0;
    double E_max = 0.0;
    if (!energy_response_autorange(E_true, E_min, E_max)) return E_true;

    const double pmax = energy_response_max_on_range(E_true, E_min, E_max);
    if (!(pmax > 0.0) || !std::isfinite(pmax)) return E_true;

    for (int it = 0; it < 100000; ++it) {
        const double x = rng.Uniform(E_min, E_max);
        const double u = rng.Uniform(0.0, pmax);
        if (u < energy_response_shape(x, E_true)) return x;
    }
    return E_true;
}

#endif // P2MEG_DETECTOR_RESOLUTION_H
