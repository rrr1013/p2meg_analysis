// src/MichelPolTemplate.cc
#include "p2meg/MichelPolConfig.h"

#include <cmath>
#include <vector>
#include <limits>

#include "p2meg/MichelSpectrum.h"
#include "p2meg/DetectorResolution.h"

// ============================================================
// 内部用: Michel 偏極測定の理論テンプレート K_i を作る
//
// 定義:
//   A_i^{th}(P_mu) = P_mu * K_i
//   K_i = |cosθ| * (V_i / U_i)
//
// ここで U_i, V_i は「真の U(E), V(E) を e+ のエネルギー応答で畳み込み」して
// 観測エネルギービン i に落とし込んだもの。
//   μ_±(i;P_mu) ∝ U_i ± (P_mu |cosθ|) V_i
//
// 近似:
//  - 数値積分（中点則）で畳み込みを評価
//  - 応答は energy_response_shape_e を用い、auto-range で正規化して使う
//    （smear_energy_trandom3_e と同じ正規化レンジ）
// ============================================================

// 外部公開しない想定だが、別 .cc から呼ぶために非staticで定義する。
std::vector<double> BuildMichelPolKi(const MichelPolConfig& cfg);

static inline bool IsFinite(double x) {
    return std::isfinite(x);
}

static inline double ClampNonNegative(double x) {
    if (!IsFinite(x)) return 0.0;
    return (x > 0.0) ? x : 0.0;
}

static double IntegrateResponseShapeOverInterval(double E_true,
                                                 double E_low,
                                                 double E_high,
                                                 int n_steps) {
    if (!(E_high > E_low)) return 0.0;
    if (!(n_steps >= 1)) n_steps = 1;

    const double h = (E_high - E_low) / static_cast<double>(n_steps);
    double sum = 0.0;
    for (int k = 0; k < n_steps; ++k) {
        const double x = E_low + (k + 0.5) * h;
        const double sh = energy_response_shape_e(x, E_true);
        sum += ClampNonNegative(sh);
    }
    return sum * h;
}

static inline double MichelUTrue(double E_true) {
    // 偏極なし成分（costh=0, P_mu=0）
    const double u = Michel_d2Shape_dE_dCosTheta(E_true, 0.0, 0.0);
    return (IsFinite(u) ? u : 0.0);
}

static inline double MichelVTrue(double E_true) {
    // 偏極項（P_mu*cosθ に掛かる係数）を差分で取り出す
    // V(E) ≡ Michel(E, costh=+1, P_mu=1) - Michel(E, costh=+1, P_mu=0)
    const double s1 = Michel_d2Shape_dE_dCosTheta(E_true, +1.0, 1.0);
    const double s0 = Michel_d2Shape_dE_dCosTheta(E_true, +1.0, 0.0);
    const double v = (IsFinite(s1) ? s1 : 0.0) - (IsFinite(s0) ? s0 : 0.0);
    return v;
}

std::vector<double> BuildMichelPolKi(const MichelPolConfig& cfg) {
    std::vector<double> K;

    if (!(cfg.Ee_max > cfg.Ee_min)) return K;
    if (!(cfg.nbins_Ee > 0)) return K;
    if (!(cfg.cos_theta_abs >= 0.0 && cfg.cos_theta_abs <= 1.0)) return K;

    const int nb = cfg.nbins_Ee;
    K.assign(nb, 0.0);

    std::vector<double> U(nb, 0.0);
    std::vector<double> V(nb, 0.0);

    const double binw = (cfg.Ee_max - cfg.Ee_min) / static_cast<double>(nb);

    // E_true 積分の設定（数値都合。必要なら後で増やす）
    const double Etrue_min = cfg.Ee_min;
    const double Etrue_max = cfg.Ee_max;
    const int N_true = 2000;  // 粗いと K_i がガタつく。重ければ下げる。
    const double dE = (Etrue_max - Etrue_min) / static_cast<double>(N_true);

    // Eobs 積分の分割数（各ビン内）
    const int N_obs_in_bin = 6;

    for (int it = 0; it < N_true; ++it) {
        const double E_true = Etrue_min + (it + 0.5) * dE;
        if (!(E_true >= 0.0) || !IsFinite(E_true)) continue;

        const double u_true = MichelUTrue(E_true);
        const double v_true = MichelVTrue(E_true);

        // 応答の正規化（auto-range）
        double Eres_min = 0.0;
        double Eres_max = 0.0;
        if (!energy_response_autorange_e(E_true, Eres_min, Eres_max)) continue;
        const double A = energy_response_integral_e(Eres_min, Eres_max, E_true);
        if (!(A > 0.0) || !IsFinite(A)) continue;

        // 応答のビン積分 Ibin[j]
        std::vector<double> Ibin(nb, 0.0);
        for (int j = 0; j < nb; ++j) {
            const double Elow  = cfg.Ee_min + j * binw;
            const double Ehigh = Elow + binw;
            const double lo = (Elow  > Eres_min ? Elow  : Eres_min);
            const double hi = (Ehigh < Eres_max ? Ehigh : Eres_max);
            if (!(hi > lo)) {
                Ibin[j] = 0.0;
                continue;
            }
            const double I = IntegrateResponseShapeOverInterval(E_true, lo, hi, N_obs_in_bin);
            Ibin[j] = I;
        }

        // p_bin = Ibin / A を使って U_i, V_i を更新
        for (int j = 0; j < nb; ++j) {
            const double p = Ibin[j] / A;
            if (!(p > 0.0) || !IsFinite(p)) continue;

            U[j] += u_true * p * dE;
            V[j] += v_true * p * dE;
        }
    }

    // K_i = |cosθ| * V_i / U_i
    for (int j = 0; j < nb; ++j) {
        if (U[j] > 0.0 && IsFinite(U[j]) && IsFinite(V[j])) {
            K[j] = cfg.cos_theta_abs * (V[j] / U[j]);
        } else {
            K[j] = 0.0;
        }
    }

    return K;
}
