// src/MichelPolFit.cc
#include "p2meg/MichelPolFit.h"

#include <cmath>
#include <vector>
#include <string>
#include <limits>

#include "p2meg/MichelEData.h"
#include "p2meg/MichelPolConfig.h"

// MichelPolTemplate.cc で定義（ヘッダは作らない方針）
std::vector<double> BuildMichelPolKi(const MichelPolConfig& cfg);

static inline bool IsFinite(double x) {
    return std::isfinite(x);
}

static inline int BinIndex(double Ee, double Emin, double Emax, int nbins) {
    if (!(Emax > Emin) || !(nbins > 0)) return -1;
    if (!(Ee >= Emin) || !(Ee < Emax)) return -1;
    const double binw = (Emax - Emin) / static_cast<double>(nbins);
    const int idx = static_cast<int>((Ee - Emin) / binw);
    if (idx < 0 || idx >= nbins) return -1;
    return idx;
}

static inline double BinCenter(int i, double Emin, double Emax, int nbins) {
    const double binw = (Emax - Emin) / static_cast<double>(nbins);
    return Emin + (i + 0.5) * binw;
}

static void FillHistCounts(const std::vector<MichelEEvent>& evs,
                           const MichelPolConfig& cfg,
                           std::vector<double>& counts,
                           long long& n_inrange,
                           long long& n_outrange) {
    n_inrange = 0;
    n_outrange = 0;
    for (const auto& ev : evs) {
        const int ib = BinIndex(ev.Ee, cfg.Ee_min, cfg.Ee_max, cfg.nbins_Ee);
        if (ib >= 0) {
            counts[ib] += 1.0;
            n_inrange++;
        } else {
            n_outrange++;
        }
    }
}

static inline bool InFitRange(double Ecenter, const MichelPolConfig& cfg) {
    return (Ecenter >= cfg.fit_Ee_min && Ecenter <= cfg.fit_Ee_max);
}

static bool ComputeCrossRatioAsymAndErr(double a, double b, double c, double d,
                                       double& A, double& sigma_A) {
    // a,b,c,d は Poisson カウント（>=0）
    // 方針: ゼロ（または min_counts 未満）が含まれたら false（このビンを除外）
    if (!(a >= 0.0 && b >= 0.0 && c >= 0.0 && d >= 0.0)) return false;
    if (!(a >= 1.0 && b >= 1.0 && c >= 1.0 && d >= 1.0)) return false;

    const double ab = a * b;
    const double cd = c * d;
    if (!(ab > 0.0 && cd > 0.0) || !IsFinite(ab) || !IsFinite(cd)) return false;

    const double r = std::sqrt(ab / cd);
    if (!(r > 0.0) || !IsFinite(r)) return false;

    A = (r - 1.0) / (r + 1.0);
    if (!IsFinite(A)) return false;

    // 誤差伝播（log を介す近似）
    // Var(ln r) ≈ (1/4)(1/a+1/b+1/c+1/d)
    // u = (1/2) ln r,  A = tanh(u)
    // σ_A^2 ≈ (1-A^2)^2 * Var(u) = (1-A^2)^2 * Var(ln r)/4
    const double var_ln_r = 0.25 * (1.0/a + 1.0/b + 1.0/c + 1.0/d);
    const double var_u = 0.25 * var_ln_r;

    const double one_minus_A2 = 1.0 - A*A;
    const double var_A = (one_minus_A2 * one_minus_A2) * var_u;

    if (!(var_A > 0.0) || !IsFinite(var_A)) return false;
    sigma_A = std::sqrt(var_A);
    return (sigma_A > 0.0 && IsFinite(sigma_A));
}

MichelPolFitResult FitMichelPolarizationFrom4Files(const std::string& path_A_plus,
                                                   const std::string& path_A_minus,
                                                   const std::string& path_B_plus,
                                                   const std::string& path_B_minus) {
    MichelPolFitResult res{};
    res.status = 1;
    res.P_mu_hat = 0.0;
    res.err_P_mu = 0.0;
    res.chi2 = 0.0;
    res.ndf = 0;
    res.n_bins_used = 0;
    res.n_read_total = 0;
    res.n_skipped_total = 0;

    const MichelPolConfig& cfg = michel_pol_config;

    if (!(cfg.Ee_max > cfg.Ee_min) || !(cfg.nbins_Ee > 0)) return res;
    if (!(cfg.fit_Ee_max > cfg.fit_Ee_min)) return res;

    // --- 読み込み ---
    long long r1=0,s1=0, r2=0,s2=0, r3=0,s3=0, r4=0,s4=0;
    const auto ev_Ap = ReadMichelEData(path_A_plus,  &r1, &s1);
    const auto ev_Am = ReadMichelEData(path_A_minus, &r2, &s2);
    const auto ev_Bp = ReadMichelEData(path_B_plus,  &r3, &s3);
    const auto ev_Bm = ReadMichelEData(path_B_minus, &r4, &s4);

    res.n_read_total = r1 + r2 + r3 + r4;
    res.n_skipped_total = s1 + s2 + s3 + s4;

    const int nb = cfg.nbins_Ee;

    // --- ヒスト（カウント配列） ---
    std::vector<double> nAp(nb, 0.0), nAm(nb, 0.0), nBp(nb, 0.0), nBm(nb, 0.0);

    long long in1=0,out1=0,in2=0,out2=0,in3=0,out3=0,in4=0,out4=0;
    FillHistCounts(ev_Ap, cfg, nAp, in1, out1);
    FillHistCounts(ev_Am, cfg, nAm, in2, out2);
    FillHistCounts(ev_Bp, cfg, nBp, in3, out3);
    FillHistCounts(ev_Bm, cfg, nBm, in4, out4);

    // --- 理論テンプレ K_i ---
    const std::vector<double> K = BuildMichelPolKi(cfg);
    if (static_cast<int>(K.size()) != nb) return res;

    // --- 各ビンの A_data, sigma, K を作り、重み付き最小二乗 ---
    double S1 = 0.0;
    double S2 = 0.0;

    struct BinObs {
        int i;
        double A;
        double s;
        double K;
    };
    std::vector<BinObs> used;

    for (int i = 0; i < nb; ++i) {
        const double Ecen = BinCenter(i, cfg.Ee_min, cfg.Ee_max, cfg.nbins_Ee);
        if (!InFitRange(Ecen, cfg)) continue;

        // cross-ratio の定義
        // a=nA+(i), b=nB+(i), c=nA-(i), d=nB-(i)
        const double a = nAp[i];
        const double b = nBp[i];
        const double c = nAm[i];
        const double d = nBm[i];

        // min_counts_each を適用（0対策）
        if (a < cfg.min_counts_each || b < cfg.min_counts_each ||
            c < cfg.min_counts_each || d < cfg.min_counts_each) {
            continue;
        }

        double A = 0.0;
        double sA = 0.0;
        if (!ComputeCrossRatioAsymAndErr(a, b, c, d, A, sA)) continue;

        const double Ki = K[i];
        if (!IsFinite(Ki)) continue;

        const double w = 1.0 / (sA * sA);
        if (!(w > 0.0) || !IsFinite(w)) continue;

        S1 += Ki * A * w;
        S2 += Ki * Ki * w;

        used.push_back(BinObs{i, A, sA, Ki});
    }

    res.n_bins_used = static_cast<int>(used.size());
    if (res.n_bins_used < 2) return res;
    if (!(S2 > 0.0) || !IsFinite(S2) || !IsFinite(S1)) return res;

    const double P_hat = S1 / S2;
    const double errP = 1.0 / std::sqrt(S2);
    if (!IsFinite(P_hat) || !IsFinite(errP) || !(errP > 0.0)) return res;

    // chi2
    double chi2 = 0.0;
    for (const auto& u : used) {
        const double pull = (u.A - P_hat * u.K) / u.s;
        chi2 += pull * pull;
    }

    res.P_mu_hat = P_hat;
    res.err_P_mu = errP;
    res.chi2 = chi2;
    res.ndf = res.n_bins_used - 1;
    res.status = 0;
    return res;
}
