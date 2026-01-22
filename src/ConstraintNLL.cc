#include "p2meg/Likelihood.h"

#include <cmath>

// ============================================================
// 制約項（NLLに加算）
//
// 制約が不要な間は常に 0 を返す。
// 後から制約を入れたくなったら、この関数の中身を編集して足し算する。
// ============================================================

double ConstraintNLL(const std::vector<double>& yields) {
    double constraint = 0.0;

    if (yields.size() < 3) return 0.0;

    // --- N_sig constraint ---
    const double N_sig = yields[0];
    if (std::isfinite(N_sig)) {
        const double mu_sig = 0.0;
        const double sigma_sig = 1000.0;
        if (sigma_sig > 0.0 && std::isfinite(sigma_sig)) {
            const double z = (N_sig - mu_sig) / sigma_sig;
            constraint += 0.5 * z * z;
        }
    }

    // --- N_rmd constraint ---
    const double N_rmd = yields[1];
    if (std::isfinite(N_rmd)) {
        const double N_rmd_pred = 3000.0;
        const double sigma_N_rmd_pred = 100000.0;
        if (sigma_N_rmd_pred > 0.0 && std::isfinite(sigma_N_rmd_pred)) {
            const double z = (N_rmd - N_rmd_pred) / sigma_N_rmd_pred;
            constraint += 0.5 * z * z;
        }
    }

    // --- N_acc constraint ---
    const double N_acc = yields[2];
    if (std::isfinite(N_acc)) {
        const double N_acc_pred = 1332.75;
        const double sigma_N_acc_pred = 18.2534;
        if (sigma_N_acc_pred > 0.0 && std::isfinite(sigma_N_acc_pred)) {
            const double z = (N_acc - N_acc_pred) / sigma_N_acc_pred;
            constraint += 0.5 * z * z;
        }
    }

    return constraint;
}

