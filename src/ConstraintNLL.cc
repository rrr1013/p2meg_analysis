#include "p2meg/Likelihood.h"

#include <cstdlib>
#include <cmath>

// ============================================================
// 制約項（NLLに加算）
//
// 制約が不要な間は常に 0 を返す。
// 後から制約を入れたくなったら、この関数の中身を編集して足し算する。
// ============================================================

// doubleonly_allpatterns に対する TSB->AW ACC 予測
//  - Eg bin: [0,10], [10,20], [20,30], [30,80] MeV
//  - fit 用時間ヒストグラムから 4D AW 内事象を除外
//  - sigma 下限: 20 ns
//  - 予測値は macros/predict_acc_from_tsb_by_eg.C の出力を固定入力として使う
static constexpr double kAccPredAw = 10.1275;
static constexpr double kAccPredAwErr = 8.75205;

static double ReadConstraintValueFromEnv(const char* name, double fallback) {
    const char* s = std::getenv(name);
    if (!s) return fallback;
    char* endptr = nullptr;
    const double v = std::strtod(s, &endptr);
    if (endptr == s || !std::isfinite(v)) return fallback;
    return v;
}

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
        const double N_rmd_pred = 0;
        const double sigma_N_rmd_pred = 10000;
        if (sigma_N_rmd_pred > 0.0 && std::isfinite(sigma_N_rmd_pred)) {
            const double z = (N_rmd - N_rmd_pred) / sigma_N_rmd_pred;
            constraint += 0.5 * z * z;
        }
    }

    // --- N_acc constraint ---
    const double N_acc = yields[2];
    if (std::isfinite(N_acc)) {
        const double N_acc_pred =
            ReadConstraintValueFromEnv("P2MEG_ACC_CONSTRAINT_MEAN", kAccPredAw);
        const double sigma_N_acc_pred =
            ReadConstraintValueFromEnv("P2MEG_ACC_CONSTRAINT_SIGMA", kAccPredAwErr);
        if (sigma_N_acc_pred > 0.0 && std::isfinite(sigma_N_acc_pred)) {
            const double z = (N_acc - N_acc_pred) / sigma_N_acc_pred;
            constraint += 0.5 * z * z;
        }
    }

    return constraint;
}
