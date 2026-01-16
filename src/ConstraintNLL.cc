#include "p2meg/Likelihood.h"

#include <cmath>

// ============================================================
// 制約項（NLLに加算）
//
// 制約が不要な間は常に 0 を返す。
// 後から制約を入れたくなったら、この関数の中身を編集して足し算する。
// ============================================================

double ConstraintNLL(const std::vector<double>& yields) {

    // N_sig 安定化のための弱いガウス制約（データが支配的ならほとんど無視される）

    if (yields.empty()) return 0.0;

    const double N_sig = yields[0];
    if (!std::isfinite(N_sig)) return 0.0;

    const double mu_sig = 0.0;
    const double sigma_sig = 1000.0; // 弱い制約（大きめ）
    if (!(sigma_sig > 0.0) || !std::isfinite(sigma_sig)) return 0.0;

    const double z = (N_sig - mu_sig) / sigma_sig;
    return 0.5 * z * z;
}
