#include "p2meg/Likelihood.h"

#include <cmath>
#include <limits>

// ============================================================
// NLL 実装（拡張尤度）
//
// NLL = (Σ_k N_k) - Σ_i log( Σ_k N_k * p_k(x_i) ) + ConstraintNLL(yields)
//
// 仕様:
//  - N_sig を負にしてよい（境界なし）
//  - ただし log の引数が非正になる点は数学的に不適なのでペナルティ
//  - p_min は数値保護（log(0)回避）としてのみ使う（物理カットではない）
// ============================================================

double NLL(const std::vector<Event>& events,
           const std::vector<PdfComponent>& components,
           const std::vector<double>& yields)
{
    // ---- 数値保護（ヘッダに出さない固定値）----
    static constexpr double p_min = 1e-300;          // log(0)回避用（物理カットではない）
    static constexpr double penalty = 1e100;         // 不適領域に返す十分大きい値

    // ---- 入力チェック ----
    if (components.empty()) return penalty;
    if (yields.size() != components.size()) return penalty;

    // ---- Σ N_k ----
    double Ntot = 0.0;
    for (double Nk : yields) {
        if (!std::isfinite(Nk)) return penalty;
        Ntot += Nk;
    }
    // 拡張尤度の Poisson 項として ΣN_k を足すため、ΣN_k が非正は不適扱いにする
    if (!(Ntot > 0.0)) return penalty;

    // ---- - Σ log( Σ N_k p_k(x_i) ) ----
    double nll = Ntot;

    for (const auto& ev : events) {
        double pi = 0.0; // Σ_k N_k * p_k(ev)

        for (std::size_t k = 0; k < components.size(); ++k) {
            const auto& comp = components[k];
            const double Nk = yields[k];

            // PDF値
            double pk = comp.eval ? comp.eval(ev, comp.ctx) : 0.0;
            if (!std::isfinite(pk) || pk < 0.0) pk = 0.0;

            pi += Nk * pk;
        }

        // log の定義域：pi > 0 が必要
        if (!(pi > 0.0) || !std::isfinite(pi)) return penalty;

        // 数値保護（物理カットではない）
        if (pi < p_min) pi = p_min;

        nll -= std::log(pi);
    }

    // ---- 制約項 ----
    const double c = ConstraintNLL(yields);
    if (!std::isfinite(c)) return penalty;
    nll += c;

    return nll;
}
