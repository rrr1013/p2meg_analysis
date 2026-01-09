#include "p2meg/NLLFit.h"

#include <cmath>
#include <memory>

// ROOT Minimizer
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"

// ============================================================
// FitNLL 実装（Minuit2）
//
// 境界なし（負の N_sig も許す）
// 56行のコメントアウトを消すことで負のNを禁止できる
// ただし NLL が不適領域（logの引数<=0など）に入ると大ペナルティを返すので、
// 最小化は自然に定義域内へ戻ることを期待する。
// ============================================================

FitResult FitNLL(const std::vector<Event>& events,
                 const std::vector<PdfComponent>& components,
                 const FitConfig& cfg)
{
    FitResult out;
    out.status = 1;
    out.nll_min = std::numeric_limits<double>::quiet_NaN();

    const std::size_t npar = components.size();
    if (npar == 0) return out;

    if (cfg.start_yields.size() != npar) {
        // 初期値が合わない場合は失敗扱い
        return out;
    }

    // Minuit2 を作る
    std::unique_ptr<ROOT::Math::Minimizer> min(
        ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad")
    );
    if (!min) return out;

    if (cfg.max_calls > 0) min->SetMaxFunctionCalls(cfg.max_calls);
    if (cfg.tol > 0.0)      min->SetTolerance(cfg.tol);

    // NLL を評価する関数（yields を渡して NLL(...) を呼ぶだけ）
    auto fcn = [&](const double* x) -> double {
        std::vector<double> yields(npar);
        for (std::size_t i = 0; i < npar; ++i) yields[i] = x[i];
        return NLL(events, components, yields);
    };

    ROOT::Math::Functor functor(fcn, static_cast<unsigned int>(npar));
    min->SetFunction(functor);

    // N_sig > 0, N_rmd > 0 の制約（下の 1 行のコメントを外すと有効化）
     //#define P2MEG_ENABLE_YIELD_BOUNDS

    // パラメータ設定（境界なし）
    for (std::size_t i = 0; i < npar; ++i) {
        const double start = cfg.start_yields[i];
        const double step  = (std::abs(start) > 0.0) ? 0.1 * std::abs(start) : 1.0;
#ifdef P2MEG_ENABLE_YIELD_BOUNDS
        // yields = {N_sig, N_rmd, ...} を仮定（先頭2つ）
        if (i == 0 || i == 1) {
            min->SetLowerLimitedVariable(static_cast<unsigned int>(i),
                                         components[i].name, start, step, 0.0);
            continue;
        }
#endif
        min->SetVariable(static_cast<unsigned int>(i), components[i].name, start, step);
    }

    const bool ok = min->Minimize();
    out.status = ok ? 0 : 2;

    const double* xs = min->X();
    out.yields_hat.assign(xs, xs + npar);
    out.nll_min = min->MinValue();

    // 誤差（取れない場合もあるので、その場合は空にする）
    out.yields_err.clear();
    const double* es = min->Errors();
    if (es) {
        out.yields_err.assign(es, es + npar);
    }

    return out;
}
