#ifndef P2MEG_NLLFIT_H
#define P2MEG_NLLFIT_H

#include <vector>
#include "p2meg/Event.h"
#include "p2meg/Likelihood.h"

// ============================================================
// p2MEG NLL フィット（最小化層）
//
// 役割:
//  - Likelihood::NLL(...) を最小化して yields（N_sig, N_rmd, ...）を推定する。
//  - 尤度関数と最小化アルゴリズムを分離する。
//  - パラメータ境界は設けない（負の N_sig を許す）。
// ============================================================

struct FitConfig {
    std::vector<double> start_yields; // 初期値（components と同順）
    int max_calls;                    // 最大評価回数（実装側で解釈）
    double tol;                       // 収束判定（実装側で解釈）
};

struct FitResult {
    int status;                       // 0:成功（例）、非0:失敗（実装側で定義）
    std::vector<double> yields_hat;   // 推定値（components と同順）
    std::vector<double> yields_err;   // 誤差（取れる場合のみ、取れないなら空でもよい）
    double nll_min;                   // 最小 NLL
};

FitResult FitNLL(
    const std::vector<Event>& events,
    const std::vector<PdfComponent>& components,
    const FitConfig& cfg
);

#endif // P2MEG_NLLFIT_H
