#ifndef P2MEG_LIKELIHOOD_H
#define P2MEG_LIKELIHOOD_H

#include <vector>
#include "p2meg/Event.h"

// ============================================================
// p2MEG 尤度（NLL）
//
// モデル:
//  - components で与えられる PDF 成分 p_k(x)（解析窓内で正規化済み）
//  - yields は解析窓内の期待事象数:
//      yields = { N_sig, N_rmd, (N_acc, ...) }
//    並びは components と一致させる。
//
// NLL（拡張尤度）:
//   NLL = (Σ_k N_k) - Σ_i log( Σ_k N_k * p_k(x_i) ) + ConstraintNLL(yields)
//
// 注意:
//  - N_sig が負になることを許す（境界を設けない）。
//  - ただし log の引数が非正になる点は数学的に不適なので、実装側で
//    大きいペナルティ値を返す仕様にする。
//  - p_min は実装(.cc)側に直書きする（ヘッダに Options は作らない）。
// ============================================================

// ---- PDF 成分評価 ----
// 返り値: 解析窓内で正規化された PDF 密度 p_k(ev) を返す。
//         窓外・未ロード・不正入力などは 0 を返す契約。
typedef double (*PdfEval)(const Event& ev, const void* ctx);

struct PdfComponent {
    const char* name; // "sig", "rmd", "acc" など（デバッグ用）
    PdfEval eval;     // p_k(ev)
    const void* ctx;  // eval に渡す任意の設定（NULL可）
};

// ---- 制約項 ----
// 制約が不要なときは 0 を返す実装にしておく。
// 後から制約を入れる場合は、この関数の中身を編集して NLL に加算する。
double ConstraintNLL(const std::vector<double>& yields);

// ---- NLL 本体 ----
double NLL(
    const std::vector<Event>& events,
    const std::vector<PdfComponent>& components,
    const std::vector<double>& yields
);

#endif // P2MEG_LIKELIHOOD_H
