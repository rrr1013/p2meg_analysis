#ifndef P2MEG_UPPER_LIMIT_H
#define P2MEG_UPPER_LIMIT_H

#include <cstddef>
#include <vector>

#include "p2meg/Event.h"
#include "p2meg/Likelihood.h"
#include "p2meg/NLLFit.h"

// ============================================================
// p2MEG 90% C.L. upper limit 用ユーティリティ
//
// 目的:
//  - 既存の extended likelihood / PdfComponent / FitNLL を流用し、
//    μ→eγ の signal yield N_sig に対する profile-likelihood 型の
//    toy MC upper limit を追加する。
//  - 既存の通常 fit は変更せず、upper limit 専用の薄い追加層として使う。
//
// 今回の最小実装での近似:
//  - toy 生成は最終解析の PDF から直接サンプリングする。
//  - 事象数は各成分ごとに Poisson fluctuation を入れる。
//  - nuisance parameter は yields のみとし、
//      {N_sig, N_rmd, N_acc, ...}
//    のうち N_sig を固定、他を profile fit で再最適化する。
//  - sideband / normalisation 由来の追加 nuisance は、
//    現在のコードに明示的な実装が無い限り新設しない。
//    従って、今回含まれる制約は NLL 内の ConstraintNLL(yields) のみ。
//
// 単位:
//  - Event の定義に従う（Ee, Eg: MeV, t: ns, phi: rad）
//  - N_mu_eff は「有効停止ミューオン数」で無次元
// ============================================================

struct ToyGeneratorConfig {
    unsigned long long seed;  // 乱数 seed
    int pmax_scan_trials;     // 棄却法 pmax 推定の試行数
    double pmax_safety;       // pmax 安全係数
    double pmax_update;       // 生成中に pmax 不足時の更新係数
    int event_pool_size_per_component; // 成分ごとに先に作る 5D event pool のサイズ
};

struct NormalizationUncertaintyConfig {
    double N_mu_eff_nom;    // 有効停止ミューオン数の公称値
    double N_mu_eff_sigma;  // 有効停止ミューオン数の絶対誤差
};

struct UpperLimitPointConfig {
    double N_sig_test;        // 固定する signal yield 仮説値
    int n_toys;               // toy 本数
    double cl;                // 例: 0.90
    FitConfig free_fit_cfg;   // 通常 fit 用初期値
    FitConfig prof_fit_cfg;   // 固定 N_sig profile fit 用初期値
    ToyGeneratorConfig toy_cfg;
};

struct UpperLimitPointResult {
    double N_sig_test;            // 入力した仮説値 s
    double q_obs;                 // 実データの q(s) = -2 ln lambda(s)
    double p_value;               // toy 分布に対する右側 p-value
    double acceptance_threshold;  // 1 - CL
    bool accepted;                // p_value >= 1-CL なら true
    int n_toys_requested;         // 要求 toy 数
    int n_toys_valid;             // q_toy が有効に得られた toy 数
    FitResult fit_free_obs;       // 実データ free fit
    FitResult fit_prof_obs;       // 実データ fixed-N_sig profile fit
};

struct UpperLimitScanConfig {
    std::vector<double> N_sig_scan; // 走査する signal yield 値
    int n_toys_per_point;           // 各点の toy 数
    double cl;                      // 例: 0.90
    double N_mu_eff;                // BR_90 = N_sig^90 / N_mu_eff
    FitConfig free_fit_cfg;         // 通常 fit 用初期値
    FitConfig prof_fit_cfg;         // 固定 N_sig profile fit 用初期値
    ToyGeneratorConfig toy_cfg;
};

struct UpperLimitScanResult {
    FitResult fit_free_obs;                    // 実データ通常 fit
    std::vector<UpperLimitPointResult> points; // 各 s 点の結果
    double N_sig_90;                           // 90% C.L. upper limit on N_sig
    double BR_90;                              // N_sig_90 / N_mu_eff
    double N_mu_eff;                           // 入力した有効停止ミューオン数
};

struct UpperLimitBRPointConfig {
    double BR_test;                      // 固定する分岐比仮説値
    int n_toys;                          // toy 本数
    double cl;                           // 例: 0.90
    FitConfig free_fit_cfg;              // 通常 fit 用初期値
    FitConfig prof_fit_cfg;              // 固定 BR（=固定 N_sig_nom）profile fit 用初期値
    ToyGeneratorConfig toy_cfg;
    NormalizationUncertaintyConfig norm_cfg;
};

struct UpperLimitBRPointResult {
    double BR_test;                 // 入力した分岐比仮説値
    double N_sig_test_nominal;      // 公称 N_mu_eff に対応する signal yield 仮説値
    double q_obs;                   // 実データの q(BR)
    double p_value;                 // toy 分布に対する右側 p-value
    double acceptance_threshold;    // 1 - CL
    bool accepted;                  // p_value >= 1-CL なら true
    int n_toys_requested;           // 要求 toy 数
    int n_toys_valid;               // q_toy が有効に得られた toy 数
    FitResult fit_free_obs;         // 実データ free fit
    FitResult fit_prof_obs;         // 実データ fixed-BR profile fit
};

struct UpperLimitBRScanConfig {
    std::vector<double> BR_scan;      // 走査する分岐比
    int n_toys_per_point;             // 各点の toy 数
    double cl;                        // 例: 0.90
    FitConfig free_fit_cfg;           // 通常 fit 用初期値
    FitConfig prof_fit_cfg;           // 固定 BR profile fit 用初期値
    ToyGeneratorConfig toy_cfg;
    NormalizationUncertaintyConfig norm_cfg;
};

struct UpperLimitBRScanResult {
    FitResult fit_free_obs;                        // 実データ通常 fit
    std::vector<UpperLimitBRPointResult> points;  // 各 BR 点の結果
    double BR_90;                                 // 90% C.L. upper limit on BR
    double N_sig_90_nominal;                      // BR_90 * N_mu_eff_nom
    NormalizationUncertaintyConfig norm_cfg;      // 入力した normalisation 設定
};

struct ProfileLikelihoodQPoint {
    double N_sig_test;     // 固定した signal yield 仮説値
    double q_value;        // q = -2 ln lambda
    FitResult fit_prof;    // fixed-N_sig profile fit
};

struct ProfileLikelihoodQScanResult {
    FitResult fit_free;                        // 同一 dataset に対する free fit
    std::vector<ProfileLikelihoodQPoint> points; // 各固定仮説値の結果
};

// N_sig を固定した profile fit
//  - yields[0] = N_sig_fixed を固定し、残り成分を最小化する
//  - 入力不正や最小化失敗時は status!=0 を返す
FitResult FitNLLFixedSignal(
    const std::vector<Event>& events,
    const std::vector<PdfComponent>& components,
    const FitConfig& cfg,
    double N_sig_fixed
);

// q(s) = -2 ln lambda(s) = 2 [ NLL_prof(s) - NLL_best ]
//  - free fit / profile fit の両方を返す
//  - 数値誤差で負になった場合は 0 に丸める
double EvaluateProfileLikelihoodQ(
    const std::vector<Event>& events,
    const std::vector<PdfComponent>& components,
    const FitConfig& free_fit_cfg,
    const FitConfig& prof_fit_cfg,
    double N_sig_fixed,
    FitResult& fit_free_out,
    FitResult& fit_prof_out
);

// 同一 dataset に対して複数の N_sig_fixed をまとめて評価する
//  - free fit は 1 回だけ行い、各 fixed-N_sig 点では cache を再利用する
//  - 統計量 q の定義は EvaluateProfileLikelihoodQ と同じ
ProfileLikelihoodQScanResult EvaluateProfileLikelihoodQScan(
    const std::vector<Event>& events,
    const std::vector<PdfComponent>& components,
    const FitConfig& free_fit_cfg,
    const FitConfig& prof_fit_cfg,
    const std::vector<double>& N_sig_scan
);

// PDF から toy dataset を 1 本生成する
//  - mean_yields[k] を平均とする Poisson fluctuation を各成分に入れる
//  - out_generated_yields を与えた場合、実際に生成した各成分個数を返す
//  - 失敗時は false
bool GenerateToyDatasetFromModel(
    const std::vector<PdfComponent>& components,
    const std::vector<double>& mean_yields,
    const ToyGeneratorConfig& cfg,
    unsigned long long toy_index,
    std::vector<Event>& out_events,
    std::vector<double>* io_pmax_cache = nullptr,
    std::vector<double>* out_generated_yields = nullptr
);

// 単一の s 点について toy MC による受容判定を行う
//  - toy 生成モデルの nuisance は「実データに対する fixed-s profile fit」の値を使う
//    （plug-in / profile construction の最小版）
UpperLimitPointResult EvaluateUpperLimitPoint(
    const std::vector<Event>& events,
    const std::vector<PdfComponent>& components,
    const UpperLimitPointConfig& cfg
);

// s を走査して 90% C.L. upper limit を返す
//  - accepted な点のうち最大の N_sig を N_sig_90 とする
//  - N_mu_eff > 0 のとき BR_90 も返す。不正値なら BR_90=0
UpperLimitScanResult EvaluateUpperLimitScan(
    const std::vector<Event>& events,
    const std::vector<PdfComponent>& components,
    const UpperLimitScanConfig& cfg
);

// 単一の BR 点について toy MC による受容判定を行う
//  - q(BR) の固定 signal yield は BR * N_mu_eff_nom を使う
//  - toy 生成時のみ N_mu_eff を Gaussian で揺らして normalisation uncertainty を入れる
UpperLimitBRPointResult EvaluateUpperLimitBRPoint(
    const std::vector<Event>& events,
    const std::vector<PdfComponent>& components,
    const UpperLimitBRPointConfig& cfg
);

// BR を走査して 90% C.L. upper limit を返す
//  - accepted な点のうち最大の BR を BR_90 とする
UpperLimitBRScanResult EvaluateUpperLimitBRScan(
    const std::vector<Event>& events,
    const std::vector<PdfComponent>& components,
    const UpperLimitBRScanConfig& cfg
);

#endif // P2MEG_UPPER_LIMIT_H
