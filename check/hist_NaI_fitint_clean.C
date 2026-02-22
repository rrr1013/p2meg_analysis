// check/hist_NaI_fitint_clean.C
//
// NaI (A1,A2,B1,B2) を独立に見て、clean を優先して
//  - baseline: 1 digitizer event 内の ADC 最頻値(mode)
//  - パルス検出: 候補抽出 → フィット（最小モデル）
//  - 大きさ: フィット関数の面積（= 積分値）
// を求め、積分値分布ヒストを PDF 保存する。
//
// テンプレは実データ駆動ではなく関数形（差の指数）を使用。
// τ は clean サンプルから推定して固定する。
//
// さらに確認用として、PDF 最後に
//  - 先頭から「受理された clean パルスが含まれる digitizer event」5個
//  - 末尾から同様に 5個
// について、波形・フィット曲線・パラメータを 2x2 で可視化したページを追加する。
//
// 出力: doc/mainexp/hist_NaI_fitint_clean_<入力データ名>.pdf
//
// 実行例（リポジトリ直下で）:
//   root -l -q 'check/hist_NaI_fitint_clean.C'

#include <iostream>
#include <fstream>
#include <vector>
#include <deque>
#include <algorithm>
#include <cmath>
#include <limits>

#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TDatime.h"
#include "TStyle.h"
#include "TString.h"
#include "TSystem.h"
#include "TF1.h"
#include "TFile.h"

// ============================================================
// 手で変えるパラメータはここだけ
// ============================================================

// 入力ディレクトリ（例: data/rawdata/mainexp/60°/run17）
static const char* kInputDir = "data/rawdata/mainexp/120°/run8";

// 入力データ名（出力ファイル名に使う） 例: "60°_run17"
static const char* kInputTag = "120°_run8";

// 入力ファイル（NaI 4ch）
static const char* kNaIFiles[4] = {
  "wave_NaI_A1_run8.txt",
  "wave_NaI_A2_run8.txt",
  "wave_NaI_B1_run8.txt",
  "wave_NaI_B2_run8.txt"
};

// 出力ディレクトリ（固定）
static const char* kOutputDir = "doc/mainexp";

// 波形→イベント切り出し
static const int    kSamplesPerEvent = 1000;
static const double kDtNs = 4.0;

// clean 優先：フィット窓（前後150bin）
static const int kFitHalfWindow = 150;
static const int kFitPreWindow  = 40;   // フィットに使う pre-trigger 幅（samples）

// 候補抽出（暫定、ヒストを見て調整する前提）
static const double kSeedThr = 60.0;  // s = baseline - ADC がこれ以上で「開始候補」
static const double kEndThr  = 30.0;  // ヒステリシス用（候補列の整形に使う場合もある）
static const int    kMinSepSamples = 100; // 候補開始点の最小分離 [samples]

// clean 判定（候補が fit 窓内で単発）
static const bool   kCleanRequireSingleInWindow = false;

// τ 推定（clean サンプルだけで推定して固定）
static const int    kNTauEstTarget = 300; // chごとにこの本数まで平均波形に足す
static const int    kNTauEstMin    = 30;  // これ未満なら τ 推定を諦めてデフォルトへ

// τ 推定フィット範囲（相対時間 dt = [-150, +150]）
static const int    kTauEstHalfWindow = 150;

// τ 初期値（samples 単位）
static const double kTauR_init = 3.0;   // 立ち上がり時定数（samples）
static const double kTauD_init = 30.0;  // 差分 d = tau_f - tau_r（samples）※面積係数になる
static const double kTauR_min_phys = 1.0; // 物理ガード（最小立ち上がり時定数）
static const double kTauD_min_phys = 5.0; // 物理ガード（最小減衰差）
static const double kTauBoundaryEps = 1.0e-3; // 境界張り付き判定

// per-pulse fit の t0 探索幅（初期値の周り）
static const double kT0FitRange = 30.0; // [samples]
static const bool   kAdaptiveT0Range = true;
static const double kT0FitRangeExpandFactor = 2.0; // 端張り付き時に探索幅を拡大
static const double kT0BoundaryFracExpand = 0.90;  // 拡大判定
static const double kT0BoundaryFracReject = 0.98; // 最終的に端に張り付く解を棄却
static const double kT0FitAbsMin = 15.0;          // event先頭近傍の不安定解を棄却
static const double kOffsetRetryChi2 = 20.0;      // これ以下なら offset 再試行しない

// フィットの baseline（s空間オフセット）を 0 に固定する
// s = baseline(mode) - ADC を使っているので、物理的には C=0 を期待する
static const bool   kFixFitBaselineToZero = true;

// chi2/ndf カット
// - まず run 冒頭のサンプルで chi2 分布を走査し、チャンネルごとに自動決定
// - 自動決定が不可能な場合のみ fallback 値を使う
static const bool   kAutoTuneChi2Cut = true;
static const int    kChi2ScanMaxEvents = 5000;       // 走査に使う最大 event 数
static const int    kChi2AutoMinFits   = 30;         // 自動決定に必要な最小 fit 成功数
static const int    kChi2ScanTargetFitsPerCh = 300;  // この本数に達したchは十分とみなす
static const double kChi2KeepFraction  = 0.95;       // 良い側から残す割合（0.95 -> 95%）
static const double kChi2NdfFallback   = 120.0;      // fallback
static const double kChi2NdfMinClamp   = 5.0;        // 下限クランプ
static const double kChi2NdfMaxClamp   = 120.0;      // 上限クランプ
static const double kRelShapeSys       = 0.02;       // 形状系統（chi2用の相対誤差）

// 候補が複数あるときの扱い
static const bool   kUsePrimaryPulseOnly = false; // true: 1event/ch で最大面積1本のみヒストに入れる
static const bool   kDebugShowLargestArea = true; // debug表示は最大面積候補を優先
static const int    kMaxCandidatesToFitPrimary = 1; // primary-only時にfitする候補数上限

// 候補の事前品質（高速フィルタ）
// - pre-trigger が静穏でない候補を除外（前段パルス尾の誤フィットを減らす）
// - fit窓内に強い第2ピークがある候補を除外（重なりパルスをclean扱いしない）
static const int    kCandPeakLookahead      = 220;   // 候補スコアのピーク探索幅 [samples]
static const int    kCandPreStart           = 80;    // pre静穏判定の開始オフセット [samples]
static const int    kCandPreEnd             = 12;    // pre静穏判定の終了オフセット [samples]
static const int    kCandPreMinSamples      = 20;    // pre静穏判定に必要な最小サンプル数
static const double kCandPreQuietAbsMax     = 50.0;  // pre静穏判定の絶対上限 [counts]
static const double kCandPreQuietFracMax    = 0.20;  // pre静穏判定の相対上限（peak比）
static const bool   kCandRejectMultiPeak    = true;  // fit窓内の多峰候補を除外する
static const int    kCandSecondPeakMinSep   = 30;    // 主峰からこの距離以上を第2峰候補とする [samples]
static const double kCandSecondPeakFracMax  = 0.75;  // 第2峰/主峰 がこれを超えたら除外
static const bool   kUsePeakSupplementCand  = true;  // 閾値cross以外に局所ピーク由来候補も使う
static const double kPeakCandMinHeight      = 70.0;  // 局所ピーク候補の最小高さ [counts]
static const double kPeakCandMinProminence  = 15.0;  // 局所ピークの最小プロミネンス [counts]
static const int    kPeakCandMinSep         = 80;    // 局所ピーク候補の最小分離 [samples]
static const int    kPeakCandBackSearch     = 220;   // peak->t0 逆探索幅 [samples]
static const int    kPeakCandFallbackShift  = 24;    // 逆探索失敗時の t0 = peak-shift [samples]
static const int    kCandMergeMinSep        = 10;    // 候補t0マージの最小分離 [samples]
static const bool   kRequirePrimaryPeakConsistency = true; // primary-only時、主ピーク整合を要求
static const double kPrimaryPeakFracMin = 0.75;            // fit対象ピーク / event最大ピーク の下限

// 多パルス分離（重なり波形を順次フィットして引く）
static const bool   kEnablePulseDecomposition = true;
static const int    kMaxAcceptedPerEventCh = 6;
static const int    kMaxTriedPerEventCh = 24;
static const double kSubtractTailFactor = 8.0; // 減衰時定数の何倍まで引くか
static const double kDecompResidualStopThr = 55.0; // 残差の最大ピークがこれ未満なら分離停止
static const double kDecompMinModelPeakAbs = 60.0; // 受理する成分の最小モデルピーク [counts]
static const double kDecompMinModelPeakSigma = 6.0; // 受理する成分の最小S/N（peak/sigma）
static const double kDecompMinDeltaChi2 = 40.0; // 成分受理に必要な最小 Δχ²
static const double kDecompMinRelImprove = 0.12; // 成分受理に必要な最小相対改善率
static const int    kDecompMinT0Sep = 28; // 既受理パルスとの最小分離[samples]
static const int    kJointRefineIters = 2; // 多成分同時フィットの t0 交互最適化回数
static const int    kJointT0ScanRange = 12; // t0 走査半幅 [samples]
static const int    kJointT0ScanStep = 1;   // t0 走査刻み [samples]
static const int    kJointT0MinSepHard = 6; // t0 交互最適化中の最小分離 [samples]
static const int    kPeakAssignMaxDist = 36; // 成分と実ピークの対応許容距離 [samples]
static const int    kDataEvPre = 24;      // 元波形での証拠確認: t0 前の探索幅
static const int    kDataEvPost = 40;     // 元波形での証拠確認: t0 後の探索幅
static const double kDataEvPeakMin = 75.0; // 元波形の局所ピーク最小
static const double kDataEvPromMin = 24.0; // 元波形の局所prominence最小
static const double kDataEvRiseMin = 8.0;  // 元波形の最小立ち上がり差分

// 積分値ヒスト設定（単位: ADC*sample）
// ※ raw の (baseline-ADC) 総和と同じ「サンプル積分」単位に合わせる
static const int    kNBinsI = 2500;
static const double kIMin   = 0.0;
static const double kIMax   = 2e5;

// 描画
static const bool   kUseLogY = true;
static const int    kCanvasW = 1200;
static const int    kCanvasH = 900;

// デバッグ波形（PDF最後に追加）
static const int    kNDebugFirst = 5;
static const int    kNDebugLast  = 5;

// 形状診断（Michel由来の山/肩を簡易に見るための窓）
static const double kShapeW1Min = 2.0e4;  // 低I
static const double kShapeW1Max = 6.0e4;
static const double kShapeW2Min = 6.0e4;  // 山候補
static const double kShapeW2Max = 1.2e5;
static const double kShapeW3Min = 1.2e5;  // 高I tail
static const double kShapeW3Max = 1.8e5;

// ラベル
static const char* kChLabel[4] = {"NaI_A1","NaI_A2","NaI_B1","NaI_B2"};

// ============================================================

struct FitPulse {
  int    t0_init = -1;
  double t0_fit  = 0.0;
  double A       = 0.0;  // モデルの振幅係数
  double C       = 0.0;  // オフセット（s空間）
  double pre_offset = 0.0; // フィット前に引いた局所 pre-trigger オフセット
  double I       = 0.0;  // 積分値（ADC*sample）= A * d (d = tau_f - tau_r)
  double chi2ndf = 0.0;
  bool   ok      = false;
};

struct DebugEvent {
  long idx = -1;
  double baseline[4] = {0,0,0,0}; // ADC 空間
  std::vector<double> adc[4];     // 元ADC（長さ kSamplesPerEvent）
  FitPulse fit[4];                // 代表1パルス（注釈用）
  std::vector<FitPulse> fit_hist[4]; // 実際にヒストへ入れたパルス集合
  int n_cand[4] = {0,0,0,0};
  int n_acc[4]  = {0,0,0,0};      // fit受理数（候補選別後）
  int n_hist[4] = {0,0,0,0};      // ヒスト投入数
};

// 候補選別/fit 呼び出しの診断カウンタ（速度・品質確認用）
static long gCandTotalInput = 0;
static long gCandRejectEdge = 0;
static long gCandRejectClean = 0;
static long gCandRejectPreQuiet = 0;
static long gCandRejectMultiPeak = 0;
static long gCandAcceptedForFit = 0;
static long gFitCalls = 0;

struct ShapeMetric {
  double n_w1 = 0.0;
  double n_w2 = 0.0;
  double n_w3 = 0.0;
  double bump_ratio = 0.0; // >1 なら中域の山/肩が相対的に強い
};

static ShapeMetric ComputeShapeMetric(const TH1D* h)
{
  ShapeMetric m;
  if (!h) return m;
  const int b1a = h->GetXaxis()->FindBin(kShapeW1Min);
  const int b1b = h->GetXaxis()->FindBin(kShapeW1Max);
  const int b2a = h->GetXaxis()->FindBin(kShapeW2Min);
  const int b2b = h->GetXaxis()->FindBin(kShapeW2Max);
  const int b3a = h->GetXaxis()->FindBin(kShapeW3Min);
  const int b3b = h->GetXaxis()->FindBin(kShapeW3Max);
  m.n_w1 = h->Integral(std::min(b1a, b1b), std::max(b1a, b1b));
  m.n_w2 = h->Integral(std::min(b2a, b2b), std::max(b2a, b2b));
  m.n_w3 = h->Integral(std::min(b3a, b3b), std::max(b3a, b3b));
  const double denom = std::sqrt(std::max(1.0, m.n_w1 * m.n_w3));
  m.bump_ratio = m.n_w2 / denom;
  return m;
}

// ---------- baseline: event内 mode ----------
static double ModeBaselineInEvent(const std::vector<double>& samples)
{
  if (samples.empty()) return 0.0;

  int vmin = (int)std::lround(samples[0]);
  int vmax = vmin;
  for (double x : samples) {
    int v = (int)std::lround(x);
    if (v < vmin) vmin = v;
    if (v > vmax) vmax = v;
  }

  const int range = vmax - vmin + 1;
  if (range <= 0) return (double)vmin;

  std::vector<int> cnt(range, 0);
  for (double x : samples) {
    int v = (int)std::lround(x);
    const int idx = v - vmin;
    if (0 <= idx && idx < range) cnt[idx]++;
  }

  int best_i = 0;
  int best_c = cnt[0];
  for (int i = 1; i < range; ++i) {
    if (cnt[i] > best_c) { best_c = cnt[i]; best_i = i; }
  }

  return (double)(vmin + best_i);
}

// ---------- s = baseline - ADC を作る ----------
static void BuildS(const std::vector<double>& adc, double baseline, std::vector<double>& s_out)
{
  s_out.resize(adc.size());
  for (size_t i = 0; i < adc.size(); ++i) {
    s_out[i] = baseline - adc[i];
  }
}

// ---------- 候補開始点（threshold crossing）を抽出 ----------
static std::vector<int> FindStartCandidates(const std::vector<double>& s,
                                            double seed_thr,
                                            double end_thr,
                                            int min_sep_samples)
{
  std::vector<int> starts;
  const int n = (int)s.size();
  if (n <= 1) return starts;

  int dead = 0;
  for (int i = 1; i < n; ++i) {
    if (dead > 0) { dead--; continue; }

    const bool crossed = (s[i-1] < seed_thr && s[i] >= seed_thr);
    if (crossed) {
      starts.push_back(i);

      // ヒステリシス:
      // いったん end_thr を下回るまで次候補を作らない（ノイズの多重トリガ回避）
      int j = i + 1;
      while (j < n && s[j] >= end_thr) j++;

      const int width_dead = std::max(0, j - i);
      dead = std::max(min_sep_samples, width_dead);
    }
  }
  return starts;
}

// ---------- clean: fit窓内で単発か判定 ----------
static bool IsSingleCandidateInWindow(const std::vector<int>& starts,
                                      int t0,
                                      int halfwin)
{
  int n_in = 0;
  for (int s0 : starts) {
    if (std::abs(s0 - t0) <= halfwin) n_in++;
    if (n_in >= 2) return false;
  }
  return (n_in == 1);
}

// 候補 t0 の「主パルスらしさ」スコア（t0 以降の局所最大）
static double CandidatePeakScore(const std::vector<double>& s, int t0)
{
  const int n = (int)s.size();
  if (t0 < 0 || t0 >= n) return 0.0;
  const int i1 = std::min(n - 1, t0 + kCandPeakLookahead);
  double smax = 0.0;
  for (int i = t0; i <= i1; ++i) {
    if (s[i] > smax) smax = s[i];
  }
  return smax;
}

// 候補直前の pre-trigger 最大値（静穏性チェック用）
// 戻り値:
//   true  : pre区間が確保でき、pre_max が有効
//   false : pre区間不足（event先頭近傍）で静穏性を判定できない
static bool CandidatePreMax(const std::vector<double>& s, int t0, double& pre_max)
{
  pre_max = 0.0;
  const int n = (int)s.size();
  if (t0 < 0 || t0 >= n) return false;

  const int i0 = std::max(0, t0 - kCandPreStart);
  const int i1 = std::min(n - 1, t0 - kCandPreEnd);
  if (i1 < i0) return false;
  if ((i1 - i0 + 1) < kCandPreMinSamples) return false;

  double vmax = 0.0;
  for (int i = i0; i <= i1; ++i) {
    if (s[i] > vmax) vmax = s[i];
  }
  pre_max = vmax;
  return true;
}

// fit窓内に強い第2ピークがあるかを判定
// 物理意図:
//   単一指数差モデルは「単発パルス」を仮定しているため、
//   同一fit窓内で主峰と同程度の第2峰がある候補は clean と見なさない。
static bool HasStrongSecondaryPeak(const std::vector<double>& s,
                                   int t0,
                                   int fit_post,
                                   double primary_peak)
{
  if (primary_peak <= 0.0) return false;

  const int n = (int)s.size();
  const int i0 = std::max(1, t0 + 2);
  const int i1 = std::min(n - 2, t0 + fit_post);
  if (i1 <= i0) return false;

  int i_primary = -1;
  double p_primary = 0.0;
  std::vector<std::pair<int, double>> peaks;
  peaks.reserve(16);
  for (int i = i0; i <= i1; ++i) {
    if (s[i] < kSeedThr) continue;
    if (!(s[i] >= s[i - 1] && s[i] > s[i + 1])) continue;
    peaks.push_back({i, s[i]});
    if (s[i] > p_primary) {
      p_primary = s[i];
      i_primary = i;
    }
  }
  if (i_primary < 0) return false;

  double p_second = 0.0;
  for (const auto& pk : peaks) {
    const int ip = pk.first;
    if (std::abs(ip - i_primary) < kCandSecondPeakMinSep) continue;
    if (pk.second > p_second) p_second = pk.second;
  }

  const double ref_peak = std::max(primary_peak, p_primary);
  if (ref_peak <= 0.0) return false;
  return (p_second > kCandSecondPeakFracMax * ref_peak);
}

// primary-only の場合、fit対象がイベント主ピークを表しているか確認する
static bool PassPrimaryPeakConsistency(const std::vector<double>& s, double t0_fit)
{
  if (!kUsePrimaryPulseOnly) return true;
  if (!kRequirePrimaryPeakConsistency) return true;
  if (s.empty()) return false;

  double s_global_peak = 0.0;
  for (double y : s) if (y > s_global_peak) s_global_peak = y;
  if (s_global_peak <= 0.0) return false;

  const int t0i = (int)std::lround(t0_fit);
  const double s_fit_peak = CandidatePeakScore(s, t0i);
  return (s_fit_peak >= kPrimaryPeakFracMin * s_global_peak);
}

// 候補時刻を昇順にして、近すぎる重複を統合する
static void MergeCandidateTimes(std::vector<int>& t0s, int min_sep)
{
  if (t0s.empty()) return;
  std::sort(t0s.begin(), t0s.end());

  std::vector<int> out;
  out.reserve(t0s.size());
  out.push_back(t0s[0]);
  for (size_t i = 1; i < t0s.size(); ++i) {
    if (t0s[i] - out.back() >= min_sep) out.push_back(t0s[i]);
  }
  t0s.swap(out);
}

// 局所ピーク由来の補助候補を追加する
// 目的:
//   閾値cross法だけでは、基線が高いまま立ち上がる主パルスを取りこぼすため、
//   局所ピークから rise 開始点を逆探索して t0 候補を補う。
static void AddPeakSupplementCandidates(const std::vector<double>& s,
                                        std::vector<int>& starts_io)
{
  if (!kUsePeakSupplementCand) return;
  const int n = (int)s.size();
  if (n < 5) return;

  std::vector<int> t0_extra;
  t0_extra.reserve(32);
  int last_peak = -kPeakCandMinSep;
  for (int i = 2; i <= n - 3; ++i) {
    if (i - last_peak < kPeakCandMinSep) continue;
    const double si = s[i];
    if (si < kPeakCandMinHeight) continue;

    if (!(si >= s[i - 1] && si > s[i + 1])) continue;

    const double shoulder = std::max(s[i - 2], s[i + 2]);
    const double prom = si - shoulder;
    if (prom < kPeakCandMinProminence) continue;

    int t0 = i - kPeakCandFallbackShift;
    if (t0 < 1) t0 = 1;
    const int jmin = std::max(1, i - kPeakCandBackSearch);
    for (int j = i; j > jmin; --j) {
      if (s[j - 1] < kEndThr && s[j] >= kEndThr) {
        t0 = j;
        break;
      }
    }

    t0_extra.push_back(t0);
    last_peak = i;
  }

  if (t0_extra.empty()) return;
  starts_io.insert(starts_io.end(), t0_extra.begin(), t0_extra.end());
  MergeCandidateTimes(starts_io, kCandMergeMinSep);
}

// イベント最大ピーク由来の候補を必ず1本追加する
// 目的:
//   閾値cross/局所ピーク補助でも主ピークが候補化できないケースを救済する。
static void AddGlobalPeakCandidate(const std::vector<double>& s,
                                   std::vector<int>& starts_io)
{
  if (s.empty()) return;
  int i_peak = -1;
  double s_peak = 0.0;
  for (int i = 0; i < (int)s.size(); ++i) {
    if (s[i] > s_peak) {
      s_peak = s[i];
      i_peak = i;
    }
  }
  if (i_peak < 0 || s_peak < kSeedThr) return;

  int t0 = std::max(1, i_peak - kPeakCandFallbackShift);
  const int jmin = std::max(1, i_peak - kPeakCandBackSearch);
  for (int j = i_peak; j > jmin; --j) {
    if (s[j - 1] < kEndThr && s[j] >= kEndThr) {
      t0 = j;
      break;
    }
  }
  starts_io.push_back(t0);
  MergeCandidateTimes(starts_io, kCandMergeMinSep);
}

// fit 対象候補を作る:
// - 幾何/clean 条件を満たす候補を集める
// - primary-only ならピーク順に上位だけ残す
static void BuildFitCandidates(const std::vector<int>& starts,
                               const std::vector<double>& s,
                               int fit_pre,
                               int fit_post,
                               std::vector<int>& out_t0)
{
  out_t0.clear();

  struct Cand {
    int t0;
    double score;
  };
  std::vector<Cand> tmp;
  tmp.reserve(starts.size());

  const int n = (int)s.size();
  double s_global_peak = 0.0;
  for (double y : s) if (y > s_global_peak) s_global_peak = y;
  const bool strict_prefilter = kUsePrimaryPulseOnly;
  gCandTotalInput += (long)starts.size();
  for (int t0 : starts) {
    if (t0 - fit_pre < 0 || t0 + fit_post >= n) {
      gCandRejectEdge++;
      continue;
    }
    if (kCleanRequireSingleInWindow) {
      if (!IsSingleCandidateInWindow(starts, t0, fit_post)) {
        gCandRejectClean++;
        continue;
      }
    }

    const double s_peak = CandidatePeakScore(s, t0);
    if (s_peak < kSeedThr) continue;
    const bool is_global_like = (s_global_peak > 0.0 && s_peak >= 0.85 * s_global_peak);

    double pre_max = 0.0;
    if (strict_prefilter) {
      if (!CandidatePreMax(s, t0, pre_max)) {
        gCandRejectPreQuiet++;
        continue;
      }
      const double pre_allow = std::max(kCandPreQuietAbsMax, kCandPreQuietFracMax * s_peak);
      if (!is_global_like && pre_max > pre_allow) {
        gCandRejectPreQuiet++;
        continue;
      }

      if (kCandRejectMultiPeak) {
        if (!is_global_like && HasStrongSecondaryPeak(s, t0, fit_post, s_peak)) {
          gCandRejectMultiPeak++;
          continue;
        }
      }
    }

    const double score = s_peak - pre_max;
    tmp.push_back({t0, score});
  }

  if (tmp.empty()) return;

  if (!kUsePrimaryPulseOnly) {
    out_t0.reserve(tmp.size());
    for (const auto& c : tmp) out_t0.push_back(c.t0);
    gCandAcceptedForFit += (long)out_t0.size();
    return;
  }

  const int nkeep = std::max(1, std::min((int)tmp.size(), kMaxCandidatesToFitPrimary));
  std::partial_sort(tmp.begin(), tmp.begin() + nkeep, tmp.end(),
                    [](const Cand& a, const Cand& b) {
                      if (a.score != b.score) return a.score > b.score;
                      return a.t0 < b.t0;
                    });

  out_t0.reserve(nkeep);
  for (int i = 0; i < nkeep; ++i) out_t0.push_back(tmp[i].t0);
  gCandAcceptedForFit += (long)out_t0.size();
}

// ---------- pre-trigger でのノイズ RMS 推定（chi2 重み用） ----------
static double EstimateNoiseSigma(const std::vector<double>& s,
                                 int t0_ref,
                                 int halfwin)
{
  if (s.empty()) return 1.0;

  const int n = (int)s.size();
  const int i0 = std::max(0, t0_ref - halfwin);
  const int i1 = std::min(n - 1, t0_ref - 10); // 立ち上がり付近を除外
  if (i1 < i0) return 1.0;

  const int m = i1 - i0 + 1;
  if (m < 20) return 1.0;

  double mean = 0.0;
  for (int i = i0; i <= i1; ++i) mean += s[i];
  mean /= (double)m;

  double v = 0.0;
  for (int i = i0; i <= i1; ++i) {
    const double d = s[i] - mean;
    v += d * d;
  }
  v /= (double)(m - 1);

  double sigma = std::sqrt(std::max(0.0, v));
  if (!std::isfinite(sigma) || sigma < 1.0) sigma = 1.0;
  return sigma;
}

// ---------- 候補直前の局所オフセットを頑健推定 ----------
static double EstimateLocalPreOffset(const std::vector<double>& s,
                                     int t0_ref)
{
  if (s.empty()) return 0.0;
  const int n = (int)s.size();

  const int i0 = std::max(0, t0_ref - 80);
  const int i1 = std::min(n - 1, t0_ref - 20);
  if (i1 < i0) return 0.0;

  std::vector<double> v;
  v.reserve(i1 - i0 + 1);
  for (int i = i0; i <= i1; ++i) v.push_back(s[i]);
  if (v.size() < 10) return 0.0;

  std::sort(v.begin(), v.end());
  const size_t lo = (size_t)std::floor(0.2 * (double)v.size());
  const size_t hi = (size_t)std::floor(0.8 * (double)v.size());
  if (hi <= lo) return 0.0;

  double sum = 0.0;
  int m = 0;
  for (size_t i = lo; i < hi; ++i) {
    sum += v[i];
    m++;
  }
  return (m > 0) ? (sum / (double)m) : 0.0;
}

// ---------- 平均波形の pre-trigger オフセットを 0 に合わせる ----------
static void SubtractPreTriggerOffset(std::vector<double>& avg,
                                     const std::vector<double>& xdt)
{
  if (avg.empty() || xdt.size() != avg.size()) return;

  double sum = 0.0;
  int n = 0;
  for (size_t i = 0; i < avg.size(); ++i) {
    if (xdt[i] <= -30.0) { // パルス立ち上がり前の領域
      sum += avg[i];
      n++;
    }
  }
  if (n <= 0) return;

  const double offset = sum / (double)n;
  for (double& y : avg) y -= offset;
}

// ---------- 分位点（0..1） ----------
static double QuantileOf(std::vector<double> v, double q)
{
  if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
  if (q < 0.0) q = 0.0;
  if (q > 1.0) q = 1.0;

  std::sort(v.begin(), v.end());
  const size_t idx = (size_t)std::llround(q * (double)(v.size() - 1));
  return v[idx];
}

// ============================================================
// テンプレ関数（最小モデル）
// s(x) = C + A * ( exp(-(x-t0)/(tr+d)) - exp(-(x-t0)/tr) ) * Theta(x-t0)
// パラメータ: [0]=A, [1]=t0, [2]=C, [3]=tr, [4]=d  (tf = tr + d)
// ============================================================

static double NaIPulseFunc_ExpDiff(double* x, double* p)
{
  const double A  = p[0];
  const double t0 = p[1];
  const double C  = p[2];
  const double tr = p[3];
  const double d  = p[4];
  const double tf = tr + d;

  if (tr <= 0.0 || d <= 0.0 || tf <= 0.0) return C;

  const double t = x[0] - t0;
  if (t < 0.0) return C;

  return C + A * (std::exp(-t/tf) - std::exp(-t/tr));
}

// g(t) = exp(-t/tf) - exp(-t/tr) の最大値（t>=0）をざっくり求める（A 初期値用）
static double Gmax_ExpDiff(double tr, double d, int tmax_samples)
{
  const double tf = tr + d;
  if (tr <= 0 || d <= 0 || tf <= 0) return 1.0;

  double gmax = 0.0;
  for (int t = 0; t <= tmax_samples; ++t) {
    const double g = std::exp(-(double)t/tf) - std::exp(-(double)t/tr);
    if (g > gmax) gmax = g;
  }
  if (gmax <= 0.0) gmax = 1.0;
  return gmax;
}

// C=0 のときの波形カーネル（A=1）
static double PulseKernelNoOffset(double x, double t0, double tr, double d)
{
  const double tf = tr + d;
  if (tr <= 0.0 || d <= 0.0 || tf <= 0.0) return 0.0;
  const double t = x - t0;
  if (t < 0.0) return 0.0;
  return std::exp(-t/tf) - std::exp(-t/tr);
}

// C=0 固定・tr,d 固定で、t0 をグリッド探索し A を解析的に求める
static bool FitFixedBaselineByGrid(const std::vector<double>& x,
                                   const std::vector<double>& y,
                                   const std::vector<double>& ey,
                                   double tr,
                                   double d,
                                   double t0_center,
                                   double t0_range,
                                   double& t0_best,
                                   double& A_best,
                                   double& chi2ndf_best)
{
  const int n = (int)x.size();
  if (n <= 4 || y.size() != x.size() || ey.size() != x.size()) return false;
  if (tr <= 0.0 || d <= 0.0) return false;

  auto EvalAtT0 = [&](double t0, double& A, double& chi2) -> bool {
    double s1 = 0.0; // Σ (y*g/σ^2)
    double s2 = 0.0; // Σ (g^2/σ^2)
    for (int i = 0; i < n; ++i) {
      const double sigma = (ey[i] > 0.0) ? ey[i] : 1.0;
      const double w = 1.0 / (sigma * sigma);
      const double g = PulseKernelNoOffset(x[i], t0, tr, d);
      s1 += y[i] * g * w;
      s2 += g * g * w;
    }
    if (s2 <= 0.0 || !std::isfinite(s2)) return false;

    A = s1 / s2;
    if (!std::isfinite(A)) return false;
    if (A < 0.0) A = 0.0;

    chi2 = 0.0;
    for (int i = 0; i < n; ++i) {
      const double sigma = (ey[i] > 0.0) ? ey[i] : 1.0;
      const double yhat = A * PulseKernelNoOffset(x[i], t0, tr, d);
      const double r = (y[i] - yhat) / sigma;
      chi2 += r * r;
    }
    if (!std::isfinite(chi2)) return false;
    return true;
  };

  // 粗密2段探索:
  // 1) 粗探索で大域的に良い領域を掴む
  // 2) その近傍を細探索して精度を戻す
  const double tmin = t0_center - t0_range;
  const double tmax = t0_center + t0_range;
  const int ngrid_coarse = 25;
  const int ngrid_fine = 19;

  bool found = false;
  double best_chi2 = 0.0;
  double best_t0 = t0_center;
  double best_A = 0.0;

  for (int ig = 0; ig < ngrid_coarse; ++ig) {
    const double t0 = tmin + (tmax - tmin) * (double)ig / (double)(ngrid_coarse - 1);
    double A = 0.0, chi2 = 0.0;
    if (!EvalAtT0(t0, A, chi2)) continue;

    if (!found || chi2 < best_chi2) {
      found = true;
      best_chi2 = chi2;
      best_t0 = t0;
      best_A = A;
    }
  }
  if (!found) return false;

  const double coarse_step = (tmax - tmin) / (double)(ngrid_coarse - 1);
  const double fmin = std::max(tmin, best_t0 - 2.0 * coarse_step);
  const double fmax = std::min(tmax, best_t0 + 2.0 * coarse_step);
  for (int ig = 0; ig < ngrid_fine; ++ig) {
    const double t0 = fmin + (fmax - fmin) * (double)ig / (double)(ngrid_fine - 1);
    double A = 0.0, chi2 = 0.0;
    if (!EvalAtT0(t0, A, chi2)) continue;
    if (chi2 < best_chi2) {
      best_chi2 = chi2;
      best_t0 = t0;
      best_A = A;
    }
  }

  const double ndf = (double)(n - 2); // A,t0
  chi2ndf_best = (ndf > 0.0) ? (best_chi2 / ndf) : best_chi2;
  t0_best = best_t0;
  A_best = best_A;
  return std::isfinite(chi2ndf_best);
}

// ============================================================
// τ 推定: clean パルスを相対時間で平均し、その平均波形をフィットして tr,d を得る
// ============================================================

struct TauEstResult {
  double tr = 0.0; // [samples]
  double d  = 0.0; // [samples]
  double t0 = 0.0; // [samples]（平均波形上の微小ずれ）
  double chi2ndf = 0.0;
  int    n_used  = 0;
  bool   ok = false;
};

static TauEstResult EstimateTauFromAverage(const std::vector<double>& avg,
                                           const std::vector<double>& xdt,
                                           int n_used)
{
  TauEstResult out;
  out.n_used = n_used;

  if (avg.empty() || xdt.size() != avg.size() || n_used < kNTauEstMin) {
    return out;
  }

  // 原因対策:
  // τ推定失敗の主因は、平均波形に pre-trigger オフセットが残ることと、
  // 初期値1点だけに依存した最小化の不安定さだった。
  // ここでは (1) pre-trigger 平均を 0 に揃える、(2) 多初期値で再試行する。
  std::vector<double> avg0 = avg;
  SubtractPreTriggerOffset(avg0, xdt);

  TGraph gr((int)avg0.size(), xdt.data(), avg0.data());

  const double xmin = xdt.front();
  const double xmax = xdt.back();
  const double t0_min = -40.0;
  const double t0_max =  40.0;

  double ymax = 0.0;
  for (double y : avg0) if (y > ymax) ymax = y;
  if (ymax <= 0.0) return out;

  TauEstResult best;
  bool has_best = false;
  int trial_id = 0;

  const double t0_seed[] = {-20.0, -10.0, 0.0, 10.0, 20.0};
  const double tr_seed[] = {2.0, 3.0, 5.0, 8.0, 12.0, 20.0};
  const double d_seed[]  = {10.0, 20.0, 30.0, 50.0, 80.0, 120.0};

  for (double t0_init : t0_seed) {
    for (double tr_init : tr_seed) {
      for (double d_init : d_seed) {
        TF1 f(Form("f_tau_%d", trial_id++), NaIPulseFunc_ExpDiff, xmin, xmax, 5);
        f.SetParNames("A","t0","C","tr","d");

        const double gmax = Gmax_ExpDiff(tr_init, d_init, kTauEstHalfWindow);
        const double A_guess = std::max(1.0, ymax / gmax);

        f.SetParameter(0, A_guess);
        f.SetParameter(1, t0_init);
        f.SetParameter(2, 0.0);
        f.SetParameter(3, tr_init);
        f.SetParameter(4, d_init);

        f.SetParLimits(0, 0.0, 1e9);
        f.SetParLimits(1, t0_min, t0_max);
        if (kFixFitBaselineToZero) {
          f.FixParameter(2, 0.0);
        } else {
          f.SetParLimits(2, -1e4, 1e4);
        }
        f.SetParLimits(3, kTauR_min_phys, 200.0);
        f.SetParLimits(4, kTauD_min_phys, 2000.0);

        const int fitStatus = gr.Fit(&f, "QNR");
        if (fitStatus != 0) continue;

        TauEstResult cand;
        cand.n_used = n_used;
        cand.tr = f.GetParameter(3);
        cand.d  = f.GetParameter(4);
        cand.t0 = f.GetParameter(1);
        const double ndf = f.GetNDF();
        cand.chi2ndf = (ndf > 0) ? (f.GetChisquare() / ndf) : 0.0;
        cand.ok = (cand.tr >= kTauR_min_phys &&
                   cand.d  >= kTauD_min_phys &&
                   std::isfinite(cand.chi2ndf));
        if (cand.ok) {
          // 下限境界に張り付いた解は不安定解として除外する
          if (cand.tr <= kTauR_min_phys + kTauBoundaryEps) cand.ok = false;
          if (cand.d  <= kTauD_min_phys + kTauBoundaryEps) cand.ok = false;
        }
        if (!cand.ok) continue;

        if (!has_best || cand.chi2ndf < best.chi2ndf) {
          best = cand;
          has_best = true;
        }
      }
    }
  }

  if (!has_best) return out;
  out = best;
  return out;
}

// ============================================================
// per-pulse fit（tr,d 固定で A,t0,C をフィット）→ I=A*d
// ============================================================

static FitPulse FitOnePulse(const std::vector<double>& s,
                            int t0_init,
                            double tr_fix,
                            double d_fix,
                            int fit_halfwin,
                            double t0_fit_range)
{
  gFitCalls++;
  FitPulse fp_best;
  fp_best.t0_init = t0_init;
  bool has_best = false;

  const int n = (int)s.size();
  if (n <= 0) return fp_best;

  // pre-trigger は短めにして、前段の残留 tail が chi2 を支配しないようにする
  const int i0 = t0_init - kFitPreWindow;
  const int i1 = t0_init + fit_halfwin;
  if (i0 < 0 || i1 >= n) return fp_best;

  const int npts = i1 - i0 + 1;
  std::vector<double> x(npts);
  for (int i = 0; i < npts; ++i) {
    const int idx = i0 + i;
    x[i] = (double)idx;   // samples index
  }

  const double pre_offset_raw = EstimateLocalPreOffset(s, t0_init);

  for (int io = 0; io < 2; ++io) {
    if (io == 1) {
      if (std::abs(pre_offset_raw) <= 1.0e-9) break;
      if (has_best && fp_best.chi2ndf <= kOffsetRetryChi2) break;
    }
    const double pre_offset = (io == 0) ? pre_offset_raw : 0.0;

    std::vector<double> y(npts);
    double s_peak = 0.0;
    for (int i = 0; i < npts; ++i) {
      const int idx = i0 + i;
      y[i] = s[idx] - pre_offset;
      if (y[i] > s_peak) s_peak = y[i];
    }
    if (s_peak <= 0.0) continue;

    const double sigma_noise = EstimateNoiseSigma(s, t0_init, fit_halfwin);
    std::vector<double> ex(npts, 0.0);
    std::vector<double> ey(npts, 0.0);
    for (int i = 0; i < npts; ++i) {
      const double ys = std::max(0.0, y[i]);
      const double sig_sys = kRelShapeSys * ys;
      ey[i] = std::sqrt(sigma_noise * sigma_noise + sig_sys * sig_sys);
      if (!std::isfinite(ey[i]) || ey[i] < 1.0) ey[i] = 1.0;
    }

    FitPulse fp;
    fp.t0_init = t0_init;
    fp.pre_offset = pre_offset;

    if (kFixFitBaselineToZero) {
      double t0_fit = 0.0;
      double A_fit = 0.0;
      double chi2ndf = 0.0;
      double used_range = t0_fit_range;
      const bool ok_grid = FitFixedBaselineByGrid(x, y, ey, tr_fix, d_fix,
                                                   (double)t0_init, used_range,
                                                   t0_fit, A_fit, chi2ndf);
      if (!ok_grid) continue;

      // 端張り付きなら探索幅を拡大して再探索
      if (kAdaptiveT0Range) {
        const double dt0 = std::abs(t0_fit - (double)t0_init);
        if (dt0 > kT0BoundaryFracExpand * used_range) {
          const double range2 = used_range * kT0FitRangeExpandFactor;
          double t0_fit2 = 0.0;
          double A_fit2 = 0.0;
          double chi2ndf2 = 0.0;
          const bool ok2 = FitFixedBaselineByGrid(x, y, ey, tr_fix, d_fix,
                                                  (double)t0_init, range2,
                                                  t0_fit2, A_fit2, chi2ndf2);
          if (ok2) {
            t0_fit = t0_fit2;
            A_fit = A_fit2;
            chi2ndf = chi2ndf2;
            used_range = range2;
          }
        }
      }

      fp.A = A_fit;
      fp.t0_fit = t0_fit;
      fp.C = 0.0;
      fp.chi2ndf = chi2ndf;
      fp.I = fp.A * d_fix;
      fp.ok = (fp.A >= 0.0 && std::isfinite(fp.chi2ndf));
      if (fp.ok) {
        const double dt0 = std::abs(fp.t0_fit - (double)t0_init);
        if (dt0 > kT0BoundaryFracReject * used_range) fp.ok = false;
        if (fp.t0_fit < kT0FitAbsMin) fp.ok = false;
        if (fp.t0_fit > (double)(n - 1)) fp.ok = false;
      }
      if (!fp.ok) continue;
    } else {
      TGraphErrors gr(npts, x.data(), y.data(), ex.data(), ey.data());
      TF1 f(Form("f_pulse_%d", io), NaIPulseFunc_ExpDiff, (double)i0, (double)i1, 5);
      f.SetParNames("A","t0","C","tr","d");
      f.FixParameter(3, tr_fix);
      f.FixParameter(4, d_fix);

      const double gmax = Gmax_ExpDiff(tr_fix, d_fix, fit_halfwin);
      const double A_guess = std::max(0.0, s_peak / gmax);

      f.SetParameter(0, A_guess);
      f.SetParameter(1, (double)t0_init);
      f.SetParameter(2, 0.0);

      f.SetParLimits(0, 0.0, 1e12);
      f.SetParLimits(1, (double)t0_init - t0_fit_range, (double)t0_init + t0_fit_range);
      f.SetParLimits(2, -1e5, 1e5);

      const int fitStatus = gr.Fit(&f, "QNR");
      if (fitStatus != 0) continue;

      fp.A = f.GetParameter(0);
      fp.t0_fit = f.GetParameter(1);
      fp.C = f.GetParameter(2);
      const double ndf = f.GetNDF();
      fp.chi2ndf = (ndf > 0) ? (f.GetChisquare() / ndf) : 0.0;
      fp.I = fp.A * d_fix;
      fp.ok = (fp.A >= 0.0 && std::isfinite(fp.chi2ndf));
      if (fp.ok) {
        const double dt0 = std::abs(fp.t0_fit - (double)t0_init);
        if (dt0 > kT0BoundaryFracReject * t0_fit_range) fp.ok = false;
        if (fp.t0_fit < kT0FitAbsMin) fp.ok = false;
        if (fp.t0_fit > (double)(n - 1)) fp.ok = false;
      }
      if (!fp.ok) continue;
    }

    if (!has_best || fp.chi2ndf < fp_best.chi2ndf) {
      fp_best = fp;
      has_best = true;
    }
  }

  if (!has_best) return FitPulse{};
  return fp_best;
}

static bool HasNearbyT0(const std::vector<int>& t0s, int t0, int sep)
{
  for (int x : t0s) {
    if (std::abs(x - t0) <= sep) return true;
  }
  return false;
}

static double ResidualGlobalPeak(const std::vector<double>& s)
{
  double mx = 0.0;
  for (double y : s) if (y > mx) mx = y;
  return mx;
}

// 受理した1パルスを s 波形から引いて、重なりの次パルスを見えやすくする
static void SubtractAcceptedPulse(std::vector<double>& s_work,
                                  const FitPulse& fp,
                                  double tr_fix,
                                  double d_fix)
{
  if (!fp.ok) return;
  if (fp.A <= 0.0 || tr_fix <= 0.0 || d_fix <= 0.0) return;

  const int n = (int)s_work.size();
  if (n <= 0) return;

  const double tf = tr_fix + d_fix;
  if (tf <= 0.0) return;

  const int i0 = std::max(0, (int)std::floor(fp.t0_fit) - 2);
  const int i1 = std::min(n - 1, (int)std::ceil(fp.t0_fit + kSubtractTailFactor * tf));
  for (int i = i0; i <= i1; ++i) {
    const double yhat = fp.A * PulseKernelNoOffset((double)i, fp.t0_fit, tr_fix, d_fix);
    s_work[i] -= yhat;
  }
}

// 分離で受理する成分の「有意性」を評価する
// - その成分で局所 χ² が十分改善するか
// - モデルピークがノイズに対して十分有意か
static bool PassDecompositionQuality(const std::vector<double>& s_work,
                                     const FitPulse& fp,
                                     double tr_fix,
                                     double d_fix)
{
  if (!fp.ok) return false;
  if (fp.A <= 0.0 || tr_fix <= 0.0 || d_fix <= 0.0) return false;

  const int n = (int)s_work.size();
  const int t0i = (int)std::lround(fp.t0_fit);
  const int i0 = std::max(0, t0i - kFitPreWindow);
  const int i1 = std::min(n - 1, t0i + kFitHalfWindow);
  if (i1 - i0 + 1 < 20) return false;

  const double sigma_noise = EstimateNoiseSigma(s_work, t0i, kFitHalfWindow);
  if (!std::isfinite(sigma_noise) || sigma_noise <= 0.0) return false;

  double chi2_null = 0.0;
  double chi2_model = 0.0;
  for (int i = i0; i <= i1; ++i) {
    const double y = s_work[i] - fp.pre_offset;
    const double yhat = fp.A * PulseKernelNoOffset((double)i, fp.t0_fit, tr_fix, d_fix);
    const double sig_sys = kRelShapeSys * std::max(0.0, y);
    double sigma = std::sqrt(sigma_noise * sigma_noise + sig_sys * sig_sys);
    if (!std::isfinite(sigma) || sigma < 1.0) sigma = 1.0;
    const double r0 = y / sigma;
    const double r1 = (y - yhat) / sigma;
    chi2_null += r0 * r0;
    chi2_model += r1 * r1;
  }
  if (!std::isfinite(chi2_null) || !std::isfinite(chi2_model)) return false;
  const double delta_chi2 = chi2_null - chi2_model;
  const double rel_improve = delta_chi2 / std::max(1.0, chi2_null);

  const double model_peak = fp.A * Gmax_ExpDiff(tr_fix, d_fix, kFitHalfWindow);
  const double peak_req = std::max(kDecompMinModelPeakAbs, kDecompMinModelPeakSigma * sigma_noise);
  if (model_peak < peak_req) return false;
  if (delta_chi2 < kDecompMinDeltaChi2) return false;
  if (rel_improve < kDecompMinRelImprove) return false;
  return true;
}

// 分離成分が「元の生波形」に実在するかを確認する
// 目的:
//   残差再探索で生じる偽ピーク（指数尾の近似誤差）を受理しない。
static bool PassDataPulseEvidence(const std::vector<double>& s_data,
                                  const FitPulse& fp)
{
  if (!fp.ok) return false;
  const int n = (int)s_data.size();
  if (n < 5) return false;

  const int t0i = (int)std::lround(fp.t0_fit);
  if (t0i < 1 || t0i >= n - 1) return false;

  const int i0 = std::max(1, t0i - kDataEvPre);
  const int i1 = std::min(n - 2, t0i + kDataEvPost);
  if (i1 <= i0) return false;

  int ip = -1;
  double sp = -1.0e9;
  for (int i = i0; i <= i1; ++i) {
    if (s_data[i] > sp) {
      sp = s_data[i];
      ip = i;
    }
  }
  if (ip < 1 || ip >= n - 1) return false;
  if (sp < kDataEvPeakMin) return false;

  double pre_min = sp;
  const int j0 = std::max(0, ip - kDataEvPre);
  const int j1 = std::max(j0, ip - 2);
  for (int j = j0; j <= j1; ++j) {
    if (s_data[j] < pre_min) pre_min = s_data[j];
  }
  const double prom = sp - pre_min;
  if (prom < kDataEvPromMin) return false;

  double rise_max = 0.0;
  for (int j = std::max(1, ip - 8); j <= ip; ++j) {
    const double ds = s_data[j] - s_data[j - 1];
    if (ds > rise_max) rise_max = ds;
  }
  if (rise_max < kDataEvRiseMin) return false;

  return true;
}

static bool SolveLinearSystem(std::vector<std::vector<double>> a,
                              std::vector<double> b,
                              std::vector<double>& x)
{
  const int n = (int)a.size();
  if (n <= 0 || (int)b.size() != n) return false;
  for (int i = 0; i < n; ++i) {
    if ((int)a[i].size() != n) return false;
  }

  for (int i = 0; i < n; ++i) {
    int piv = i;
    double vmax = std::abs(a[i][i]);
    for (int r = i + 1; r < n; ++r) {
      const double v = std::abs(a[r][i]);
      if (v > vmax) { vmax = v; piv = r; }
    }
    if (vmax <= 1.0e-12 || !std::isfinite(vmax)) return false;
    if (piv != i) {
      std::swap(a[piv], a[i]);
      std::swap(b[piv], b[i]);
    }

    const double diag = a[i][i];
    for (int c = i; c < n; ++c) a[i][c] /= diag;
    b[i] /= diag;

    for (int r = 0; r < n; ++r) {
      if (r == i) continue;
      const double f = a[r][i];
      if (f == 0.0) continue;
      for (int c = i; c < n; ++c) a[r][c] -= f * a[i][c];
      b[r] -= f * b[i];
    }
  }

  x = b;
  for (double xi : x) {
    if (!std::isfinite(xi)) return false;
  }
  return true;
}

static double ComputeChi2ForSet(const std::vector<double>& s_data,
                                const std::vector<FitPulse>& pulses,
                                double tr_fix,
                                double d_fix,
                                double sigma_noise,
                                int skip_idx)
{
  if (sigma_noise < 1.0) sigma_noise = 1.0;
  const double inv = 1.0 / (sigma_noise * sigma_noise);
  double chi2 = 0.0;
  const int n = (int)s_data.size();
  for (int i = 0; i < n; ++i) {
    double yhat = 0.0;
    for (int j = 0; j < (int)pulses.size(); ++j) {
      if (j == skip_idx) continue;
      if (!pulses[j].ok || pulses[j].A <= 0.0) continue;
      yhat += pulses[j].A * PulseKernelNoOffset((double)i, pulses[j].t0_fit, tr_fix, d_fix);
    }
    const double r = s_data[i] - yhat;
    chi2 += r * r * inv;
  }
  return chi2;
}

static bool SolveAmplitudesAndChi2(const std::vector<double>& s_data,
                                   const std::vector<FitPulse>& pulses,
                                   double tr_fix,
                                   double d_fix,
                                   double sigma_noise,
                                   std::vector<double>& A_out,
                                   double& chi2_out)
{
  const int m = (int)pulses.size();
  const int n = (int)s_data.size();
  if (m <= 0 || n <= 0) return false;
  if (sigma_noise < 1.0) sigma_noise = 1.0;
  const double inv = 1.0 / (sigma_noise * sigma_noise);

  std::vector<std::vector<double>> M(m, std::vector<double>(m, 0.0));
  std::vector<double> b(m, 0.0);

  for (int i = 0; i < n; ++i) {
    std::vector<double> g(m, 0.0);
    for (int j = 0; j < m; ++j) {
      g[j] = PulseKernelNoOffset((double)i, pulses[j].t0_fit, tr_fix, d_fix);
    }
    for (int j = 0; j < m; ++j) {
      b[j] += s_data[i] * g[j] * inv;
      for (int k = j; k < m; ++k) {
        M[j][k] += g[j] * g[k] * inv;
      }
    }
  }
  for (int j = 0; j < m; ++j) {
    for (int k = 0; k < j; ++k) M[j][k] = M[k][j];
  }

  std::vector<double> A;
  if (!SolveLinearSystem(M, b, A)) return false;
  for (double a : A) {
    if (!std::isfinite(a)) return false;
  }

  chi2_out = 0.0;
  for (int i = 0; i < n; ++i) {
    double yhat = 0.0;
    for (int j = 0; j < m; ++j) {
      yhat += A[j] * PulseKernelNoOffset((double)i, pulses[j].t0_fit, tr_fix, d_fix);
    }
    const double r = s_data[i] - yhat;
    chi2_out += r * r * inv;
  }
  if (!std::isfinite(chi2_out)) return false;
  A_out.swap(A);
  return true;
}

static void ExtractEvidencePeaks(const std::vector<double>& s_data,
                                 std::vector<int>& peaks_out)
{
  peaks_out.clear();
  const int n = (int)s_data.size();
  if (n < 5) return;

  int last = -kDecompMinT0Sep;
  for (int i = 2; i <= n - 3; ++i) {
    if (i - last < kDecompMinT0Sep) continue;
    if (!(s_data[i] >= s_data[i - 1] && s_data[i] > s_data[i + 1])) continue;
    if (s_data[i] < kDataEvPeakMin) continue;

    double pre_min = s_data[i];
    const int j0 = std::max(0, i - kDataEvPre);
    const int j1 = std::max(j0, i - 2);
    for (int j = j0; j <= j1; ++j) {
      if (s_data[j] < pre_min) pre_min = s_data[j];
    }
    const double prom = s_data[i] - pre_min;
    if (prom < kDataEvPromMin) continue;

    peaks_out.push_back(i);
    last = i;
  }
}

static bool PruneByUniquePeakAssignment(const std::vector<double>& s_data,
                                        std::vector<FitPulse>& pulses)
{
  if (pulses.size() <= 1) return false;

  std::vector<int> peaks;
  ExtractEvidencePeaks(s_data, peaks);
  if (peaks.empty()) return false;

  std::vector<int> p2k(pulses.size(), -1);
  std::vector<int> p2d(pulses.size(), 1e9);
  for (size_t ip = 0; ip < pulses.size(); ++ip) {
    const int t0 = (int)std::lround(pulses[ip].t0_fit);
    for (size_t ik = 0; ik < peaks.size(); ++ik) {
      const int d = std::abs(t0 - peaks[ik]);
      if (d < p2d[ip]) {
        p2d[ip] = d;
        p2k[ip] = (int)ik;
      }
    }
  }

  std::vector<bool> keep(pulses.size(), true);
  for (size_t ik = 0; ik < peaks.size(); ++ik) {
    int best = -1;
    double best_score = -1.0;
    for (size_t ip = 0; ip < pulses.size(); ++ip) {
      if (p2k[ip] != (int)ik) continue;
      if (p2d[ip] > kPeakAssignMaxDist) continue;
      const double score = pulses[ip].I;
      if (score > best_score) {
        best_score = score;
        best = (int)ip;
      }
    }
    for (size_t ip = 0; ip < pulses.size(); ++ip) {
      if (p2k[ip] == (int)ik && p2d[ip] <= kPeakAssignMaxDist) {
        if ((int)ip != best) keep[ip] = false;
      }
    }
  }

  bool removed = false;
  for (int ip = (int)pulses.size() - 1; ip >= 0; --ip) {
    if (!keep[ip]) {
      pulses.erase(pulses.begin() + ip);
      removed = true;
    }
  }
  return removed;
}

static void RefineAcceptedPulsesJoint(const std::vector<double>& s_data,
                                      double tr_fix,
                                      double d_fix,
                                      std::vector<FitPulse>& pulses)
{
  if (pulses.size() <= 1) return;
  if (tr_fix <= 0.0 || d_fix <= 0.0) return;
  const int n = (int)s_data.size();
  if (n < 20) return;

  std::sort(pulses.begin(), pulses.end(),
            [](const FitPulse& a, const FitPulse& b) { return a.t0_fit < b.t0_fit; });

  const int guard_loops = 10;
  for (int loop = 0; loop < guard_loops; ++loop) {
    const int m = (int)pulses.size();
    if (m <= 1) break;

    const double sigma_noise = EstimateNoiseSigma(s_data, (int)std::lround(pulses[0].t0_fit), kFitHalfWindow);

    // 1) t0 固定で A を解く
    std::vector<double> A(m, 0.0);
    double chi2_now = 0.0;
    if (!SolveAmplitudesAndChi2(s_data, pulses, tr_fix, d_fix, sigma_noise, A, chi2_now)) break;
    bool removed_nonpos = false;
    for (int j = m - 1; j >= 0; --j) {
      if (!(A[j] > 0.0)) {
        pulses.erase(pulses.begin() + j);
        removed_nonpos = true;
      }
    }
    if (removed_nonpos) continue;
    for (int j = 0; j < m; ++j) {
      pulses[j].A = A[j];
      pulses[j].I = pulses[j].A * d_fix;
      pulses[j].ok = (pulses[j].A > 0.0);
    }

    // 2) 交互最適化で t0 も調整（pileupピーク重なり対策）
    bool restart_outer = false;
    for (int it = 0; it < kJointRefineIters; ++it) {
      bool moved = false;
      for (int j = 0; j < m; ++j) {
        const double t0_org = pulses[j].t0_fit;
        double best_t0 = t0_org;
        double best_chi2 = chi2_now;

        for (int dt = -kJointT0ScanRange; dt <= kJointT0ScanRange; dt += kJointT0ScanStep) {
          const double t0_try = t0_org + (double)dt;
          if (t0_try < kT0FitAbsMin || t0_try > (double)(n - 1)) continue;

          bool too_close = false;
          for (int k = 0; k < m; ++k) {
            if (k == j) continue;
            if (std::abs(t0_try - pulses[k].t0_fit) < kJointT0MinSepHard) {
              too_close = true;
              break;
            }
          }
          if (too_close) continue;

          pulses[j].t0_fit = t0_try;
          std::vector<double> A_try;
          double chi2_try = 0.0;
          if (!SolveAmplitudesAndChi2(s_data, pulses, tr_fix, d_fix, sigma_noise, A_try, chi2_try)) {
            continue;
          }
          bool all_pos = true;
          for (double a : A_try) {
            if (!(a > 0.0)) { all_pos = false; break; }
          }
          if (!all_pos) continue;
          if (chi2_try < best_chi2) {
            best_chi2 = chi2_try;
            best_t0 = t0_try;
          }
        }

        pulses[j].t0_fit = best_t0;
        if (std::abs(best_t0 - t0_org) > 1.0e-9) moved = true;
      }

      std::vector<double> A_ref;
      double chi2_ref = 0.0;
      if (!SolveAmplitudesAndChi2(s_data, pulses, tr_fix, d_fix, sigma_noise, A_ref, chi2_ref)) break;
      bool ref_nonpos = false;
      for (int j = (int)A_ref.size() - 1; j >= 0; --j) {
        if (!(A_ref[j] > 0.0)) {
          pulses.erase(pulses.begin() + j);
          ref_nonpos = true;
        }
      }
      if (ref_nonpos) {
        restart_outer = true;
        break;
      }
      for (int j = 0; j < m; ++j) {
        pulses[j].A = A_ref[j];
        pulses[j].I = pulses[j].A * d_fix;
        pulses[j].ok = (pulses[j].A > 0.0);
      }
      chi2_now = chi2_ref;
      if (!moved) break;
    }
    if (restart_outer) continue;

    // 3) 同一実ピークへ複数成分が張り付く場合は1本化
    if (PruneByUniquePeakAssignment(s_data, pulses)) continue;

    const double peak_req = std::max(kDecompMinModelPeakAbs, kDecompMinModelPeakSigma * std::max(1.0, sigma_noise));
    const double chi2_all = ComputeChi2ForSet(s_data, pulses, tr_fix, d_fix, sigma_noise, -1);

    bool removed = false;
    for (int j = m - 1; j >= 0; --j) {
      const double model_peak = pulses[j].A * Gmax_ExpDiff(tr_fix, d_fix, kFitHalfWindow);
      if (model_peak < peak_req) {
        pulses.erase(pulses.begin() + j);
        removed = true;
        continue;
      }
      if (!PassDataPulseEvidence(s_data, pulses[j])) {
        pulses.erase(pulses.begin() + j);
        removed = true;
        continue;
      }
      const double chi2_wo = ComputeChi2ForSet(s_data, pulses, tr_fix, d_fix, sigma_noise, j);
      const double delta = chi2_wo - chi2_all;
      if (delta < kDecompMinDeltaChi2) {
        pulses.erase(pulses.begin() + j);
        removed = true;
        continue;
      }
    }
    if (!removed) break;
  }

  // t0 が近すぎる成分を統合（強い方を残す）
  if (pulses.size() >= 2) {
    std::sort(pulses.begin(), pulses.end(),
              [](const FitPulse& a, const FitPulse& b) { return a.t0_fit < b.t0_fit; });
    std::vector<FitPulse> keep;
    keep.reserve(pulses.size());
    keep.push_back(pulses[0]);
    for (size_t i = 1; i < pulses.size(); ++i) {
      FitPulse cur = pulses[i];
      FitPulse& prev = keep.back();
      if (std::abs(cur.t0_fit - prev.t0_fit) < kDecompMinT0Sep) {
        if (cur.I > prev.I) prev = cur;
      } else {
        keep.push_back(cur);
      }
    }
    pulses.swap(keep);
  }
}

// 多パルス分離:
//   候補を順次 fit → 受理パルスを引く → 残差から候補再探索
static void FitPulsesWithDecomposition(const std::vector<double>& s_input,
                                       const std::vector<int>& starts_seed,
                                       double tr_fix,
                                       double d_fix,
                                       double chi2_cut,
                                       bool apply_chi2_cut,
                                       std::vector<FitPulse>& accepted_out)
{
  accepted_out.clear();
  if (s_input.empty()) return;

  std::vector<double> s_work = s_input;
  std::vector<int> tried_t0;
  std::vector<int> accepted_t0;
  std::vector<int> pending;

  auto EnqueueFromWave = [&](const std::vector<double>& swave) {
    std::vector<int> starts = FindStartCandidates(swave, kSeedThr, kEndThr, kMinSepSamples);
    AddPeakSupplementCandidates(swave, starts);
    AddGlobalPeakCandidate(swave, starts);

    std::vector<int> fit_cands;
    BuildFitCandidates(starts, swave, kFitPreWindow, kFitHalfWindow, fit_cands);
    for (int t0 : fit_cands) {
      if (HasNearbyT0(tried_t0, t0, kCandMergeMinSep)) continue;
      if (HasNearbyT0(accepted_t0, t0, kCandMergeMinSep)) continue;
      if (HasNearbyT0(pending, t0, kCandMergeMinSep)) continue;
      pending.push_back(t0);
    }
  };

  // 初期候補（元波形）
  {
    std::vector<int> fit_cands;
    BuildFitCandidates(starts_seed, s_work, kFitPreWindow, kFitHalfWindow, fit_cands);
    for (int t0 : fit_cands) pending.push_back(t0);
  }

  for (int iter = 0; iter < kMaxTriedPerEventCh; ++iter) {
    if (pending.empty()) break;
    if ((int)accepted_out.size() >= kMaxAcceptedPerEventCh) break;
    if (kEnablePulseDecomposition) {
      if (ResidualGlobalPeak(s_work) < kDecompResidualStopThr) break;
    }

    // 現在の残差で最も強い候補を優先する
    int ibest = 0;
    double best_score = -1.0;
    for (size_t i = 0; i < pending.size(); ++i) {
      const double sc = CandidatePeakScore(s_work, pending[i]);
      if (sc > best_score) {
        best_score = sc;
        ibest = (int)i;
      }
    }
    const int t0 = pending[ibest];
    pending.erase(pending.begin() + ibest);

    if (HasNearbyT0(tried_t0, t0, kCandMergeMinSep)) continue;
    tried_t0.push_back(t0);

    FitPulse fp = FitOnePulse(s_work, t0, tr_fix, d_fix, kFitHalfWindow, kT0FitRange);
    if (!fp.ok) continue;
    if (apply_chi2_cut && fp.chi2ndf > chi2_cut) continue;
    if (!PassPrimaryPeakConsistency(s_input, fp.t0_fit)) continue;
    if (HasNearbyT0(accepted_t0, (int)std::lround(fp.t0_fit), kDecompMinT0Sep)) continue;
    if (!PassDecompositionQuality(s_work, fp, tr_fix, d_fix)) continue;
    if (!PassDataPulseEvidence(s_input, fp)) continue;

    accepted_out.push_back(fp);
    accepted_t0.push_back((int)std::lround(fp.t0_fit));

    if (kEnablePulseDecomposition) {
      SubtractAcceptedPulse(s_work, fp, tr_fix, d_fix);
      pending.clear();
      EnqueueFromWave(s_work);
    }
  }

  if (accepted_out.size() >= 2) {
    RefineAcceptedPulsesJoint(s_input, tr_fix, d_fix, accepted_out);
  }
}

// ============================================================
// 描画ユーティリティ
// ============================================================

static void SetupPad(TPad* p)
{
  p->SetFillStyle(0);
  p->SetLeftMargin(0.12);
  p->SetRightMargin(0.05);
  p->SetTopMargin(0.08);
  p->SetBottomMargin(0.12);
}

static void DrawLabelTopLeft(const char* label)
{
  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextSize(0.06);
  lat.DrawLatex(0.14, 0.88, label);
}

static void DrawDebugPage(const DebugEvent& ev,
                          const double tr_fix[4],
                          const double d_fix[4],
                          const TString& out_pdf,
                          bool is_last_page)
{
  TCanvas* c = new TCanvas(Form("c_dbg_nai_%ld", ev.idx),
                           Form("NaI debug event %ld", ev.idx),
                           kCanvasW, kCanvasH);

  // 2x2 pads
  TPad* p1 = new TPad(Form("p_dbg1_%ld", ev.idx), "", 0.00, 0.50, 0.50, 1.00);
  TPad* p2 = new TPad(Form("p_dbg2_%ld", ev.idx), "", 0.50, 0.50, 1.00, 1.00);
  TPad* p3 = new TPad(Form("p_dbg3_%ld", ev.idx), "", 0.00, 0.00, 0.50, 0.50);
  TPad* p4 = new TPad(Form("p_dbg4_%ld", ev.idx), "", 0.50, 0.00, 1.00, 0.50);
  TPad* pads[4] = {p1,p2,p3,p4};

  for (int i = 0; i < 4; ++i) { SetupPad(pads[i]); pads[i]->Draw(); }

  // 各chで s(t)=baseline-ADC を描く＋フィット曲線
  for (int ch = 0; ch < 4; ++ch) {
    pads[ch]->cd();

    const int n = (int)ev.adc[ch].size();
    std::vector<double> tx(n), sy(n);
    for (int i = 0; i < n; ++i) {
      tx[i] = i * kDtNs;                 // ns
      sy[i] = ev.baseline[ch] - ev.adc[ch][i]; // s
    }

    TGraph* gr = new TGraph(n, tx.data(), sy.data());
    gr->SetTitle(Form("%s  (digitizer event %ld);time [ns];s = baseline - ADC [counts]",
                      kChLabel[ch], ev.idx));
    gr->Draw("AL");

    // 閾値線（s空間）
    TLine* l0 = new TLine(tx.front(), 0.0, tx.back(), 0.0);
    l0->SetLineStyle(2); l0->SetLineWidth(2); l0->Draw();

    TLine* lseed = new TLine(tx.front(), kSeedThr, tx.back(), kSeedThr);
    lseed->SetLineStyle(3); lseed->SetLineWidth(2); lseed->Draw();

    TLine* lend = new TLine(tx.front(), kEndThr, tx.back(), kEndThr);
    lend->SetLineStyle(3); lend->SetLineWidth(1); lend->Draw();

    // fit curve（ヒストへ入れた全パルスを描画）
    const std::vector<FitPulse>& fps = ev.fit_hist[ch];
    if (!fps.empty()) {
      // x_ns を samples に戻して関数評価する
      struct Local {
        static double Eval(double* x, double* p) {
          const double x_ns = x[0];
          const double dt_samples = x_ns / kDtNs;
          double pp[5];
          pp[0]=p[0]; pp[1]=p[1]; pp[2]=p[2]; pp[3]=p[3]; pp[4]=p[4];
          double xs = dt_samples;
          return NaIPulseFunc_ExpDiff(&xs, pp);
        }
      };

      for (size_t jf = 0; jf < fps.size(); ++jf) {
        const FitPulse& fpj = fps[jf];
        if (!fpj.ok) continue;

        TF1 fdraw(Form("fdraw_%ld_%d_%zu", ev.idx, ch, jf), Local::Eval, 0.0, (n-1)*kDtNs, 5);
        fdraw.SetParameters(fpj.A, fpj.t0_fit, fpj.C, tr_fix[ch], d_fix[ch]);
        fdraw.SetLineWidth((jf == 0) ? 2 : 1);
        fdraw.SetLineColor(kRed + 1);
        fdraw.DrawCopy("SAME");

        TLine* lt0 = new TLine(fpj.t0_fit * kDtNs, 0.0, fpj.t0_fit * kDtNs, std::max(1.0, gr->GetYaxis()->GetXmax()));
        lt0->SetLineStyle(3);
        lt0->SetLineWidth(1);
        lt0->SetLineColor(kRed + 1);
        lt0->Draw();
      }

      // 全成分の和（見やすい太線）
      std::vector<double> ysum(n, 0.0);
      for (int i = 0; i < n; ++i) {
        double y = 0.0;
        for (const auto& fpj : fps) {
          if (!fpj.ok) continue;
          y += fpj.A * PulseKernelNoOffset((double)i, fpj.t0_fit, tr_fix[ch], d_fix[ch]);
        }
        ysum[i] = y;
      }
      TGraph* gsum = new TGraph(n, tx.data(), ysum.data());
      gsum->SetLineColor(kRed + 2);
      gsum->SetLineWidth(3);
      gsum->Draw("L SAME");
    }

    // 注釈
    TLatex t;
    t.SetNDC(true);
    t.SetTextSize(0.040);
    t.DrawLatex(0.14, 0.84, Form("%s", kChLabel[ch]));
    t.SetTextSize(0.035);
    t.DrawLatex(0.14, 0.78, Form("baseline(mode) = %.2f", ev.baseline[ch]));
    t.DrawLatex(0.14, 0.72, Form("cand=%d  fit_ok=%d  hist=%d", ev.n_cand[ch], ev.n_acc[ch], ev.n_hist[ch]));
    t.DrawLatex(0.14, 0.66, Form("tr=%.2f samp (%.1f ns)", tr_fix[ch], tr_fix[ch]*kDtNs));
    t.DrawLatex(0.14, 0.60, Form("d =%.2f samp (%.1f ns)", d_fix[ch],  d_fix[ch]*kDtNs));
    const FitPulse& fp = ev.fit[ch];
    t.DrawLatex(0.14, 0.54, Form("pre-offset(rep) = %.2f", fp.pre_offset));

    if (ev.n_hist[ch] > 0 && fp.ok) {
      t.DrawLatex(0.14, 0.48, Form("rep t0=%.2f samp (%.1f ns)", fp.t0_fit, fp.t0_fit*kDtNs));
      t.DrawLatex(0.14, 0.42, Form("rep I=A*d = %.1f (ADC*samp)  = %.1f (ADC*ns)", fp.I, fp.I*kDtNs));
      t.DrawLatex(0.14, 0.36, Form("rep chi2/ndf = %.2f", fp.chi2ndf));
    } else {
      t.DrawLatex(0.14, 0.48, "fit: (no accepted pulse for histogram)");
    }
  }

  // PDFへ
  if (is_last_page) {
    TString out_close = out_pdf; out_close += ")";
    c->SaveAs(out_close.Data());
  } else {
    c->SaveAs(out_pdf.Data());
  }
}

// ============================================================

void hist_NaI_fitint_clean()
{
  gSystem->mkdir(kOutputDir, true);
  TString out_pdf = Form("%s/hist_NaI_fitint_clean_%s.pdf", kOutputDir, kInputTag);

  // 診断カウンタ初期化
  gCandTotalInput = 0;
  gCandRejectEdge = 0;
  gCandRejectClean = 0;
  gCandRejectPreQuiet = 0;
  gCandRejectMultiPeak = 0;
  gCandAcceptedForFit = 0;
  gFitCalls = 0;

  // 入力ストリーム（pass1, pass2 で開き直す）
  TString fpath[4];
  for (int ch = 0; ch < 4; ++ch) fpath[ch] = Form("%s/%s", kInputDir, kNaIFiles[ch]);

  // ============================================================
  // Pass 1: clean サンプルで τ 推定（平均波形→フィット）
  // ============================================================

  // 平均波形用のアキュムレータ
  const int n_dt = 2*kTauEstHalfWindow + 1;
  std::vector<double> xdt(n_dt);
  for (int i = 0; i < n_dt; ++i) xdt[i] = (double)(i - kTauEstHalfWindow); // samples

  std::vector<double> sum[4];
  int n_used[4] = {0,0,0,0};
  for (int ch = 0; ch < 4; ++ch) sum[ch].assign(n_dt, 0.0);

  {
    std::ifstream in[4];
    for (int ch = 0; ch < 4; ++ch) in[ch].open(fpath[ch].Data());
    bool ok = true;
    for (int ch = 0; ch < 4; ++ch) if (!in[ch].is_open()) ok = false;
    if (!ok) {
      std::cerr << "ERROR: could not open NaI files (pass1)\n";
      for (int ch = 0; ch < 4; ++ch) std::cerr << "  " << fpath[ch] << "\n";
      return;
    }

    std::vector<double> buf_adc[4];
    for (int ch = 0; ch < 4; ++ch) { buf_adc[ch].reserve(kSamplesPerEvent); }

    while (true) {
      double v[4];
      bool read_ok = true;
      for (int ch = 0; ch < 4; ++ch) {
        if (!(in[ch] >> v[ch])) read_ok = false;
      }
      if (!read_ok) break;

      for (int ch = 0; ch < 4; ++ch) buf_adc[ch].push_back(v[ch]);

      bool full = true;
      for (int ch = 0; ch < 4; ++ch) if ((int)buf_adc[ch].size() != kSamplesPerEvent) full = false;
      if (!full) continue;

      // イベント処理：各chで clean 候補を探し、平均波形に足す
      for (int ch = 0; ch < 4; ++ch) {
        if (n_used[ch] >= kNTauEstTarget) continue;

        const double baseline = ModeBaselineInEvent(buf_adc[ch]);

        std::vector<double> s;
        BuildS(buf_adc[ch], baseline, s);

        std::vector<int> starts = FindStartCandidates(s, kSeedThr, kEndThr, kMinSepSamples);
        AddPeakSupplementCandidates(s, starts);
        AddGlobalPeakCandidate(s, starts);

        // τ推定では「イベント内で唯一」ではなく「局所窓で単発(clean)」を使う。
        // A2 のようにイベント内に複数候補がある channel でも、孤立パルスを活用できる。
        for (int t0 : starts) {
          if (n_used[ch] >= kNTauEstTarget) break;

          // fit窓が取れない端は捨てる
          if (t0 - kTauEstHalfWindow < 0) continue;
          if (t0 + kTauEstHalfWindow >= (int)s.size()) continue;
          if (!IsSingleCandidateInWindow(starts, t0, kTauEstHalfWindow)) continue;

          // 平均へ加算（相対時間で）
          for (int i = 0; i < n_dt; ++i) {
            const int dt = i - kTauEstHalfWindow;
            sum[ch][i] += s[t0 + dt];
          }
          n_used[ch]++;
        }
      }

      // 次イベント
      for (int ch = 0; ch < 4; ++ch) buf_adc[ch].clear();
      bool all_done = true;
      for (int ch = 0; ch < 4; ++ch) if (n_used[ch] < kNTauEstTarget) all_done = false;
      if (all_done) break;
    }

    for (int ch = 0; ch < 4; ++ch) in[ch].close();
  }

  // 平均波形を作って τ 推定
  TauEstResult tauRes[4];
  double tr_fix[4] = {kTauR_init, kTauR_init, kTauR_init, kTauR_init};
  double d_fix [4] = {kTauD_init, kTauD_init, kTauD_init, kTauD_init};

  std::vector<double> avg[4];
  for (int ch = 0; ch < 4; ++ch) {
    avg[ch].assign(n_dt, 0.0);
    if (n_used[ch] >= kNTauEstMin) {
      for (int i = 0; i < n_dt; ++i) avg[ch][i] = sum[ch][i] / (double)n_used[ch];
      SubtractPreTriggerOffset(avg[ch], xdt);
      tauRes[ch] = EstimateTauFromAverage(avg[ch], xdt, n_used[ch]);
      if (tauRes[ch].ok) {
        tr_fix[ch] = tauRes[ch].tr;
        d_fix[ch]  = tauRes[ch].d;
      }
    }
  }

  // 一部チャネルで τ 推定が失敗した場合、
  // 成功チャネルの平均 τ を fallback として使って形状崩壊を防ぐ
  double tr_mean_ok = 0.0;
  double d_mean_ok = 0.0;
  int n_ok = 0;
  for (int ch = 0; ch < 4; ++ch) {
    if (tauRes[ch].ok) {
      tr_mean_ok += tr_fix[ch];
      d_mean_ok  += d_fix[ch];
      n_ok++;
    }
  }
  if (n_ok > 0) {
    tr_mean_ok /= (double)n_ok;
    d_mean_ok  /= (double)n_ok;
    for (int ch = 0; ch < 4; ++ch) {
      if (!tauRes[ch].ok) {
        tr_fix[ch] = tr_mean_ok;
        d_fix[ch]  = d_mean_ok;
      }
    }
  }

  // ============================================================
  // Pass 2: 固定τで clean fit → 積分値ヒスト
  // ============================================================

  // ------------------------------------------------------------
  // Pass 2a: chi2/ndf カットの自動探索
  // ------------------------------------------------------------
  double chi2_cut_ch[4] = {
    kChi2NdfFallback, kChi2NdfFallback, kChi2NdfFallback, kChi2NdfFallback
  };
  long n_fit_scan[4] = {0,0,0,0};
  long n_events_scan = 0;

  if (kAutoTuneChi2Cut) {
    std::vector<double> chi2_scan[4];

    std::ifstream in[4];
    for (int ch = 0; ch < 4; ++ch) in[ch].open(fpath[ch].Data());
    bool ok = true;
    for (int ch = 0; ch < 4; ++ch) if (!in[ch].is_open()) ok = false;
    if (!ok) {
      std::cerr << "ERROR: could not open NaI files (pass2 scan)\n";
      for (int ch = 0; ch < 4; ++ch) std::cerr << "  " << fpath[ch] << "\n";
      return;
    }

    std::vector<double> buf_adc[4];
    for (int ch = 0; ch < 4; ++ch) buf_adc[ch].reserve(kSamplesPerEvent);

    while (true) {
      if (kChi2ScanMaxEvents > 0 && n_events_scan >= kChi2ScanMaxEvents) break;

      double v[4];
      bool read_ok = true;
      for (int ch = 0; ch < 4; ++ch) {
        if (!(in[ch] >> v[ch])) read_ok = false;
      }
      if (!read_ok) break;

      for (int ch = 0; ch < 4; ++ch) buf_adc[ch].push_back(v[ch]);

      bool full = true;
      for (int ch = 0; ch < 4; ++ch) if ((int)buf_adc[ch].size() != kSamplesPerEvent) full = false;
      if (!full) continue;

      bool scan_done_all = true;
      for (int c = 0; c < 4; ++c) {
        if ((int)chi2_scan[c].size() < kChi2ScanTargetFitsPerCh) {
          scan_done_all = false;
          break;
        }
      }
      if (scan_done_all) break;

      for (int ch = 0; ch < 4; ++ch) {
        const double baseline = ModeBaselineInEvent(buf_adc[ch]);
        std::vector<double> s;
        BuildS(buf_adc[ch], baseline, s);
        std::vector<int> starts = FindStartCandidates(s, kSeedThr, kEndThr, kMinSepSamples);
        AddPeakSupplementCandidates(s, starts);
        AddGlobalPeakCandidate(s, starts);
        std::vector<int> fit_cands;
        BuildFitCandidates(starts, s, kFitPreWindow, kFitHalfWindow, fit_cands);
        for (int t0 : fit_cands) {
          FitPulse fp = FitOnePulse(s, t0, tr_fix[ch], d_fix[ch], kFitHalfWindow, kT0FitRange);
          if (!fp.ok) continue;
          if (!std::isfinite(fp.chi2ndf)) continue;
          if (!PassDecompositionQuality(s, fp, tr_fix[ch], d_fix[ch])) continue;
          chi2_scan[ch].push_back(fp.chi2ndf);
        }
      }

      for (int ch = 0; ch < 4; ++ch) buf_adc[ch].clear();
      n_events_scan++;
    }

    for (int ch = 0; ch < 4; ++ch) in[ch].close();

    for (int ch = 0; ch < 4; ++ch) {
      n_fit_scan[ch] = (long)chi2_scan[ch].size();
      if (n_fit_scan[ch] >= kChi2AutoMinFits) {
        double cut = QuantileOf(chi2_scan[ch], kChi2KeepFraction);
        if (!std::isfinite(cut)) cut = kChi2NdfFallback;
        if (cut < kChi2NdfMinClamp) cut = kChi2NdfMinClamp;
        if (cut > kChi2NdfMaxClamp) cut = kChi2NdfMaxClamp;
        chi2_cut_ch[ch] = cut;
      }
    }
  }

  TH1D* hI[4];
  for (int ch = 0; ch < 4; ++ch) {
    hI[ch] = new TH1D(Form("hI_%s", kChLabel[ch]),
                      Form("%s integral (fit area);I = area [ADC*samples];Pulses", kChLabel[ch]),
                      kNBinsI, kIMin, kIMax);
  }

  long n_events = 0;
  long n_pulses_acc[4] = {0,0,0,0};

  std::vector<DebugEvent> first_dbg;
  first_dbg.reserve(kNDebugFirst);
  std::deque<DebugEvent> last_dbg;

  {
    std::ifstream in[4];
    for (int ch = 0; ch < 4; ++ch) in[ch].open(fpath[ch].Data());
    bool ok = true;
    for (int ch = 0; ch < 4; ++ch) if (!in[ch].is_open()) ok = false;
    if (!ok) {
      std::cerr << "ERROR: could not open NaI files (pass2)\n";
      for (int ch = 0; ch < 4; ++ch) std::cerr << "  " << fpath[ch] << "\n";
      return;
    }

    std::vector<double> buf_adc[4];
    for (int ch = 0; ch < 4; ++ch) buf_adc[ch].reserve(kSamplesPerEvent);

    while (true) {
      double v[4];
      bool read_ok = true;
      for (int ch = 0; ch < 4; ++ch) {
        if (!(in[ch] >> v[ch])) read_ok = false;
      }
      if (!read_ok) break;

      for (int ch = 0; ch < 4; ++ch) buf_adc[ch].push_back(v[ch]);

      bool full = true;
      for (int ch = 0; ch < 4; ++ch) if ((int)buf_adc[ch].size() != kSamplesPerEvent) full = false;
      if (!full) continue;

      DebugEvent dbg;
      dbg.idx = n_events;
      bool any_accepted_in_event = false;

      for (int ch = 0; ch < 4; ++ch) {
        dbg.adc[ch] = buf_adc[ch];

        const double baseline = ModeBaselineInEvent(buf_adc[ch]);
        dbg.baseline[ch] = baseline;

        std::vector<double> s;
        BuildS(buf_adc[ch], baseline, s);

        std::vector<int> starts = FindStartCandidates(s, kSeedThr, kEndThr, kMinSepSamples);
        AddPeakSupplementCandidates(s, starts);
        AddGlobalPeakCandidate(s, starts);
        dbg.n_cand[ch] = (int)starts.size();

        std::vector<FitPulse> accepted;
        FitPulsesWithDecomposition(s, starts, tr_fix[ch], d_fix[ch],
                                   chi2_cut_ch[ch], true, accepted);

        dbg.n_acc[ch] = (int)accepted.size();
        dbg.fit_hist[ch].clear();
        dbg.n_hist[ch] = 0;
        if (!accepted.empty()) {
          int ibest = 0;
          if (kDebugShowLargestArea) {
            for (size_t i = 1; i < accepted.size(); ++i) {
              if (accepted[i].I > accepted[ibest].I) ibest = (int)i;
            }
          }
          dbg.fit[ch] = accepted[ibest];
          any_accepted_in_event = true;

          if (kUsePrimaryPulseOnly) {
            dbg.fit_hist[ch].push_back(accepted[ibest]);
          } else {
            dbg.fit_hist[ch] = accepted;
          }

          dbg.n_hist[ch] = (int)dbg.fit_hist[ch].size();
          for (const auto& fp : dbg.fit_hist[ch]) {
            hI[ch]->Fill(fp.I);
            n_pulses_acc[ch]++;
          }
        } else {
          dbg.fit[ch] = FitPulse{};
          dbg.fit_hist[ch].clear();
          dbg.n_hist[ch] = 0;
        }
      }

      // デバッグイベント保存（受理が1つでもあれば）
      if (any_accepted_in_event) {
        if ((int)first_dbg.size() < kNDebugFirst) first_dbg.push_back(dbg);
        last_dbg.push_back(dbg);
        while ((int)last_dbg.size() > kNDebugLast) last_dbg.pop_front();
      }

      // 次イベント
      for (int ch = 0; ch < 4; ++ch) buf_adc[ch].clear();
      n_events++;
    }

    for (int ch = 0; ch < 4; ++ch) in[ch].close();
  }

  std::cout << "[p2MEG NaI] chi2 scan events = " << n_events_scan << "\n";
  for (int ch = 0; ch < 4; ++ch) {
    std::cout << "  " << kChLabel[ch]
              << "  tau(tr,d)=(" << tr_fix[ch] << ", " << d_fix[ch] << ")"
              << "  tau_ok=" << (int)tauRes[ch].ok
              << "  scan_fit_ok=" << n_fit_scan[ch]
              << "  chi2_cut=" << chi2_cut_ch[ch]
              << "  accepted=" << n_pulses_acc[ch]
              << "\n";
  }
  for (int ch = 0; ch < 4; ++ch) {
    const ShapeMetric sm = ComputeShapeMetric(hI[ch]);
    std::cout << "  " << kChLabel[ch]
              << "  shape[w1,w2,w3]=("
              << sm.n_w1 << ", " << sm.n_w2 << ", " << sm.n_w3 << ")"
              << "  bump_ratio=" << sm.bump_ratio
              << "\n";
  }
  std::cout << "[p2MEG NaI] candidate diagnostics"
            << "  input=" << gCandTotalInput
            << "  pass_for_fit=" << gCandAcceptedForFit
            << "  reject(edge/clean/pre/multi)=("
            << gCandRejectEdge << "/"
            << gCandRejectClean << "/"
            << gCandRejectPreQuiet << "/"
            << gCandRejectMultiPeak << ")"
            << "  fit_calls=" << gFitCalls
            << "\n";

  // 数値比較用に ROOT ファイルへも保存
  {
    TString out_root = Form("%s/hist_NaI_fitint_clean_%s.root", kOutputDir, kInputTag);
    TFile fout(out_root.Data(), "RECREATE");
    if (fout.IsOpen()) {
      for (int ch = 0; ch < 4; ++ch) hI[ch]->Write();
      fout.Close();
    }
  }

  // ============================================================
  // PDF 作成
  // ============================================================

  // 1ページ目: 積分ヒスト 2x2
  gStyle->SetOptStat(1110);

  TCanvas* c1 = new TCanvas("c_nai_I", "NaI integral hist", kCanvasW, kCanvasH);

  TPad* p11 = new TPad("pI11","", 0.00, 0.50, 0.50, 1.00);
  TPad* p12 = new TPad("pI12","", 0.50, 0.50, 1.00, 1.00);
  TPad* p21 = new TPad("pI21","", 0.00, 0.00, 0.50, 0.50);
  TPad* p22 = new TPad("pI22","", 0.50, 0.00, 1.00, 0.50);
  TPad* pads[4] = {p11,p12,p21,p22};
  for (int i = 0; i < 4; ++i) { SetupPad(pads[i]); pads[i]->Draw(); }

  for (int ch = 0; ch < 4; ++ch) {
    pads[ch]->cd();
    if (kUseLogY) gPad->SetLogy(1);
    hI[ch]->Draw();
    DrawLabelTopLeft(kChLabel[ch]);
  }

  TString out_open = out_pdf; out_open += "(";
  c1->SaveAs(out_open.Data());

  // 2ページ目: メタ情報
  gStyle->SetOptStat(0);

  TCanvas* c2 = new TCanvas("c_nai_meta", "meta", kCanvasW, kCanvasH);
  c2->cd();

  TPaveText* pt = new TPaveText(0.04, 0.06, 0.96, 0.94, "NDC");
  pt->SetFillStyle(0);
  pt->SetBorderSize(0);
  pt->SetTextAlign(12);
  pt->SetTextSize(0.035);

  TDatime now;
  pt->AddText("=== NaI clean fit integral meta ===");
  pt->AddText(Form("InputDir : %s", kInputDir));
  pt->AddText(Form("Tag      : %s", kInputTag));
  pt->AddText(Form("Output   : %s", out_pdf.Data()));
  pt->AddText(" ");
  pt->AddText("Baseline: mode(ADC) within each digitizer event (per channel).");
  pt->AddText("Model (s = baseline - ADC):  C + A*(exp(-(t-t0)/(tr+d)) - exp(-(t-t0)/tr))*Theta(t-t0)");
  pt->AddText(Form("Baseline in fit: C %s", kFixFitBaselineToZero ? "fixed to 0" : "free"));
  pt->AddText("Integral (area): I = A * d   (units: ADC*samples; multiply by dt to get ADC*ns)");
  pt->AddText(" ");
  pt->AddText(Form("samples_per_event=%d, dt=%.1f ns", kSamplesPerEvent, kDtNs));
  pt->AddText(Form("fit window: pre=%d samp (%.0f ns), post=%d samp (%.0f ns)",
                   kFitPreWindow, kFitPreWindow*kDtNs,
                   kFitHalfWindow, kFitHalfWindow*kDtNs));
  pt->AddText(Form("candidate: seed_thr=%.1f, end_thr=%.1f, min_sep=%d samples", kSeedThr, kEndThr, kMinSepSamples));
  pt->AddText(Form("clean: require single candidate in fit window = %d", (int)kCleanRequireSingleInWindow));
  pt->AddText(Form("fit: t0 range +-%.1f samples (adaptive=%d, x%.1f)",
                   kT0FitRange, (int)kAdaptiveT0Range, kT0FitRangeExpandFactor));
  pt->AddText(Form("fit quality guard: |dt0| edge reject=%.2f*range, t0_fit>=%.1f, chi2_cut<=%.1f",
                   kT0BoundaryFracReject, kT0FitAbsMin, kChi2NdfMaxClamp));
  pt->AddText(Form("candidate handling: primary-only=%d, max-fit-cands=%d, debug-largest-area=%d",
                   (int)kUsePrimaryPulseOnly, kMaxCandidatesToFitPrimary, (int)kDebugShowLargestArea));
  pt->AddText(Form("candidate prefilter: pre-quiet(abs<=%.1f, frac<=%.2f), multi-peak-veto=%d (2nd/1st<=%.2f)",
                   kCandPreQuietAbsMax, kCandPreQuietFracMax,
                   (int)kCandRejectMultiPeak, kCandSecondPeakFracMax));
  pt->AddText(Form("candidate supplement: peak-based=%d (h>=%.1f, prom>=%.1f, min_sep=%d)",
                   (int)kUsePeakSupplementCand, kPeakCandMinHeight, kPeakCandMinProminence, kPeakCandMinSep));
  pt->AddText(Form("primary consistency: enabled=%d (fit_peak >= %.2f * event_peak, primary-only mode)",
                   (int)kRequirePrimaryPeakConsistency, kPrimaryPeakFracMin));
  pt->AddText(Form("candidate diagnostics: input=%ld, pass_for_fit=%ld, reject(edge/clean/pre/multi)=(%ld/%ld/%ld/%ld), fit_calls=%ld",
                   gCandTotalInput, gCandAcceptedForFit,
                   gCandRejectEdge, gCandRejectClean, gCandRejectPreQuiet, gCandRejectMultiPeak,
                   gFitCalls));
  pt->AddText(" ");
  pt->AddText("Tau estimation: average of locally-isolated clean pulses, then fit tr,d and fix.");
  for (int ch = 0; ch < 4; ++ch) {
    if (tauRes[ch].ok) {
      pt->AddText(Form("  %-6s : used=%d  tr=%.2f samp (%.1f ns)  d=%.2f samp (%.1f ns)  chi2/ndf=%.2f",
                       kChLabel[ch], tauRes[ch].n_used,
                       tr_fix[ch], tr_fix[ch]*kDtNs,
                       d_fix[ch],  d_fix[ch]*kDtNs,
                       tauRes[ch].chi2ndf));
    } else {
      pt->AddText(Form("  %-6s : used=%d  -> fallback tr=%.2f samp, d=%.2f samp",
                       kChLabel[ch], tauRes[ch].n_used, tr_fix[ch], d_fix[ch]));
    }
  }
  pt->AddText(" ");
  if (kAutoTuneChi2Cut) {
    pt->AddText(Form("chi2 tuning: scan first %ld events (target/ch=%d), keep best %.0f%% fits",
                     n_events_scan, kChi2ScanTargetFitsPerCh, 100.0 * kChi2KeepFraction));
  } else {
    pt->AddText("chi2 tuning: disabled (fallback cut only)");
  }
  for (int ch = 0; ch < 4; ++ch) {
    pt->AddText(Form("  %-6s : scan fit_ok=%ld  cut chi2/ndf <= %.2f",
                     kChLabel[ch], n_fit_scan[ch], chi2_cut_ch[ch]));
  }
  pt->AddText(" ");
  pt->AddText(Form("processed digitizer events = %ld", n_events));
  for (int ch = 0; ch < 4; ++ch) {
    pt->AddText(Form("  %-6s : accepted pulses = %ld", kChLabel[ch], n_pulses_acc[ch]));
  }
  pt->AddText(Form("Generated: %04d-%02d-%02d %02d:%02d:%02d",
                   now.GetYear(), now.GetMonth(), now.GetDay(),
                   now.GetHour(), now.GetMinute(), now.GetSecond()));
  pt->Draw();
  c2->SaveAs(out_pdf.Data());

  // 3ページ目: 平均波形と τ フィット（2x2）
  TCanvas* c3 = new TCanvas("c_nai_tau", "tau estimation", kCanvasW, kCanvasH);

  TPad* q11 = new TPad("q11","", 0.00, 0.50, 0.50, 1.00);
  TPad* q12 = new TPad("q12","", 0.50, 0.50, 1.00, 1.00);
  TPad* q21 = new TPad("q21","", 0.00, 0.00, 0.50, 0.50);
  TPad* q22 = new TPad("q22","", 0.50, 0.00, 1.00, 0.50);
  TPad* qpad[4] = {q11,q12,q21,q22};
  for (int i = 0; i < 4; ++i) { SetupPad(qpad[i]); qpad[i]->Draw(); }

  for (int ch = 0; ch < 4; ++ch) {
    qpad[ch]->cd();

    if (n_used[ch] < kNTauEstMin) {
      TLatex t;
      t.SetNDC(true);
      t.SetTextSize(0.06);
      t.DrawLatex(0.12, 0.80, Form("%s", kChLabel[ch]));
      t.SetTextSize(0.05);
      t.DrawLatex(0.12, 0.70, Form("tau estimation: insufficient clean pulses (used=%d)", n_used[ch]));
      continue;
    }

    TGraph* gr = new TGraph(n_dt, xdt.data(), avg[ch].data());
    gr->SetTitle(Form("%s avg(clean pulses);dt [samples];avg s [counts]", kChLabel[ch]));
    gr->Draw("AL");

    TF1 f("f_tau_draw", NaIPulseFunc_ExpDiff, xdt.front(), xdt.back(), 5);
    f.SetParameters(1.0, 0.0, 0.0, tr_fix[ch], d_fix[ch]); // Aは適当にスケールして見せる
    // A はグラフに合わせて（最大値/ gmax）にする
    double ymax = 0.0;
    for (double y : avg[ch]) if (y > ymax) ymax = y;
    const double gmax = Gmax_ExpDiff(tr_fix[ch], d_fix[ch], kTauEstHalfWindow);
    const double A_show = (gmax > 0) ? (ymax / gmax) : ymax;
    f.SetParameter(0, A_show);
    f.SetParameter(1, 0.0);
    f.SetParameter(2, 0.0);
    f.SetParameter(3, tr_fix[ch]);
    f.SetParameter(4, d_fix[ch]);
    f.SetLineWidth(2);
    f.Draw("SAME");

    TLatex t;
    t.SetNDC(true);
    t.SetTextSize(0.05);
    t.DrawLatex(0.14, 0.84, Form("%s", kChLabel[ch]));
    t.SetTextSize(0.04);
    t.DrawLatex(0.14, 0.78, Form("used=%d  tr=%.2f samp  d=%.2f samp",
                                 n_used[ch], tr_fix[ch], d_fix[ch]));
  }

  c3->SaveAs(out_pdf.Data());

  // ============================================================
  // デバッグページ（PDF最後）：先頭5 + 末尾5
  // ============================================================

  std::vector<DebugEvent> dbg_list;
  dbg_list.reserve(kNDebugFirst + kNDebugLast);
  for (const auto& e : first_dbg) dbg_list.push_back(e);
  for (const auto& e : last_dbg)  dbg_list.push_back(e);

  if (dbg_list.empty()) {
    // デバッグ無しでもPDFを閉じる
    TString out_close = out_pdf; out_close += ")";
    c3->SaveAs(out_close.Data());
    std::cout << "Saved: " << out_pdf << "\n";
    return;
  }

  for (size_t i = 0; i < dbg_list.size(); ++i) {
    const bool is_last = (i == dbg_list.size() - 1);
    DrawDebugPage(dbg_list[i], tr_fix, d_fix, out_pdf, is_last);
  }

  std::cout << "Saved: " << out_pdf << "\n";
}
