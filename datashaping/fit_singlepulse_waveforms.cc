#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <cctype>
#include <cerrno>
#include <cstdlib>
#include <random>
#include <iomanip>
#include <filesystem>

#include "TCanvas.h"
#include "TGraph.h"
#include "TF1.h"
#include "TH1D.h"

// ============================================================
// p2MEG: 単パルス・テンプレート構築 (ROOT)
//
// 目的:
//  1) 単パルス選別
//  2) ベースライン推定（MAD最小窓の中央値）
//  3) leading-edge 相互相関で整列
//  4) 平均テンプレート作成（min=-1 で正規化）
//  5) 指数パルス + ガウス畳み込みフィット
//  6) 可視化
//
// 注意:
//  - 入力波形は負パルスを想定
//  - サンプル間隔は一定
//  - 速度より頑健性を優先
// ============================================================

// ------------------------------
// 既定パラメータ
// ------------------------------
static const int kBaselineWinSamples = 100;   // ベースライン探索窓の長さ [sample]
static const double kSnrThreshold = 8.0;       // S/N カット
// 相互相関窓: t0_index 基準（[t0-kAlignWinPreSamples, t0+kAlignWinPostSamples)）
static const int kAlignWinPreSamples = 15;
static const int kAlignWinPostSamples = 15;
static const int kMaxShift = 15;               // 許容シフト
static const double kDerivThreshold = 2.5;     // 2つ目パルスの導関数閾値 (×MAD換算σ)
static const double kRecoveryDip = 2.5;        // 回復後再ディップ閾値 (×MAD換算σ)
static const double kSaturationMargin = 2.0;   // 飽和判定の余裕 [ADC]
static const bool kUseSaturationCheck = false; // 既知の ADC 範囲がない場合は false
static const double kFitFrac = 0.20;           // フィット窓 (min の割合)
static const int kConvHalfWin = 60;            // 畳み込みの半幅 [sample]
static const double kSnrBottomFrac = 0.15;     // S/N カット通過後に下位から除外する割合
static const int kMinPrePulseSeparation = 200;  // [sample] 前パルスが主パルスに近すぎる判定
static const double kPrePulseThrSigma = 6.0;   // 前パルス検出しきい値（×MAD換算σ）
static const int kBaselineRecoverLen = 15;     // 回復判定: baseline近傍が連続するサンプル数
static const double kBaselineRecoverSigma = 2.5; // baseline近傍判定: |y| < kBaselineRecoverSigma*(MAD換算σ)
static const int kPreIgnoreBeforeMin =50;     // 最小値直前の立ち下がり領域は前パルス判定から除外
static const int kMinPreSamples = 50;         // 最小値の前に必要な最小サンプル数
static const int kMinPostSamples = 200;        // 最小値の後ろに必要な最小サンプル数
static const double kCoarseStartSigma = 2.5;   // 立ち下がり開始検出しきい値（×MAD換算σ）
static const int kCoarseLookback = 120;        // 最小点からどのくらい前まで立ち下がり探索するか [sample]
static const int kCoarseTargetIndex = 120;     // 立ち下がり開始位置を揃える目標サンプル
static const int kCoarseStartConsec = 3;       // 立ち下がり開始: 連続N点で閾値超え判定
static const double kLeadReboundSigma = 1.0; // 立ち下がり途中の“戻り”判定（×MAD換算σ）
static const double kLeadDipSigma = 1.0;     // 立ち下がり途中の“再ディップ”判定（×MAD換算σ）
static const double kSampleDt = 1.0;          // サンプル間隔（未定義の場合は 1.0 とする）
static const char* kNormDescription = "min=-1 (event-wise, then averaged, final renorm)";

// ------------------------------
// ユーティリティ
// ------------------------------
static double Median(std::vector<double> v) {
    std::nth_element(v.begin(), v.begin() + v.size()/2, v.end());
    return v[v.size()/2];
}

static double ComputeRMS(const std::vector<double>& v) {
    if (v.empty()) return 0.0;
    double mean = std::accumulate(v.begin(), v.end(), 0.0) / v.size();
    double s2 = 0.0;
    for (auto &x : v) {
        double d = x - mean;
        s2 += d * d;
    }
    return std::sqrt(s2 / v.size());
}

static double ComputeMadSigma(const std::vector<double>& v) {
    if (v.empty()) return 0.0;
    std::vector<double> tmp = v;
    double med = Median(tmp);
    std::vector<double> dev;
    dev.reserve(v.size());
    for (auto x : v) dev.push_back(std::abs(x - med));
    double mad = Median(dev);
    return 1.4826 * mad; // 正規分布換算（MAD→σ相当）
}

// ------------------------------
// 入力ファイル名から run / module を推定
// 例: wave_PS_A_run7.txt, wave_NaI_B1_run7.txt
// ------------------------------
static bool ParseInputMetaFromPath(const std::string& path, int& run_no, std::string& ch_name) {
    run_no = -1;
    ch_name = "unknown";
    std::string base = std::filesystem::path(path).filename().string();

    // run 番号: "run" の後ろの数字
    size_t prun = base.rfind("run");
    if (prun != std::string::npos) {
        size_t p = prun + 3;
        size_t q = p;
        while (q < base.size() && std::isdigit(static_cast<unsigned char>(base[q]))) ++q;
        if (q > p) {
            run_no = std::atoi(base.substr(p, q - p).c_str());
        }
    }

    // CH: "wave_" または "Michel_" の直後〜"_run" の直前
    size_t pmod = std::string::npos;
    if (base.rfind("wave_", 0) == 0) {
        pmod = std::string("wave_").size();
    } else if (base.rfind("Michel_", 0) == 0) {
        pmod = std::string("Michel_").size();
    }
    if (pmod != std::string::npos) {
        size_t end = base.find("_run", pmod);
        if (end != std::string::npos && end > pmod) {
            ch_name = base.substr(pmod, end - pmod);
        }
    }

    return (run_no >= 0);
}

// 前方宣言（t0 推定で使用）
static int FindFallingStartIndex(const std::vector<double>& y, double noise_mad, int idx_min);

// ------------------------------
// 波形CSV出力
// ------------------------------
static void WriteWaveformsCsv(const std::string& path,
                              const std::vector<std::vector<double>>& events) {
    std::error_code ec;
    std::filesystem::create_directories("data/shapeddata", ec);
    if (ec) {
        std::cerr << "Failed to create output dir: data/shapeddata (" << ec.message() << ")\n";
        return;
    }
    std::ofstream ofs(path);
    if (!ofs) {
        std::cerr << "Failed to open output file: " << path << std::endl;
        return;
    }
    for (size_t i = 0; i < events.size(); ++i) {
        ofs << i;
        const auto& y = events[i];
        for (size_t j = 0; j < y.size(); ++j) {
            ofs << "," << y[j];
        }
        ofs << "\n";
    }
}

// ------------------------------
// 入力読み込み
// ------------------------------
static std::vector<double> LoadBins(const std::string& path) {
    std::ifstream ifs(path);
    if (!ifs) throw std::runtime_error("Failed to open input file: " + path);

    std::vector<double> bins;
    bins.reserve(1 << 20);

    std::string line;
    long long line_no = 0;
    while (std::getline(ifs, line)) {
        ++line_no;
        auto is_space = [](unsigned char c) { return std::isspace(c); };
        line.erase(line.begin(), std::find_if(line.begin(), line.end(),
                    [&](unsigned char c){ return !is_space(c); }));
        line.erase(std::find_if(line.rbegin(), line.rend(),
                    [&](unsigned char c){ return !is_space(c); }).base(), line.end());
        if (line.empty()) continue;

        char* endptr = nullptr;
        errno = 0;
        double v = std::strtod(line.c_str(), &endptr);
        if (errno != 0 || endptr == line.c_str()) {
            throw std::runtime_error(
                "Non-numeric line in " + path + " at line " + std::to_string(line_no) + ": " + line
            );
        }
        bins.push_back(v);
    }
    return bins;
}

static std::vector<std::vector<double>> SplitEvents(const std::vector<double>& all_bins, int bins_per_event) {
    if (bins_per_event <= 0) throw std::runtime_error("bins_per_event must be > 0");

    const size_t n = all_bins.size();
    const size_t N = static_cast<size_t>(bins_per_event);
    if (n % N != 0) {
        throw std::runtime_error(
            "Total bins (" + std::to_string(n) + ") is not divisible by bins_per_event (" +
            std::to_string(bins_per_event) + ")."
        );
    }

    const size_t nevt = n / N;
    std::vector<std::vector<double>> events;
    events.reserve(nevt);

    for (size_t e = 0; e < nevt; ++e) {
        auto start = all_bins.begin() + static_cast<long long>(e * N);
        auto end   = start + static_cast<long long>(N);
        events.emplace_back(start, end);
    }
    return events;
}

// ------------------------------
// ベースライン推定（MAD最小窓を探し、その窓の中央値を使う）
// ------------------------------
static double EstimateBaseline(const std::vector<double>& p, double *noise_rms, double *noise_mad) {
    const int N = static_cast<int>(p.size());
    const int nwin = std::min(N, std::max(1, kBaselineWinSamples));

    // MAD が最小の「静かな窓」を探す（信号が無い区間を想定）
    int best_start = 0;
    double best_mad = 1e99;
    for (int start = 0; start + nwin <= N; ++start) {
        std::vector<double> w(p.begin() + start, p.begin() + start + nwin);
        double mad = ComputeMadSigma(w);
        if (mad < best_mad) {
            best_mad = mad;
            best_start = start;
        }
    }
    std::vector<double> head(p.begin() + best_start, p.begin() + best_start + nwin);

    // 選ばれた窓の RMS / MAD(換算σ)
    const double rms = ComputeRMS(head);
    const double mad = ComputeMadSigma(head);
    if (noise_rms) *noise_rms = rms;
    if (noise_mad) *noise_mad = mad;

    // 窓内の中央値を baseline とする
    return Median(head);
}

// ------------------------------
// 相互相関（leading-edge の短窓）
// ------------------------------
static int FindShiftByCorrelation(
    const std::vector<double>& y,
    const std::vector<double>& ref,
    int win_start,
    int win_end,
    int max_shift
) {
    double best = -1e99;
    int best_shift = 0;

    for (int s = -max_shift; s <= max_shift; ++s) {
        double c = 0.0;
        for (int i = win_start; i < win_end; ++i) {
            int j = i + s;
            if (j < 0 || j >= (int)ref.size()) continue;
            c += y[i] * ref[j];
        }
        if (c > best) {
            best = c;
            best_shift = s;
        }
    }
    return best_shift;
}

// 波形シフト
static std::vector<double> ShiftWave(const std::vector<double>& y, int shift) {
    std::vector<double> out(y.size(), 0.0);
    for (size_t i = 0; i < y.size(); ++i) {
        int j = (int)i + shift;
        if (j >= 0 && j < (int)y.size()) out[j] = y[i];
    }
    return out;
}

static std::pair<double, double> MinMaxVec(const std::vector<double>& y) {
    if (y.empty()) return {0.0, 0.0};
    auto mm = std::minmax_element(y.begin(), y.end());
    return {*mm.first, *mm.second};
}

static std::pair<double, double> ExpandRange(double ymin, double ymax, double frac = 0.05) {
    double span = ymax - ymin;
    if (span <= 0) span = 1.0;
    double pad = span * frac;
    return {ymin - pad, ymax + pad};
}

// t0_index 基準で相関窓 [start, end) を作る
static void BuildAlignWindowFromT0(int t0_index, int n_samples, int& win_start, int& win_end) {
    win_start = t0_index - kAlignWinPreSamples;
    win_end = t0_index + kAlignWinPostSamples;
    if (win_start < 0) win_start = 0;
    if (win_end > n_samples) win_end = n_samples;
    if (win_end <= win_start + 1) {
        win_start = 0;
        win_end = std::min(n_samples, 2);
    }
}

// 傾きが負に最大（最も急な立ち下がり）の点を t0 とする
static int FindTemplateT0ByMaxNegSlope(const std::vector<double>& templ) {
    if (templ.size() < 2) return 0;
    double min_dy = 1e99;
    int idx = 0;
    for (int i = 0; i + 1 < static_cast<int>(templ.size()); ++i) {
        double dy = templ[i + 1] - templ[i];
        if (dy < min_dy) {
            min_dy = dy;
            idx = i + 1;
        }
    }
    return idx;
}
// ------------------------------
// 出力（テンプレート JSON / ベースライン CSV）
// ------------------------------
static void WriteTemplateJson(
    const std::string& path,
    const std::vector<double>& templ,
    double dt,
    const char* norm_desc,
    int t0_index,
    TF1* f,
    int run_no,
    const std::string& ch_name
) {
    std::ofstream ofs(path);
    if (!ofs) {
        std::cerr << "Failed to open output file: " << path << std::endl;
        return;
    }

    double sum = 0.0;
    double sum2 = 0.0;
    for (auto v : templ) {
        sum += v;
        sum2 += v * v;
    }
    double areaT = sum * dt;
    double sumT2 = sum2;

    auto mm = MinMaxVec(templ);
    const char* polarity = (mm.first < 0.0) ? "negative" : "positive";

    ofs << std::setprecision(12);
    ofs << "{\n";
    ofs << "  \"dt\": " << dt << ",\n";
    ofs << "  \"run\": " << run_no << ",\n";
    ofs << "  \"CH\": \"" << ch_name << "\",\n";
    ofs << "  \"polarity\": \"" << polarity << "\",\n";
    ofs << "  \"norm\": \"" << norm_desc << "\",\n";
    ofs << "  \"t0_index\": " << t0_index << ",\n";
    ofs << "  \"template\": [";
    for (size_t i = 0; i < templ.size(); ++i) {
        if (i != 0) ofs << ", ";
        ofs << templ[i];
    }
    ofs << "],\n";
    ofs << "  \"areaT\": " << areaT << ",\n";
    ofs << "  \"sumT2\": " << sumT2 << ",\n";
    ofs << "  \"pulseModelFit\": {\n";
    if (f) {
        ofs << "    \"A\": " << f->GetParameter(0) << ",\n";
        ofs << "    \"t0\": " << f->GetParameter(1) << ",\n";
        ofs << "    \"tau_r\": " << f->GetParameter(2) << ",\n";
        ofs << "    \"tau_d\": " << f->GetParameter(3) << ",\n";
        ofs << "    \"sigma\": " << f->GetParameter(4) << "\n";
    } else {
        ofs << "    \"A\": 0,\n";
        ofs << "    \"t0\": 0,\n";
        ofs << "    \"tau_r\": 0,\n";
        ofs << "    \"tau_d\": 0,\n";
        ofs << "    \"sigma\": 0\n";
    }
    ofs << "  }\n";
    ofs << "}\n";
}

static void WriteBaselineCsv(
    const std::string& path,
    const std::vector<double>& baseline,
    const std::vector<double>& noise_mad,
    int run_no,
    const std::string& ch_name
) {
    std::ofstream ofs(path);
    if (!ofs) {
        std::cerr << "Failed to open output file: " << path << std::endl;
        return;
    }
    ofs << std::setprecision(12);
    ofs << "run,event_id,CH,baseline,noise\n";
    const size_t n = std::min(baseline.size(), noise_mad.size());
    for (size_t i = 0; i < n; ++i) {
        ofs << run_no << "," << i << "," << ch_name
            << "," << baseline[i] << "," << noise_mad[i] << "\n";
    }
}

// 立ち下がり開始（負パルス）を検出して粗い揃え用インデックスを返す
// idx_min に紐づけて主パルス側（idx_min より前）の範囲で探す
static int FindFallingStartIndex(const std::vector<double>& y, double noise_mad, int idx_min) {
    if (noise_mad <= 0) return -1;
    if (idx_min <= 0) return -1;
    const double thr = -kCoarseStartSigma * noise_mad;
    const int N = static_cast<int>(y.size());
    const int end = std::min(N - 1, idx_min);
    const int start = std::max(0, end - kCoarseLookback);
    // idx_min から逆向きに「閾値クロス」を探す
    // 直前が閾値以上で、次が閾値未満になった点を立ち下がり開始とみなす
    const int kConsec = kCoarseStartConsec; // 連続N点で閾値超え判定
    for (int i = end - 1; i >= start; --i) {
        bool ok = true;
        for (int j = 1; j <= kConsec; ++j) {
            if (i + j > end) { ok = false; break; }
            if (!(y[i + j] < thr)) { ok = false; break; }
        }
        if (ok && y[i] >= thr) return i + 1;
    }
    return -1;
}

// 粗シフト後に相関窓が有効範囲内かどうか
static bool IsCorrelationWindowValid(int shift, int n_samples, int t0_index) {
    int win_start = 0;
    int win_end = 0;
    BuildAlignWindowFromT0(t0_index, n_samples, win_start, win_end);
    const int win_last = win_end - 1;
    const int src_start = win_start - shift;
    const int src_end = win_last - shift;
    return (src_start >= 0 && src_end < n_samples);
}

// ------------------------------
// 単パルス判定
// ------------------------------
static bool IsSinglePulse(
    const std::vector<double>& y,
    double noise_rms,
    double noise_mad,
    double adc_min,
    double adc_max,
    int* fail_code = nullptr
) {
    if (y.empty()) return false;

    auto it_min = std::min_element(y.begin(), y.end());
    double minv = *it_min;
    int idx_min = static_cast<int>(it_min - y.begin());

    // -2) 最小値より前のサンプル数不足を除外
    if (idx_min < kMinPreSamples) {
        if (fail_code) *fail_code = 10;
        return false;
    }

    // -1) 最小値以降のサンプル数不足を除外
    if (static_cast<int>(y.size()) - 1 - idx_min < kMinPostSamples) {
        if (fail_code) *fail_code = 8;
        return false;
    }

    // a) S/N カット（最も軽量）
    if (std::abs(minv) <= kSnrThreshold * noise_rms) {
        if (fail_code) *fail_code = 1;
        return false;
    }

    // テンプレフィット相当の範囲（y<kFitFrac*minv の最初/最後）
    double thr = kFitFrac * minv; // minv は負
    int fit_start = 0;
    int fit_end = static_cast<int>(y.size()) - 1;
    for (int i = 0; i < static_cast<int>(y.size()); ++i) {
        if (y[i] < thr) { fit_start = i; break; }
    }
    for (int i = static_cast<int>(y.size()) - 1; i >= 0; --i) {
        if (y[i] < thr) { fit_end = i; break; }
    }

    // 0) 前パルス許容判定（十分離れて回復していればOK）
    if (noise_mad > 0) {
        int last_pre_idx = -1;
        double pre_thr = -kPrePulseThrSigma * noise_mad;
        int pre_end = std::max(0, idx_min - kPreIgnoreBeforeMin);
        for (int i = 0; i < pre_end; ++i) {
            if (y[i] < pre_thr) last_pre_idx = i;
        }

        if (last_pre_idx >= 0) {
            int dt = idx_min - last_pre_idx;
            if (dt < kMinPrePulseSeparation) {
                if (fail_code) *fail_code = 6;
                return false;
            }

            // 回復判定: baseline近傍が連続する区間があるか
            int consec = 0;
            double rec_thr = kBaselineRecoverSigma * noise_mad;
            for (int i = last_pre_idx + 1; i < pre_end; ++i) {
                if (std::abs(y[i]) < rec_thr) {
                    ++consec;
                    if (consec >= kBaselineRecoverLen) break;
                } else {
                    consec = 0;
                }
            }
            if (consec < kBaselineRecoverLen) {
                if (fail_code) *fail_code = 7;
                return false;
            }
        }
    }

    // 0.5) 立ち下がり途中の「戻り」＋「再ディップ」チェック
    if (noise_mad > 0) {
        int t_le = FindFallingStartIndex(y, noise_mad, idx_min);
        if (t_le >= 0 && idx_min > t_le + 2) {
            double running_min = y[t_le];
            double max_rebound = 0.0;
            bool rebound_seen = false;
            double rebound_peak = y[t_le];
            for (int i = t_le + 1; i <= idx_min; ++i) {
                if (y[i] < running_min) running_min = y[i];
                double rebound = y[i] - running_min; // 上に戻る量（正）
                if (rebound > max_rebound) max_rebound = rebound;
                if (rebound > kLeadReboundSigma * noise_mad) {
                    rebound_seen = true;
                    if (y[i] > rebound_peak) rebound_peak = y[i];
                }
                if (rebound_seen && (rebound_peak - y[i]) > kLeadDipSigma * noise_mad) {
                    if (fail_code) *fail_code = 12;
                    return false;
                }
            }
            if (max_rebound > kLeadReboundSigma * noise_mad) {
                if (fail_code) *fail_code = 9;
                return false;
            }
        }
    }

    // b) 飽和除外（ADC 範囲が既知のときのみ有効化）
    if (kUseSaturationCheck) {
        if (minv <= adc_min + kSaturationMargin) {
            if (fail_code) *fail_code = 2;
            return false;
        }
    }

    // c) 二つ目パルスの導関数チェック
    for (int i = idx_min + 1; i + 1 <= fit_end; ++i) {
        double dy = y[i + 1] - y[i];
        if (dy < -kDerivThreshold * noise_mad) {
            if (fail_code) *fail_code = 3;
            return false;
        }
    }

    // d) 回復後の再ディップ
    double max_after = y[idx_min];
    for (int i = idx_min + 1; i <= fit_end; ++i) {
        if (y[i] > max_after) max_after = y[i];
        if (y[i] < max_after - kRecoveryDip * noise_mad) {
            if (fail_code) *fail_code = 4;
            return false;
        }
    }


    // 上側の飽和も除外
    double maxv = *std::max_element(y.begin(), y.end());
    if (kUseSaturationCheck) {
        if (maxv >= adc_max - kSaturationMargin) {
            if (fail_code) *fail_code = 5;
            return false;
        }
    }

    return true;
}

// ------------------------------
// 単指数 + ガウス畳み込み
// ------------------------------
static double pulseModel(double *x, double *p) {
    // p[0]=A, p[1]=t0, p[2]=tau_r, p[3]=tau_d, p[4]=sigma
    double t = x[0];
    double A = p[0];
    double t0 = p[1];
    double tau_r = p[2];
    double tau_d = p[3];
    double sigma = p[4];

    // 数値畳み込み（サンプル間隔=1 を仮定）
    double sum = 0.0;
    double norm = 0.0;
    const double sig2 = sigma * sigma;

    for (int k = -kConvHalfWin; k <= kConvHalfWin; ++k) {
        double tt = t - k;
        double h = 0.0;
        if (tt >= t0) {
            h = -(std::exp(-(tt - t0)/tau_d) - std::exp(-(tt - t0)/tau_r));
        }
        double g = (sig2 > 0) ? std::exp(-0.5 * k * k / sig2) : (k == 0 ? 1.0 : 0.0);
        sum += h * g;
        norm += g;
    }

    if (norm <= 0.0) return 0.0;
    return A * sum / norm;
}

// ------------------------------
// メイン処理
// ------------------------------
void singlePulseTemplateDemo(const std::vector<std::vector<double>>& pulses,
                             int run_no,
                             const std::string& ch_name) {
    if (pulses.empty()) return;
    const int N = static_cast<int>(pulses[0].size());

    // ADC の上限/下限を観測データから推定（飽和検出用）
    double adc_min = 1e99;
    double adc_max = -1e99;
    for (auto &p : pulses) {
        auto mm = std::minmax_element(p.begin(), p.end());
        if (*mm.first < adc_min) adc_min = *mm.first;
        if (*mm.second > adc_max) adc_max = *mm.second;
    }

    std::vector<std::vector<double>> cleaned;   // baseline 除去済み
    std::vector<double> cleaned_rms;            // cleaned に対応する noise_rms
    std::vector<double> cleaned_mad;            // cleaned に対応する noise_mad
    std::vector<std::vector<double>> aligned;   // 整列済み（元スケール）
    std::vector<std::vector<double>> aligned_norm; // 整列後にイベント毎正規化（min=-1）
    std::vector<std::vector<double>> rejected;  // 除外波形
    std::vector<double> noise_rms_list;
    std::vector<double> noise_mad_list;
    std::vector<double> minv_list;
    std::vector<double> snr_list;
    long long fail_snr = 0;
    long long fail_snr_bottom = 0;
    long long fail_sat = 0;
    long long fail_deriv = 0;
    long long fail_recover = 0;
    long long fail_sat_hi = 0;
    long long fail_pre_close = 0;
    long long fail_pre_norec = 0;
    long long fail_post_short = 0;
    long long fail_le_bump = 0;
    long long fail_pre_short = 0;
    long long fail_le_dip = 0;

    // Step1-1: ベースライン推定 + S/N の事前計算
    std::vector<std::vector<double>> y_all;
    std::vector<double> noise_rms_all;
    std::vector<double> noise_mad_all;
    std::vector<double> snr_all;
    std::vector<double> minv_all;
    std::vector<double> baseline_all;
    y_all.reserve(pulses.size());
    noise_rms_all.reserve(pulses.size());
    noise_mad_all.reserve(pulses.size());
    snr_all.reserve(pulses.size());
    minv_all.reserve(pulses.size());
    baseline_all.reserve(pulses.size());

    for (auto &p : pulses) {
        double noise_rms = 0.0;
        double noise_mad = 0.0;
        double base = EstimateBaseline(p, &noise_rms, &noise_mad);

        std::vector<double> y = p;
        for (auto &v : y) v -= base;

        double minv = *std::min_element(y.begin(), y.end());
        double snr = (noise_mad > 0) ? std::abs(minv) / noise_mad : 0.0;
        noise_rms_list.push_back(noise_rms);
        noise_mad_list.push_back(noise_mad);
        minv_list.push_back(minv);
        snr_list.push_back(snr);

        y_all.push_back(std::move(y));
        noise_rms_all.push_back(noise_rms);
        noise_mad_all.push_back(noise_mad);
        snr_all.push_back(snr);
        minv_all.push_back(minv);
        baseline_all.push_back(base);
    }

    // 全イベントの baseline / noise を出力
    WriteBaselineCsv("data/shapeddata/event_baseline_noise.csv", baseline_all, noise_mad_all,
                     run_no, ch_name);

    // Step1-2: S/N カット通過後の下位10%を除外する閾値を決める
    double snr_bottom_cut = -1.0;
    std::vector<double> snr_pass;
    snr_pass.reserve(snr_all.size());
    for (double s : snr_all) {
        if (s > kSnrThreshold) snr_pass.push_back(s);
    }
    if (!snr_pass.empty()) {
        std::sort(snr_pass.begin(), snr_pass.end());
        size_t idx = static_cast<size_t>(std::floor(kSnrBottomFrac * (snr_pass.size() - 1)));
        snr_bottom_cut = snr_pass[idx];
    }

    // Step1-3: 単パルス選別
    for (size_t i = 0; i < y_all.size(); ++i) {
        const std::vector<double>& y = y_all[i];
        double noise_rms = noise_rms_all[i];
        double noise_mad = noise_mad_all[i];
        double snr = snr_all[i];

        // まず S/N カット
        if (snr <= kSnrThreshold) {
            ++fail_snr;
            rejected.push_back(y);
            continue;
        }
        // 次に下位10%を除外
        if (snr_bottom_cut >= 0 && snr < snr_bottom_cut) {
            ++fail_snr_bottom;
            rejected.push_back(y);
            continue;
        }

        int fail_code = 0;
        if (IsSinglePulse(y, noise_rms, noise_mad, adc_min, adc_max, &fail_code)) {
            cleaned.push_back(y);
            cleaned_rms.push_back(noise_rms);
            cleaned_mad.push_back(noise_mad);
        } else {
            rejected.push_back(y);
            if (fail_code == 1) ++fail_snr;
            else if (fail_code == 2) ++fail_sat;
            else if (fail_code == 3) ++fail_deriv;
            else if (fail_code == 4) ++fail_recover;
            else if (fail_code == 5) ++fail_sat_hi;
            else if (fail_code == 6) ++fail_pre_close;
            else if (fail_code == 7) ++fail_pre_norec;
            else if (fail_code == 8) ++fail_post_short;
            else if (fail_code == 9) ++fail_le_bump;
            else if (fail_code == 10) ++fail_pre_short;
            else if (fail_code == 12) ++fail_le_dip;
        }
    }

    std::cout << "Total events: " << pulses.size()
              << "  selected: " << cleaned.size()
              << "  rejected: " << rejected.size() << std::endl;
    if (!noise_rms_list.empty()) {
        double rms_med = Median(noise_rms_list);
        double mad_med = Median(noise_mad_list);
        double minv_med = Median(minv_list);
        double snr_med = Median(snr_list);
        std::cout << "Median noise RMS: " << rms_med
                  << "  median MAD: " << mad_med
                  << "  median min: " << minv_med
                  << "  median SNR: " << snr_med << std::endl;
        std::cout << "Reject counts: "
                  << "SNR=" << fail_snr
                  << "  SNR_bottom=" << fail_snr_bottom
                  << "  sat_low=" << fail_sat
                  << "  deriv=" << fail_deriv
                  << "  recover=" << fail_recover
                  << "  sat_high=" << fail_sat_hi
                  << "  pre_close=" << fail_pre_close
                  << "  pre_norec=" << fail_pre_norec
                  << "  post_short=" << fail_post_short
                  << "  le_bump=" << fail_le_bump
                  << "  pre_short=" << fail_pre_short
                  << "  le_dip=" << fail_le_dip
                  << std::endl;
    }

    if (cleaned.empty()) {
        std::cerr << "No single-pulse events passed selection. "
                  << "Consider loosening thresholds." << std::endl;
        return;
    }

    // Step3-1: 粗い揃え（立ち下がり開始位置を合わせる）
    std::vector<std::vector<double>> coarse_aligned;
    coarse_aligned.reserve(cleaned.size());
    for (size_t ie = 0; ie < cleaned.size(); ++ie) {
        const std::vector<double>& y = cleaned[ie];
        double noise_mad = cleaned_mad[ie];
        std::vector<double> ycoarse = y;
        int idx_min = static_cast<int>(std::min_element(y.begin(), y.end()) - y.begin());
        int fall_idx = FindFallingStartIndex(y, noise_mad, idx_min);
        if (fall_idx >= 0) {
            int s_coarse = kCoarseTargetIndex - fall_idx;
            if (IsCorrelationWindowValid(s_coarse, N, kCoarseTargetIndex)) {
                ycoarse = ShiftWave(y, s_coarse);
            }
        }
        coarse_aligned.push_back(ycoarse);
    }

    // Step3-2: 仮テンプレ（粗い揃え後の単純平均）
    std::vector<double> ref(N, 0.0);
    for (auto &y : coarse_aligned)
        for (int i = 0; i < N; ++i)
            ref[i] += y[i];
    for (auto &v : ref) v /= coarse_aligned.size();

    // Step3-3: 相互相関で微調整（leading-edge の短窓のみ）
    int align_win_start = 0;
    int align_win_end = 0;
    BuildAlignWindowFromT0(kCoarseTargetIndex, N, align_win_start, align_win_end);
    for (auto &ycoarse : coarse_aligned) {
        int s = FindShiftByCorrelation(ycoarse, ref, align_win_start, align_win_end, kMaxShift);
        std::vector<double> yshift = ShiftWave(ycoarse, -s);
        aligned.push_back(yshift);

        // イベント毎に正規化（min=-1）してテンプレ用に保存
        std::vector<double> ynorm = yshift;
        double ymin = *std::min_element(ynorm.begin(), ynorm.end());
        if (ymin < 0) {
            for (auto &v : ynorm) v /= -ymin;
        }
        aligned_norm.push_back(ynorm);
    }

    // Step4: 最終テンプレ（イベント毎正規化後に平均）
    std::vector<double> templ(N, 0.0);
    for (auto &y : aligned_norm)
        for (int i = 0; i < N; ++i)
            templ[i] += y[i];
    for (auto &v : templ) v /= aligned_norm.size();

    // min = -1 に正規化（念のため）
    double tmin = *std::min_element(templ.begin(), templ.end());
    if (tmin < 0) {
        for (auto &v : templ) v /= -tmin;
    }

    // ------------------------------
    // 可視化
    // ------------------------------
    TCanvas *c1 = new TCanvas("c1", "Single pulse alignment", 1200, 800);
    c1->Divide(2,1);

    // 左：整列前（選別済み）
    c1->cd(1);
    double y_min_left = 1e99;
    double y_max_left = -1e99;
    for (auto &y : cleaned) {
        auto mm = MinMaxVec(y);
        if (mm.first < y_min_left) y_min_left = mm.first;
        if (mm.second > y_max_left) y_max_left = mm.second;
    }
    auto range_left = ExpandRange(y_min_left, y_max_left);
    bool first_left = true;
    for (auto &y : cleaned) {
        TGraph *g = new TGraph(N);
        for (int i = 0; i < N; ++i) g->SetPoint(i, i, y[i]);
        g->SetLineColorAlpha(kGray+1, 0.3);
        if (first_left) {
            g->SetMinimum(range_left.first);
            g->SetMaximum(range_left.second);
            g->GetXaxis()->SetLimits(0, N - 1);
            g->Draw("AL");
            first_left = false;
        } else {
            g->Draw("L SAME");
        }
    }

    // 右：整列後＋テンプレ（正規化表示）
    c1->cd(2);
    double y_min_right = 1e99;
    double y_max_right = -1e99;
    for (auto &y : aligned_norm) {
        auto mm = MinMaxVec(y);
        if (mm.first < y_min_right) y_min_right = mm.first;
        if (mm.second > y_max_right) y_max_right = mm.second;
        TGraph *g = new TGraph(N);
        for (int i = 0; i < N; ++i) g->SetPoint(i, i, y[i]);
        g->SetLineColorAlpha(kBlue, 0.2);
        g->Draw("L SAME");
    }
    auto mm_t = MinMaxVec(templ);
    if (mm_t.first < y_min_right) y_min_right = mm_t.first;
    if (mm_t.second > y_max_right) y_max_right = mm_t.second;
    auto range_right = ExpandRange(y_min_right, y_max_right);

    TGraph *gtempl = new TGraph(N);
    for (int i = 0; i < N; ++i) gtempl->SetPoint(i, i, templ[i]);
    gtempl->SetLineColor(kRed);
    gtempl->SetLineWidth(3);
    gtempl->SetMinimum(range_right.first);
    gtempl->SetMaximum(range_right.second);
    gtempl->GetXaxis()->SetLimits(0, N - 1);
    gtempl->Draw("AL");

    // 右側の整列後波形を再描画（軸を gtempl に任せる）
    for (auto &y : aligned_norm) {
        TGraph *g = new TGraph(N);
        for (int i = 0; i < N; ++i) g->SetPoint(i, i, y[i]);
        g->SetLineColorAlpha(kBlue, 0.2);
        g->Draw("L SAME");
    }

    // 参考: 除外パルスの表示は省略

    // ------------------------------
    // フィット（c2 は表示しないが、c4 用に実行）
    // ------------------------------
    double thr = kFitFrac * (*std::min_element(templ.begin(), templ.end())); // min は負
    int fit_start = 0;
    int fit_end = N - 1;
    for (int i = 0; i < N; ++i) {
        if (templ[i] < thr) { fit_start = i; break; }
    }
    for (int i = N - 1; i >= 0; --i) {
        if (templ[i] < thr) { fit_end = i; break; }
    }

    // 最終テンプレの t0_index を先に推定し、以後の時刻基準をそろえる
    int t0_slope = FindTemplateT0ByMaxNegSlope(templ);

    TF1 *f = new TF1("f", pulseModel, fit_start, fit_end, 5);
    f->SetParameters(1.0, static_cast<double>(t0_slope), 5.0, 80.0, 2.0);
    f->SetParNames("A","t0","tau_r","tau_d","sigma");

    // フィットは実行（c4 の重ね描きで使用）
    TGraph *gfit = new TGraph(N);
    for (int i = 0; i < N; ++i) gfit->SetPoint(i, i, templ[i]);
    gfit->Fit(f, "R");

    // テンプレート JSON を出力（t0_index は負の最大傾き位置）
    WriteTemplateJson("data/shapeddata/template.json", templ, kSampleDt, kNormDescription,
                      t0_slope, f, run_no, ch_name);

    // ------------------------------
    // 4例＋ランダム4例の可視化（生データ + フィット曲線）
    // ------------------------------
    const int n_examples = std::min(4, static_cast<int>(cleaned.size()));
    if (n_examples > 0) {
        // 先頭4 + ランダム4（重複は許容しない）
        std::vector<int> indices;
        indices.reserve(cleaned.size());
        for (int i = 0; i < static_cast<int>(cleaned.size()); ++i) indices.push_back(i);

        std::random_device rd;
        std::mt19937 gen(rd());
        std::shuffle(indices.begin(), indices.end(), gen);

        std::vector<int> show;
        show.reserve(8);
        for (int i = 0; i < n_examples; ++i) show.push_back(i);
        int added = 0;
        for (int idx : indices) {
            if (added >= 4) break;
            bool exists = false;
            for (int v : show) {
                if (v == idx) { exists = true; break; }
            }
            if (!exists) {
                show.push_back(idx);
                ++added;
            }
        }

        const int n_show = static_cast<int>(show.size());
        TCanvas *c4 = new TCanvas("c4", "Examples: raw + fit", 1400, 900);
        c4->Divide(4, 2);

        for (int k = 0; k < n_show; ++k) {
            int ie = show[k];
            const std::vector<double>& y = cleaned[ie]; // 生データ（整列なし）
            auto it_min = std::min_element(y.begin(), y.end());
            double minv = *it_min;

            // 生データに合わせたモデルシフト量を推定
            double noise_mad = cleaned_mad[ie];
            int idx_min = static_cast<int>(it_min - y.begin());
            int fall_idx = FindFallingStartIndex(y, noise_mad, idx_min);
            int s_coarse = 0;
            if (fall_idx >= 0) {
                s_coarse = kCoarseTargetIndex - fall_idx;
                if (!IsCorrelationWindowValid(s_coarse, N, kCoarseTargetIndex)) {
                    s_coarse = 0;
                }
            }
            std::vector<double> ycoarse = (s_coarse != 0) ? ShiftWave(y, s_coarse) : y;
            int draw_win_start = 0;
            int draw_win_end = 0;
            BuildAlignWindowFromT0(t0_slope, N, draw_win_start, draw_win_end);
            int s_fine = FindShiftByCorrelation(ycoarse, templ, draw_win_start, draw_win_end, kMaxShift);
            int s_total = s_coarse - s_fine;
            double scale = (minv < 0) ? (-minv) : 1.0;

            TGraph *graw = new TGraph(N);
            TGraph *gmodel = new TGraph(N);
            for (int i = 0; i < N; ++i) {
                graw->SetPoint(i, i, y[i]);
                double t = i + s_total; // モデルを raw 座標に合わせて重ねる
                double mv = scale * f->Eval(t);
                gmodel->SetPoint(i, i, mv);
            }

            c4->cd(k + 1);
            graw->SetLineColor(kBlack);
            auto mm_raw = MinMaxVec(y);
            double model_min = 1e99;
            double model_max = -1e99;
            for (int i = 0; i < N; ++i) {
                double t = i + s_total;
                double mv = scale * f->Eval(t);
                if (mv < model_min) model_min = mv;
                if (mv > model_max) model_max = mv;
            }
            double ymin = std::min(mm_raw.first, model_min);
            double ymax = std::max(mm_raw.second, model_max);
            auto range_ex = ExpandRange(ymin, ymax);
            graw->SetMinimum(range_ex.first);
            graw->SetMaximum(range_ex.second);
            graw->GetXaxis()->SetLimits(0, N - 1);
            graw->Draw("AL");
            gmodel->SetLineColor(kRed);
            gmodel->SetLineWidth(2);
            gmodel->Draw("L SAME");
        }
    }
}

// ------------------------------
// ファイルから実行する簡易エントリ
// ------------------------------
void run_singlepulse_from_file(const char* path, int bins_per_event = 1000, int max_events = 0) {
    std::vector<double> all_bins = LoadBins(path);
    std::vector<std::vector<double>> events = SplitEvents(all_bins, bins_per_event);

    if (max_events > 0 && static_cast<size_t>(max_events) < events.size()) {
        events.resize(static_cast<size_t>(max_events));
    }

    int run_no = -1;
    std::string ch_name = "unknown";
    ParseInputMetaFromPath(path, run_no, ch_name);

    std::cout << "Loaded events: " << events.size() << " (bins/event=" << bins_per_event << ")" << std::endl;
    std::cout << "Input meta: run=" << run_no << " CH=" << ch_name << std::endl;
    WriteWaveformsCsv("data/shapeddata/pileup_waveforms.csv", events);
    singlePulseTemplateDemo(events, run_no, ch_name);
}

// ROOT のマクロ呼び出し互換（ファイル名 = 関数名）
void fit_singlepulse(const char* path, int bins_per_event = 1000, int max_events = 0) {
    run_singlepulse_from_file(path, bins_per_event, max_events);
}
