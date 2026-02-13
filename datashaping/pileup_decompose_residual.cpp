#include <fstream>
#include <algorithm>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <iostream>
#include <filesystem>
#include <unordered_map>

// 既定パラメータ
// nSigma: 相関ピーク閾値の係数（sigmaCorr に掛ける）
static const double kDefaultNSigma = 2.0;
// minSep: 候補の最小間隔。負値なら L*ratio で決める
static const double kDefaultMinSepRatio = 0.0; // L/4
// maxK: 保持する候補の最大数（逐次追加の上限）
static const int kDefaultMaxK = 30;
// positiveArea: 面積の符号を正に強制
static const bool kDefaultPositiveArea = false;
// selftest: 内部の簡易テスト波形で実行
static const bool kDefaultSelfTest = false;
// txt_input: 1サンプル/行のテキスト入力を使う
static const bool kDefaultTxtInput = false;
// bins_per_event: txt入力時の1イベントのサンプル数
static const int kDefaultBinsPerEvent = 1000;
// debug: デバッグ出力を有効化
static const bool kDefaultDebug = false;
// no_sign_cut: 極性の符号カットを無効化
static const bool kDefaultNoSignCut = false;
// ridge: 正規方程式に足すリッジ係数（0で無効）
static const double kDefaultRidge = 0.0;
// no_crop: テンプレの切り出しを無効化
static const bool kDefaultNoCrop = false;
// crop_pre/crop_post: t0_index を中心とした切り出し幅
static const int kDefaultCropPre =20;
static const int kDefaultCropPost = 20;
// rss_improve_min: RSS改善率の打ち切り閾値
static const double kDefaultRssImproveMin = 0.005;
// thr0_factor: 候補生成の相関閾値（sigmaCorr に掛ける）
static const double kDefaultThr0Factor = 1.0;
// nms_window: NMSの半窓幅（サンプル数）
static const int kDefaultNmsWindow = 10;
// debug_event: デバッグ出力対象イベント（負値なら全イベント）
static const long long kDefaultDebugEvent = -1;
// 既定パス
static const char* kDefaultTemplatePath = "data/shapeddata/template.json";
static const char* kDefaultBaselinePath = "data/shapeddata/event_baseline_noise.csv";
static const char* kDefaultWavePath = "data/shapeddata/pileup_waveforms.csv";

// ============================================================
// pileup_decompose_residual.cpp
//
// ビルド例:
//   g++ -O2 -std=c++17 -o pileup_decompose datashaping/pileup_decompose.cpp
//
// 入力:
//   1) template.json
//      keys: dt, polarity, template[], areaT, sumT2, t0_index, ...
//   2) event_baseline_noise.csv
//      run,event_id,CH,baseline,noise
//      (旧形式: event_id,baseline,sigma も許容)
//   3) pileup_waveforms.csv
//      1行1イベント: event_id, y0, y1, ..., y(N-1)
//
// 出力:
//   data/shapeddata/pulses.csv
//      run,event_id,pulse_index,area,t_fall_start_sample_raw,amplitude,baseline_fit
//      t_fall_start_sample_raw = (t_start_sample_full + t0_index)
//
// コマンドライン:
//   pileup_decompose [template.json event_baseline_noise.csv pileup_waveforms.csv]
//     [--nSigma 4.0] [--minSep 50] [--maxK 12] [--positiveArea] [--no-sign-cut]
//     [--ridge 1e-6] [--rss-improve-min 0.01]
//     [--crop-pre 10] [--crop-post 20] [--no-crop] [--selftest]
//     [--debug-event 123]
//
//   txt入力(1サンプル/行)の場合:
//   pileup_decompose [template.json event_baseline_noise.csv wave.txt]
//     --txt --bins-per-event 1000
//     [--nSigma 4.0] [--minSep 50] [--maxK 12] [--positiveArea] [--no-sign-cut]
//     [--ridge 1e-6] [--rss-improve-min 0.01]
//     [--crop-pre 10] [--crop-post 20] [--no-crop]
//     [--debug-event 123]
//
// 時刻の定義:
//   - 生データのサンプル番号: i = 0..N-1（入力波形のインデックス）
//   - テンプレの基準点 t0_index: 立ち下がり最大傾き（テンプレ作成側で定義）
//   - フィットは「生データ時刻 i を固定」し、テンプレ開始位置 t_start_sample_full を動かす
//   - 生データ基準の立ち下がり開始サンプル:
//       t_fall_start_sample_raw = t_start_sample_full + t0_index
//
// 注意:
//   - 相関は crop 版テンプレ（短い）で実施
//   - 最終フィットは full テンプレ（長い）で実施
//   - 正規方程式が特異/悪条件ならイベントをスキップ
// ============================================================

struct TemplateData {
    double dt;
    double polarity;
    double areaT;
    double sumT2;
    int t0_index;
    int crop_start;
    int run_no;
    std::vector<double> T;
};

struct BaselineNoise {
    long long event_id;
    double baseline;
    double sigma;
};

struct Pulse {
    int t_start;
    double amplitude;
};

static bool ReadFileToString(const std::string& path, std::string& out) {
    std::ifstream ifs(path.c_str());
    if (!ifs) return false;
    std::ostringstream ss;
    ss << ifs.rdbuf();
    out = ss.str();
    return true;
}

static bool ExtractNumberAfterKey(const std::string& s, const std::string& key, double& out) {
    std::string pat = "\"" + key + "\"";
    size_t p = s.find(pat);
    if (p == std::string::npos) return false;
    p = s.find(':', p);
    if (p == std::string::npos) return false;
    ++p;
    while (p < s.size() && (s[p] == ' ' || s[p] == '\n' || s[p] == '\r' || s[p] == '\t')) ++p;
    size_t end = p;
    while (end < s.size() && (s[end] == '-' || s[end] == '+' || s[end] == '.' ||
                              (s[end] >= '0' && s[end] <= '9') || s[end] == 'e' || s[end] == 'E')) {
        ++end;
    }
    if (end == p) return false;
    out = std::atof(s.substr(p, end - p).c_str());
    return true;
}

static bool ExtractStringAfterKey(const std::string& s, const std::string& key, std::string& out) {
    std::string pat = "\"" + key + "\"";
    size_t p = s.find(pat);
    if (p == std::string::npos) return false;
    p = s.find(':', p);
    if (p == std::string::npos) return false;
    ++p;
    while (p < s.size() && (s[p] == ' ' || s[p] == '\n' || s[p] == '\r' || s[p] == '\t')) ++p;
    if (p >= s.size() || s[p] != '\"') return false;
    ++p;
    size_t end = s.find('\"', p);
    if (end == std::string::npos) return false;
    out = s.substr(p, end - p);
    return true;
}

static bool ParseLongLong(const std::string& s, long long& out) {
    char* endp = nullptr;
    long long v = std::strtoll(s.c_str(), &endp, 10);
    if (endp == s.c_str() || *endp != '\0') return false;
    out = v;
    return true;
}

static bool ParseDouble(const std::string& s, double& out) {
    char* endp = nullptr;
    double v = std::strtod(s.c_str(), &endp);
    if (endp == s.c_str() || *endp != '\0') return false;
    out = v;
    return true;
}

static bool ExtractIntAfterKey(const std::string& s, const std::string& key, int& out) {
    double v = 0.0;
    if (!ExtractNumberAfterKey(s, key, v)) return false;
    out = static_cast<int>(std::floor(v + 0.5));
    return true;
}

static bool ExtractArrayAfterKey(const std::string& s, const std::string& key, std::vector<double>& out) {
    std::string pat = "\"" + key + "\"";
    size_t p = s.find(pat);
    if (p == std::string::npos) return false;
    p = s.find('[', p);
    if (p == std::string::npos) return false;
    size_t q = s.find(']', p);
    if (q == std::string::npos) return false;
    std::string body = s.substr(p + 1, q - p - 1);
    std::stringstream ss(body);
    std::string token;
    out.clear();
    while (std::getline(ss, token, ',')) {
        std::stringstream ts(token);
        double v = 0.0;
        ts >> v;
        out.push_back(v);
    }
    return !out.empty();
}

static bool LoadTemplate(const std::string& path, TemplateData& t) {
    std::string s;
    if (!ReadFileToString(path, s)) return false;
    if (!ExtractNumberAfterKey(s, "dt", t.dt)) return false;
    if (!ExtractNumberAfterKey(s, "polarity", t.polarity)) {
        std::string pol;
        if (!ExtractStringAfterKey(s, "polarity", pol)) return false;
        if (pol == "negative") t.polarity = -1.0;
        else if (pol == "positive") t.polarity = 1.0;
        else return false;
    }
    if (!ExtractNumberAfterKey(s, "areaT", t.areaT)) return false;
    if (!ExtractNumberAfterKey(s, "sumT2", t.sumT2)) return false;
    if (!ExtractIntAfterKey(s, "t0_index", t.t0_index)) return false;
    if (!ExtractArrayAfterKey(s, "template", t.T)) return false;
    t.crop_start = 0;
    t.run_no = -1;
    ExtractIntAfterKey(s, "run", t.run_no);
    return true;
}

static void CropTemplateAroundT0(TemplateData& t, int pre, int post) {
    if (t.T.empty()) return;
    const int N = static_cast<int>(t.T.size());
    int start = t.t0_index - pre;
    if (start < 0) start = 0;
    int end = t.t0_index + post;
    if (end > N - 1) end = N - 1;
    if (start == 0 && end == N - 1) return;

    std::vector<double> out;
    out.reserve(end - start + 1);
    for (int i = start; i <= end; ++i) out.push_back(t.T[i]);
    t.crop_start = start;
    t.t0_index = t.t0_index - start;
    t.T.swap(out);

    // areaT / sumT2 は切り出し後に再計算
    double sum = 0.0;
    double sum2 = 0.0;
    for (double v : t.T) {
        sum += v;
        sum2 += v * v;
    }
    t.areaT = sum * t.dt;
    t.sumT2 = sum2;
}

static bool LoadBaselineNoise(const std::string& path, std::vector<BaselineNoise>& out) {
    std::ifstream ifs(path.c_str());
    if (!ifs) return false;
    out.clear();
    std::string line;
    while (std::getline(ifs, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string col;
        BaselineNoise bn;
        if (!std::getline(ss, col, ',')) continue;
        if (!ParseLongLong(col, bn.event_id)) {
            // 新形式: run,event_id,CH,baseline,noise
            // 旧形式: event_id,baseline,sigma
            // ヘッダ or 非数値行はスキップ
            // new format: 先頭は run なので event_id は次の列
            if (!std::getline(ss, col, ',')) continue;
            if (!ParseLongLong(col, bn.event_id)) continue;
            // CH
            if (!std::getline(ss, col, ',')) continue;
            // baseline
            if (!std::getline(ss, col, ',')) continue;
            if (!ParseDouble(col, bn.baseline)) continue;
            // noise
            if (!std::getline(ss, col, ',')) continue;
            if (!ParseDouble(col, bn.sigma)) continue;
            out.push_back(bn);
            continue;
        }
        // 旧形式: event_id,baseline,sigma
        if (!std::getline(ss, col, ',')) continue;
        if (!ParseDouble(col, bn.baseline)) continue;
        if (!std::getline(ss, col, ',')) continue;
        if (!ParseDouble(col, bn.sigma)) continue;
        out.push_back(bn);
    }
    return true;
}

static void BuildBaselineMap(const std::vector<BaselineNoise>& v,
                             std::unordered_map<long long, BaselineNoise>& out) {
    out.clear();
    out.reserve(v.size());
    for (size_t i = 0; i < v.size(); ++i) {
        out[v[i].event_id] = v[i];
    }
}

static bool FindBaselineNoise(const std::unordered_map<long long, BaselineNoise>& m,
                              long long event_id, double& baseline, double& sigma) {
    auto it = m.find(event_id);
    if (it == m.end()) return false;
    baseline = it->second.baseline;
    sigma = it->second.sigma;
    return true;
}

static void EstimateBaselineSigma(const std::vector<double>& y, double& baseline, double& sigma) {
    const int n0 = (y.size() < 200) ? static_cast<int>(y.size()) : 200;
    if (n0 <= 1) { baseline = 0.0; sigma = 1.0; return; }
    std::vector<double> tmp;
    tmp.reserve(n0);
    for (int i = 0; i < n0; ++i) tmp.push_back(y[i]);
    std::sort(tmp.begin(), tmp.end());
    double med = tmp[n0 / 2];
    // MAD
    for (int i = 0; i < n0; ++i) tmp[i] = std::fabs(tmp[i] - med);
    std::sort(tmp.begin(), tmp.end());
    double mad = tmp[n0 / 2];
    baseline = med;
    sigma = 1.4826 * mad;
    if (sigma <= 0.0) sigma = 1.0;
}

static void ComputeCorrelation(const std::vector<double>& x, const std::vector<double>& T,
                               int offset, std::vector<double>& c, std::vector<double>& a_hat, double sumT2) {
    const int N = static_cast<int>(x.size());
    const int L = static_cast<int>(T.size());
    const int S = N - (offset + L) + 1;
    if (S <= 0) {
        c.clear();
        a_hat.clear();
        return;
    }
    c.assign(S, 0.0);
    a_hat.assign(S, 0.0);
    for (int s = 0; s < S; ++s) {
        double acc = 0.0;
        for (int j = 0; j < L; ++j) acc += x[s + offset + j] * T[j];
        c[s] = acc;
        a_hat[s] = acc / sumT2;
    }
}

static std::vector<int> SelectCandidates(const std::vector<double>& c, double sigmaCorr,
                                         double thr0_factor, int nms_window, int maxK,
                                         int dom_sign) {
    const int S = static_cast<int>(c.size());
    std::vector<int> local;
    const double thr0 = thr0_factor * sigmaCorr;
    if (S == 1) {
        if (std::fabs(c[0]) > thr0 && (dom_sign > 0 ? (c[0] > 0.0) : (c[0] < 0.0))) {
            local.push_back(0);
        }
    } else if (S == 2) {
        if (std::fabs(c[0]) >= std::fabs(c[1]) && std::fabs(c[0]) > thr0 &&
            (dom_sign > 0 ? (c[0] > 0.0) : (c[0] < 0.0))) {
            local.push_back(0);
        }
        if (std::fabs(c[1]) > std::fabs(c[0]) && std::fabs(c[1]) > thr0 &&
            (dom_sign > 0 ? (c[1] > 0.0) : (c[1] < 0.0))) {
            local.push_back(1);
        }
    } else if (S >= 3) {
        for (int s = 0; s < S; ++s) {
            if (std::fabs(c[s]) <= thr0) continue;
            if (dom_sign > 0) {
                if (c[s] <= 0.0) continue;
            } else {
                if (c[s] >= 0.0) continue;
            }
            int lo = s - nms_window;
            int hi = s + nms_window;
            if (lo < 0) lo = 0;
            if (hi > S - 1) hi = S - 1;
            bool is_max = true;
            for (int j = lo; j <= hi; ++j) {
                if (std::fabs(c[j]) > std::fabs(c[s])) { is_max = false; break; }
            }
            if (is_max) local.push_back(s);
        }
    }
    // sort by |c| descending (simple selection)
    std::vector<int> selected;
    while (!local.empty() && static_cast<int>(selected.size()) < maxK) {
        int best_idx = 0;
        double best_val = std::fabs(c[local[0]]);
        for (int i = 1; i < static_cast<int>(local.size()); ++i) {
            if (std::fabs(c[local[i]]) > best_val) {
                best_val = std::fabs(c[local[i]]);
                best_idx = i;
            }
        }
        int s = local[best_idx];
        local.erase(local.begin() + best_idx);
        selected.push_back(s);
    }
    return selected;
}

static bool SolveNormalEq(const std::vector< std::vector<double> >& A,
                          const std::vector<double>& b,
                          std::vector<double>& x,
                          double ridge) {
    const int n = static_cast<int>(b.size());
    std::vector< std::vector<double> > M = A;
    std::vector<double> v = b;
    x.assign(n, 0.0);
    if (ridge > 0.0 && n > 0) {
        double scale = std::fabs(M[0][0]);
        if (scale <= 0.0) scale = 1.0;
        double lambda = ridge * scale;
        for (int i = 0; i < n; ++i) M[i][i] += lambda;
    }
    for (int i = 0; i < n; ++i) {
        // pivot
        int piv = i;
        double amax = std::fabs(M[i][i]);
        for (int r = i + 1; r < n; ++r) {
            double av = std::fabs(M[r][i]);
            if (av > amax) { amax = av; piv = r; }
        }
        if (amax < 1e-12) return false;
        if (piv != i) {
            std::vector<double> tmp = M[i];
            M[i] = M[piv];
            M[piv] = tmp;
            double tv = v[i];
            v[i] = v[piv];
            v[piv] = tv;
        }
        double diag = M[i][i];
        for (int c = i; c < n; ++c) M[i][c] /= diag;
        v[i] /= diag;
        for (int r = 0; r < n; ++r) {
            if (r == i) continue;
            double f = M[r][i];
            if (f == 0.0) continue;
            for (int c = i; c < n; ++c) M[r][c] -= f * M[i][c];
            v[r] -= f * v[i];
        }
    }
    for (int i = 0; i < n; ++i) x[i] = v[i];
    return true;
}

static bool FitAmplitudes(const std::vector<double>& x, const std::vector<double>& T,
                          const std::vector<int>& tks, double& b_out, std::vector<double>& a_out,
                          double ridge) {
    const int N = static_cast<int>(x.size());
    const int L = static_cast<int>(T.size());
    const int M = static_cast<int>(tks.size());
    const int P = 1 + M;
    std::vector< std::vector<double> > ATA(P, std::vector<double>(P, 0.0));
    std::vector<double> ATy(P, 0.0);

    for (int i = 0; i < N; ++i) {
        double u0 = 1.0;
        ATy[0] += u0 * x[i];
        ATA[0][0] += u0 * u0;
        for (int k = 0; k < M; ++k) {
            int j = i - tks[k];
            double uk = (j >= 0 && j < L) ? T[j] : 0.0;
            ATy[k + 1] += uk * x[i];
            ATA[0][k + 1] += u0 * uk;
            ATA[k + 1][0] += u0 * uk;
            for (int l = k; l < M; ++l) {
                int j2 = i - tks[l];
                double ul = (j2 >= 0 && j2 < L) ? T[j2] : 0.0;
                ATA[k + 1][l + 1] += uk * ul;
                if (l != k) ATA[l + 1][k + 1] = ATA[k + 1][l + 1];
            }
        }
    }

    std::vector<double> sol;
    if (!SolveNormalEq(ATA, ATy, sol, ridge)) return false;
    b_out = sol[0];
    a_out.assign(M, 0.0);
    for (int k = 0; k < M; ++k) a_out[k] = sol[k + 1];
    return true;
}

static double ComputeRSS(const std::vector<double>& x, const std::vector<double>& T,
                         const std::vector<int>& tks, double b_fit, const std::vector<double>& a_fit) {
    const int N = static_cast<int>(x.size());
    const int L = static_cast<int>(T.size());
    const int M = static_cast<int>(tks.size());
    double rss = 0.0;
    for (int i = 0; i < N; ++i) {
        double yhat = b_fit;
        for (int k = 0; k < M; ++k) {
            int j = i - tks[k];
            if (j < 0 || j >= L) continue;
            yhat += a_fit[k] * T[j];
        }
        double r = x[i] - yhat;
        rss += r * r;
    }
    return rss;
}

static void ComputeResidual(const std::vector<double>& x, const std::vector<double>& T,
                            const std::vector<int>& tks, double b_fit, const std::vector<double>& a_fit,
                            std::vector<double>& r) {
    const int N = static_cast<int>(x.size());
    const int L = static_cast<int>(T.size());
    const int M = static_cast<int>(tks.size());
    r.assign(N, 0.0);
    for (int i = 0; i < N; ++i) {
        double yhat = b_fit;
        for (int k = 0; k < M; ++k) {
            int j = i - tks[k];
            if (j < 0 || j >= L) continue;
            yhat += a_fit[k] * T[j];
        }
        r[i] = x[i] - yhat;
    }
}

static bool FindFall50(const std::vector<double>& T, int t_start, double a,
                       double polarity, double& fall_sample, double& fall_time,
                       double dt, int t0_index) {
    const int L = static_cast<int>(T.size());
    // ピーク探索（polarityの符号で極大/極小を選択）
    int peak_j = 0;
    double peak = a * T[0];
    for (int j = 1; j < L; ++j) {
        double v = a * T[j];
        if (polarity < 0) {
            if (v < peak) { peak = v; peak_j = j; }
        } else {
            if (v > peak) { peak = v; peak_j = j; }
        }
    }
    double target = 0.5 * std::fabs(peak);
    int j_found = -1;
    for (int j = peak_j + 1; j < L; ++j) {
        double v = std::fabs(a * T[j]);
        if (v <= target) { j_found = j; break; }
    }
    if (j_found < 0) return false;
    int j0 = j_found - 1;
    double v0 = std::fabs(a * T[j0]);
    double v1 = std::fabs(a * T[j_found]);
    double frac = 0.0;
    if (v0 != v1) frac = (target - v0) / (v1 - v0);
    double sample = (t_start + j0) + frac;
    fall_sample = sample;
    fall_time = (sample - t0_index) * dt;
    return true;
}

static void WritePulsesHeader(std::ofstream& ofs) {
    ofs << "run,event_id,pulse_index,area,t_fall_start_sample_raw,amplitude,baseline_fit\n";
}

static void ProcessEvent(long long event_id, const std::vector<double>& y,
                         const TemplateData& t_corr, const TemplateData& t_full,
                         const std::unordered_map<long long, BaselineNoise>& bnmap,
                         double nSigma, int minSep, int maxK, bool positiveArea,
                         bool no_sign_cut,
                         bool debug,
                         double ridge,
                         double rss_improve_min,
                         long long debug_event,
                         std::ofstream& ofs) {
    const bool dbg = debug && (debug_event < 0 || event_id == debug_event);
    if (y.size() < t_corr.T.size()) return;
    double baseline = 0.0, sigma = 0.0;
    if (!FindBaselineNoise(bnmap, event_id, baseline, sigma)) {
        EstimateBaselineSigma(y, baseline, sigma);
        if (dbg) {
            std::cerr << "[event " << event_id << "] baseline/sigma estimated: "
                      << baseline << ", " << sigma << "\n";
        }
    }
    std::vector<double> x(y.size(), 0.0);
    for (size_t i = 0; i < y.size(); ++i) x[i] = y[i] - baseline;

    double sigmaCorr = sigma * std::sqrt(t_corr.sumT2);
    const double sigmaA = (t_full.sumT2 > 0.0) ? (sigma / std::sqrt(t_full.sumT2)) : sigma;
    const double aMin = 3.0 * sigmaA;
    const double thr = nSigma * sigmaCorr;

    std::vector<int> tks_full;
    std::vector<int> tks_crop;
    std::vector<double> c_at_accept;
    tks_full.reserve(maxK);
    tks_crop.reserve(maxK);
    c_at_accept.reserve(maxK);
    double b_fit = 0.0;
    std::vector<double> a_fit;
    double rss_prev = -1.0;

    std::vector<double> residual = x;
    std::vector<double> c_res, a_hat_res;

    for (int iter = 0; iter < maxK; ++iter) {
        ComputeCorrelation(residual, t_corr.T, 0, c_res, a_hat_res, t_corr.sumT2);
        if (c_res.empty()) break;

        double cmax = c_res[0];
        double cmin = c_res[0];
        for (size_t i = 1; i < c_res.size(); ++i) {
            if (c_res[i] > cmax) cmax = c_res[i];
            if (c_res[i] < cmin) cmin = c_res[i];
        }
        int dom_sign = (std::fabs(cmax) >= std::fabs(cmin)) ? 1 : -1;
        std::vector<int> cand = SelectCandidates(c_res, sigmaCorr, kDefaultThr0Factor,
                                                 kDefaultNmsWindow, maxK, dom_sign);
        if (dbg) {
            std::cerr << "[event " << event_id << "] iter=" << iter
                      << " candidates=" << cand.size()
                      << " cmax=" << cmax << " cmin=" << cmin
                      << " sigmaCorr=" << sigmaCorr
                      << " thr=" << (nSigma * sigmaCorr)
                      << " dom_sign=" << (dom_sign > 0 ? "+" : "-") << "\n";
            if (!cand.empty()) {
                std::cerr << "[event " << event_id << "] ranked_candidates:";
                for (size_t i = 0; i < cand.size(); ++i) {
                    int s = cand[i];
                    std::cerr << " " << s << "(" << c_res[s] << ")";
                }
                std::cerr << "\n";
            }
        }
        if (cand.empty()) break;

        bool added = false;
        for (size_t idx = 0; idx < cand.size(); ++idx) {
            int s = cand[idx];
            int s_full = s - t_corr.crop_start;
            bool ok_sep = true;
            for (size_t k = 0; k < tks_full.size(); ++k) {
                if (std::abs(tks_full[k] - s_full) < minSep) { ok_sep = false; break; }
            }
            if (!ok_sep) continue;

            std::vector<int> tks_full_try = tks_full;
            tks_full_try.push_back(s_full);
            std::vector<double> a_try;
            double b_try = 0.0;
            if (!FitAmplitudes(x, t_full.T, tks_full_try, b_try, a_try, ridge)) {
                if (dbg) {
                    std::cerr << "[event " << event_id << "] fit failed (singular/ill-conditioned)\n";
                }
                continue;
            }
            // 逐次追加の採否判定：a*c>0 を満たさない候補はスキップ
            if (!no_sign_cut && !a_try.empty()) {
                if (a_try.back() * c_res[s] <= 0.0) {
                    continue;
                }
            }
            // 1発目以降のみ：追加した成分が弱すぎるならスキップ
            if (tks_full_try.size() > 1 && !a_try.empty() && std::fabs(a_try.back()) < aMin) {
                continue;
            }
            // 2発目以降のみ：相関閾値より小さい候補ならスキップ
            if (tks_full_try.size() > 1 && std::fabs(c_res[s]) < thr) {
                continue;
            }

            double rss = ComputeRSS(x, t_full.T, tks_full_try, b_try, a_try);
            if (rss_prev > 0.0) {
                double improve = (rss_prev - rss) / rss_prev;
                if (improve < rss_improve_min) {
                    continue;
                }
            }

            tks_full.swap(tks_full_try);
            tks_crop.push_back(s);
            c_at_accept.push_back(c_res[s]);
            a_fit.swap(a_try);
            b_fit = b_try;
            rss_prev = rss;
            added = true;
            if (dbg) {
                std::cerr << "[event " << event_id << "] accept t_start=" << s
                          << " t_start_full=" << s_full
                          << " a=" << (a_fit.empty() ? 0.0 : a_fit.back()) << "\n";
            }
            break;
        }
        if (!added) break;

        ComputeResidual(x, t_full.T, tks_full, b_fit, a_fit, residual);
    }
    if (tks_full.empty()) return;

    // b_fit, a_fit は逐次追加で決定済み

    const double baseline_fit = baseline + b_fit;
    int pulse_index = 0;
    for (size_t k = 0; k < tks_full.size(); ++k) {
        double a = a_fit[k];
        bool keep = true;
        if (!no_sign_cut) {
            if (k < c_at_accept.size()) {
                keep = (a * c_at_accept[k] > 0.0);
            } else {
                keep = false;
            }
        }
        if (!keep) {
            if (dbg) {
                std::cerr << "[event " << event_id << "] reject k=" << k
                          << " (sign) a=" << a << "\n";
            }
            continue;
        }
        if (std::fabs(a) < aMin) {
            if (dbg) {
                std::cerr << "[event " << event_id << "] reject k=" << k
                          << " (|a|<aMin) a=" << a << " aMin=" << aMin << "\n";
            }
            continue;
        }

        int t_start = tks_crop[k];
        double t_start_full = static_cast<double>(tks_full[k]);
        double t_start_time = (t_start_full - t_full.t0_index) * t_full.dt;
        double area = a * t_full.areaT;
        if (positiveArea && area < 0.0) area = -area;

        double fall_sample = -1.0;
        double fall_time = 0.0;
        bool ok_fall = FindFall50(t_full.T, static_cast<int>(t_start_full), a, t_full.polarity, fall_sample, fall_time,
                                  t_full.dt, t_full.t0_index);
        if (!ok_fall) {
            fall_sample = static_cast<double>(t_start_full + static_cast<int>(t_full.T.size()) - 1);
            fall_time = (fall_sample - t_full.t0_index) * t_full.dt;
        }

        if (dbg) {
            double cv = (k < c_at_accept.size()) ? c_at_accept[k] : 0.0;
            std::cerr << "WRITE t=" << t_start << " t_full=" << t_start_full
                      << " a=" << a << " c=" << cv << "\n";
        }
        double fall_start_sample_raw =
            (t_start_full + t_full.t0_index);
        ofs << t_full.run_no << "," << event_id << "," << pulse_index << "," << area
            << "," << fall_start_sample_raw << "," << a << "," << baseline_fit << "\n";
        ++pulse_index;
    }
    if (dbg) {
        std::cerr << "[event " << event_id << "] pulses_written=" << pulse_index << "\n";
    }
}

static void ParseWaveformLine(const std::string& line, long long& event_id, std::vector<double>& y) {
    y.clear();
    std::stringstream ss(line);
    std::string col;
    if (!std::getline(ss, col, ',')) { event_id = -1; return; }
    event_id = std::atoll(col.c_str());
    while (std::getline(ss, col, ',')) {
        y.push_back(std::atof(col.c_str()));
    }
}

static void RunSelfTest(const TemplateData& t_corr, const TemplateData& t_full,
                        const std::unordered_map<long long, BaselineNoise>& bnmap,
                        double nSigma, int minSep, int maxK, bool positiveArea,
                        bool no_sign_cut,
                        bool debug,
                        double ridge,
                        double rss_improve_min,
                        long long debug_event,
                        std::ofstream& ofs) {
    const int L = static_cast<int>(t_full.T.size());
    const int N = L * 4;
    std::vector<double> y(N, 100.0);
    int t1 = L / 2;
    int t2 = L / 2 + L / 3;
    double a1 = 1.2;
    double a2 = 0.8;
    for (int i = 0; i < N; ++i) {
        int j1 = i - t1;
        int j2 = i - t2;
        if (j1 >= 0 && j1 < L) y[i] += a1 * t_full.T[j1];
        if (j2 >= 0 && j2 < L) y[i] += a2 * t_full.T[j2];
        // 擬似ノイズ（決定論的）
        double noise = 0.5 * std::sin(0.123 * i) + 0.3 * std::sin(0.077 * i);
        y[i] += noise;
    }
    ProcessEvent(0, y, t_corr, t_full, bnmap, nSigma, minSep, maxK,
                 positiveArea, no_sign_cut, debug, ridge, rss_improve_min, debug_event, ofs);
}

static bool ReadTxtWaveforms(const std::string& path, int bins_per_event,
                             std::vector<long long>& event_ids,
                             std::vector< std::vector<double> >& waves) {
    std::ifstream ifs(path.c_str());
    if (!ifs) return false;
    event_ids.clear();
    waves.clear();
    std::vector<double> cur;
    cur.reserve(bins_per_event);
    std::string line;
    long long event_id = 0;
    while (std::getline(ifs, line)) {
        if (line.empty()) continue;
        double v = std::atof(line.c_str());
        cur.push_back(v);
        if (static_cast<int>(cur.size()) == bins_per_event) {
            event_ids.push_back(event_id);
            waves.push_back(cur);
            cur.clear();
            ++event_id;
        }
    }
    if (!cur.empty()) {
        std::cerr << "warn: trailing samples dropped ("
                  << cur.size() << " < bins_per_event=" << bins_per_event << ")\n";
    }
    return !waves.empty();
}

int main(int argc, char** argv) {
    std::string template_path = kDefaultTemplatePath;
    std::string baseline_path = kDefaultBaselinePath;
    std::string wave_path = kDefaultWavePath;
    double nSigma = kDefaultNSigma;
    int minSep = -1; // 負値なら L*ratio
    int maxK = kDefaultMaxK;
    bool positiveArea = kDefaultPositiveArea;
    bool selftest = kDefaultSelfTest;
    bool txt_input = kDefaultTxtInput;
    int bins_per_event = kDefaultBinsPerEvent;
    bool debug = kDefaultDebug;
    bool no_sign_cut = kDefaultNoSignCut;
    double ridge = kDefaultRidge;
    double rss_improve_min = kDefaultRssImproveMin;
    bool no_crop = kDefaultNoCrop;
    int crop_pre = kDefaultCropPre;
    int crop_post = kDefaultCropPost;
    long long debug_event = kDefaultDebugEvent;
    std::vector<std::string> positional;
    positional.reserve(argc);
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--nSigma" && i + 1 < argc) {
            nSigma = std::atof(argv[++i]);
        } else if (a == "--minSep" && i + 1 < argc) {
            minSep = std::atoi(argv[++i]);
        } else if (a == "--maxK" && i + 1 < argc) {
            maxK = std::atoi(argv[++i]);
        } else if (a == "--positiveArea") {
            positiveArea = true;
        } else if (a == "--selftest") {
            selftest = true;
        } else if (a == "--txt") {
            txt_input = true;
        } else if (a == "--bins-per-event" && i + 1 < argc) {
            bins_per_event = std::atoi(argv[++i]);
        } else if (a == "--debug") {
            debug = true;
        } else if (a == "--debug-event" && i + 1 < argc) {
            debug_event = std::atoll(argv[++i]);
        } else if (a == "--no-sign-cut") {
            no_sign_cut = true;
        } else if (a == "--ridge" && i + 1 < argc) {
            ridge = std::atof(argv[++i]);
        } else if (a == "--rss-improve-min" && i + 1 < argc) {
            rss_improve_min = std::atof(argv[++i]);
        } else if (a == "--no-crop") {
            no_crop = true;
        } else if (a == "--crop-pre" && i + 1 < argc) {
            crop_pre = std::atoi(argv[++i]);
        } else if (a == "--crop-post" && i + 1 < argc) {
            crop_post = std::atoi(argv[++i]);
        } else if (a.rfind("--", 0) == 0) {
            // unknown option -> ignore
        } else {
            positional.push_back(a);
        }
    }
    if (positional.size() >= 1) template_path = positional[0];
    if (positional.size() >= 2) baseline_path = positional[1];
    if (positional.size() >= 3) wave_path = positional[2];

    TemplateData t_full;
    if (!LoadTemplate(template_path, t_full)) {
        std::cerr << "error: failed to load template: " << template_path << "\n";
        return 1;
    }
    TemplateData t_corr = t_full;
    if (!no_crop) {
        CropTemplateAroundT0(t_corr, crop_pre, crop_post);
    }
    if (minSep < 0) {
        minSep = static_cast<int>(std::floor(t_corr.T.size() * kDefaultMinSepRatio + 0.5));
        if (minSep < 1) minSep = 1;
    }

    std::vector<BaselineNoise> bnlist;
    if (!LoadBaselineNoise(baseline_path, bnlist)) {
        std::cerr << "warn: failed to load baseline list: " << baseline_path << "\n";
    }
    std::unordered_map<long long, BaselineNoise> bnmap;
    BuildBaselineMap(bnlist, bnmap);

    std::error_code ec;
    std::filesystem::create_directories("data/shapeddata", ec);
    if (ec) {
        std::cerr << "error: failed to create output dir data/shapeddata: " << ec.message() << "\n";
        return 1;
    }
    std::ofstream ofs("data/shapeddata/pulses.csv");
    if (!ofs) {
        std::cerr << "error: failed to open output: data/shapeddata/pulses.csv\n";
        return 1;
    }
    WritePulsesHeader(ofs);

    if (selftest) {
        RunSelfTest(t_corr, t_full, bnmap, nSigma, minSep, maxK, positiveArea, no_sign_cut,
                    debug, ridge, rss_improve_min, debug_event, ofs);
        return 0;
    }

    if (txt_input) {
        std::vector<long long> event_ids;
        std::vector< std::vector<double> > waves;
        if (!ReadTxtWaveforms(wave_path, bins_per_event, event_ids, waves)) {
            std::cerr << "error: failed to read txt waveform: " << wave_path << "\n";
            return 1;
        }
        if (debug) {
            std::cerr << "info: txt events=" << waves.size()
                      << " bins_per_event=" << bins_per_event << "\n";
        }
        for (size_t i = 0; i < waves.size(); ++i) {
            ProcessEvent(event_ids[i], waves[i], t_corr, t_full, bnmap, nSigma, minSep, maxK,
                         positiveArea, no_sign_cut, debug, ridge, rss_improve_min, debug_event, ofs);
        }
    } else {
        std::ifstream ifs(wave_path.c_str());
        if (!ifs) {
            std::cerr << "error: failed to open waveform csv: " << wave_path << "\n";
            return 1;
        }
        std::string line;
        while (std::getline(ifs, line)) {
            if (line.empty()) continue;
            long long event_id = -1;
            std::vector<double> y;
            ParseWaveformLine(line, event_id, y);
            if (event_id < 0 || y.empty()) continue;
            ProcessEvent(event_id, y, t_corr, t_full, bnmap, nSigma, minSep, maxK,
                         positiveArea, no_sign_cut, debug, ridge, rss_improve_min, debug_event, ofs);
        }
    }

    return 0;
}
