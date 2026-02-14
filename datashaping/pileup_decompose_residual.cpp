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
// minSep: 候補の最小間隔 [sample]
static const int kDefaultMinSepSamples = 20;
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
static const double kDefaultRidge = 0;
// no_crop: テンプレの切り出しを無効化
static const bool kDefaultNoCrop = false;
// crop_pre/crop_post: t0_index を中心とした切り出し幅
static const int kDefaultCropPre =20;
static const int kDefaultCropPost = 20;
// rss_improve_min: RSS改善率の打ち切り閾値
static const double kDefaultRssImproveMin = 0.005;
// aMin_sigma_factor: 振幅下限 aMin を決める係数（aMin = factor * sigmaA）
static const double kDefaultAMinSigmaFactor = 3.0;
// thr0_factor: 候補生成の相関閾値（sigmaCorr に掛ける）
static const double kDefaultThr0Factor = 2.0;
// nms_window: NMSの半窓幅（サンプル数）
static const int kDefaultNmsWindow = 10;
// local_refine_half_width: 候補採用後に再探索する時刻の半幅（±この値）
static const int kDefaultLocalRefineHalfWidth = 2;
// post_accept_refine_half_width: 採用後に全パルス時刻を再調整する半幅（±この値）
static const int kDefaultPostAcceptRefineHalfWidth = 10;
// iter0_rank_first: iter=0 は相関順位を優先し、最初に通る候補を採用する
static const bool kDefaultIter0RankFirst = true;
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
        BaselineNoise bn;
        std::vector<std::string> cols;
        {
            std::stringstream ss(line);
            std::string col;
            while (std::getline(ss, col, ',')) cols.push_back(col);
        }
        if (cols.size() >= 5) {
            // 新形式: run,event_id,CH,baseline,noise
            if (!ParseLongLong(cols[1], bn.event_id)) continue;
            if (!ParseDouble(cols[3], bn.baseline)) continue;
            if (!ParseDouble(cols[4], bn.sigma)) continue;
            out.push_back(bn);
            continue;
        }
        if (cols.size() >= 3) {
            // 旧形式: event_id,baseline,sigma
            if (!ParseLongLong(cols[0], bn.event_id)) continue;
            if (!ParseDouble(cols[1], bn.baseline)) continue;
            if (!ParseDouble(cols[2], bn.sigma)) continue;
            out.push_back(bn);
            continue;
        }
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

// 相関値 c[s] = sum_j x[s+j] T[j] を単一点で評価
static bool ComputeCorrelationAtStart(const std::vector<double>& x, const std::vector<double>& T,
                                      int s, double& c_out) {
    const int N = static_cast<int>(x.size());
    const int L = static_cast<int>(T.size());
    if (s < 0 || s + L > N) return false;
    double acc = 0.0;
    for (int j = 0; j < L; ++j) acc += x[s + j] * T[j];
    c_out = acc;
    return true;
}

static std::vector<int> SelectCandidates(const std::vector<double>& c, double sigmaCorr,
                                         double thr0_factor, int nms_window, int maxK) {
    const int S = static_cast<int>(c.size());
    std::vector<int> local;
    const double thr0 = thr0_factor * sigmaCorr;
    if (S == 1) {
        if (std::fabs(c[0]) > thr0 && c[0] > 0.0) {
            local.push_back(0);
        }
    } else if (S == 2) {
        if (std::fabs(c[0]) >= std::fabs(c[1]) && std::fabs(c[0]) > thr0 &&
            c[0] > 0.0) {
            local.push_back(0);
        }
        if (std::fabs(c[1]) > std::fabs(c[0]) && std::fabs(c[1]) > thr0 &&
            c[1] > 0.0) {
            local.push_back(1);
        }
    } else if (S >= 3) {
        for (int s = 0; s < S; ++s) {
            if (std::fabs(c[s]) <= thr0) continue;
            if (c[s] <= 0.0) continue;
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

// ベースライン b_fixed を固定して、振幅 a だけ最小二乗で解く
static bool FitAmplitudesFixedBaseline(const std::vector<double>& x, const std::vector<double>& T,
                                       const std::vector<int>& tks, double b_fixed,
                                       std::vector<double>& a_out,
                                       double ridge) {
    const int N = static_cast<int>(x.size());
    const int L = static_cast<int>(T.size());
    const int M = static_cast<int>(tks.size());
    if (M <= 0) {
        a_out.clear();
        return true;
    }

    std::vector< std::vector<double> > ATA(M, std::vector<double>(M, 0.0));
    std::vector<double> ATy(M, 0.0);

    for (int i = 0; i < N; ++i) {
        const double yi = x[i] - b_fixed;
        for (int k = 0; k < M; ++k) {
            int j = i - tks[k];
            double uk = (j >= 0 && j < L) ? T[j] : 0.0;
            ATy[k] += uk * yi;
            for (int l = k; l < M; ++l) {
                int j2 = i - tks[l];
                double ul = (j2 >= 0 && j2 < L) ? T[j2] : 0.0;
                ATA[k][l] += uk * ul;
                if (l != k) ATA[l][k] = ATA[k][l];
            }
        }
    }

    std::vector<double> sol;
    if (!SolveNormalEq(ATA, ATy, sol, ridge)) return false;
    a_out.swap(sol);
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

// 各パルスの相関最大アンカーを固定し、その近傍だけで時刻を再調整する
static bool RefineTimesAroundAnchors(const std::vector<double>& x,
                                     const std::vector<double>& T,
                                     const std::vector<int>& anchors,
                                     int minSep,
                                     int half_width,
                                     double ridge,
                                     std::vector<int>& tks_out,
                                     double& b_out,
                                     std::vector<double>& a_out,
                                     double& rss_out) {
    if (anchors.empty()) return false;

    std::vector<int> tks = anchors;
    std::vector<double> a_fit;
    double b_fit = 0.0;
    if (!FitAmplitudes(x, T, tks, b_fit, a_fit, ridge)) return false;
    double rss = ComputeRSS(x, T, tks, b_fit, a_fit);

    const int kMaxPass = 3;
    for (int pass = 0; pass < kMaxPass; ++pass) {
        bool changed = false;
        for (size_t k = 0; k < tks.size(); ++k) {
            int best_tk = tks[k];
            std::vector<double> best_a = a_fit;
            double best_b = b_fit;
            double best_rss = rss;

            for (int d = -half_width; d <= half_width; ++d) {
                int tk_try = anchors[k] + d;
                bool ok_sep = true;
                for (size_t l = 0; l < tks.size(); ++l) {
                    if (l == k) continue;
                    if (std::abs(tks[l] - tk_try) < minSep) { ok_sep = false; break; }
                }
                if (!ok_sep) continue;

                std::vector<int> tks_try = tks;
                tks_try[k] = tk_try;
                std::vector<double> a_try;
                double b_try = 0.0;
                if (!FitAmplitudes(x, T, tks_try, b_try, a_try, ridge)) continue;
                double rss_try = ComputeRSS(x, T, tks_try, b_try, a_try);
                if (rss_try < best_rss) {
                    best_tk = tk_try;
                    best_a.swap(a_try);
                    best_b = b_try;
                    best_rss = rss_try;
                }
            }

            if (best_tk != tks[k]) {
                tks[k] = best_tk;
                a_fit.swap(best_a);
                b_fit = best_b;
                rss = best_rss;
                changed = true;
            }
        }
        if (!changed) break;
    }

    tks_out.swap(tks);
    a_out.swap(a_fit);
    b_out = b_fit;
    rss_out = rss;
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
    const bool has_baseline_csv = FindBaselineNoise(bnmap, event_id, baseline, sigma);
    if (!has_baseline_csv) {
        EstimateBaselineSigma(y, baseline, sigma);
    }
    if (dbg) {
        std::cerr << "[event " << event_id << "] baseline/sigma used: "
                  << baseline << ", " << sigma
                  << " (source=" << (has_baseline_csv ? "csv" : "estimated") << ")\n";
    }
    std::vector<double> x(y.size(), 0.0);
    for (size_t i = 0; i < y.size(); ++i) x[i] = y[i] - baseline;

    double sigmaCorr = sigma * std::sqrt(t_corr.sumT2);
    const double sigmaA = (t_full.sumT2 > 0.0) ? (sigma / std::sqrt(t_full.sumT2)) : sigma;
    const double aMin = kDefaultAMinSigmaFactor * sigmaA;
    const double thr = nSigma * sigmaCorr;

    std::vector<int> tks_full;
    std::vector<int> tks_anchor_full;
    std::vector<int> tks_crop;
    tks_full.reserve(maxK);
    tks_anchor_full.reserve(maxK);
    tks_crop.reserve(maxK);
    double b_fit = 0.0;
    if (!has_baseline_csv && !x.empty()) {
        // baseline CSV が無い場合のみ、残差平均で b_fit 初期値を与える
        double sx = 0.0;
        for (double v : x) sx += v;
        b_fit = sx / static_cast<double>(x.size());
    }
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
        std::vector<int> cand = SelectCandidates(c_res, sigmaCorr, kDefaultThr0Factor,
                                                 kDefaultNmsWindow, maxK);
        if (dbg) {
            std::cerr << "[event " << event_id << "] iter=" << iter
                      << " candidates=" << cand.size()
                      << " cmax=" << cmax << " cmin=" << cmin
                      << " sigmaCorr=" << sigmaCorr
                      << " thr=" << (nSigma * sigmaCorr) << "\n";
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
        bool found_best = false;
        int best_anchor = -1;
        double best_rss_iter = 0.0;
        double best_b_fit = 0.0;
        int best_s_full = 0;
        std::vector<int> best_tks_full;
        std::vector<double> best_a_fit;

        const size_t cand_end = cand.size();
        for (size_t idx = 0; idx < cand_end; ++idx) {
            int s_center = cand[idx];
            bool found_local = false;
            double local_best_rss = 0.0;
            double local_best_b = 0.0;
            std::vector<int> local_best_tks;
            std::vector<double> local_best_a;

            for (int delta = -kDefaultLocalRefineHalfWidth;
                 delta <= kDefaultLocalRefineHalfWidth; ++delta) {
                const int s_try = s_center + delta;
                if (s_try < 0 || s_try >= static_cast<int>(c_res.size())) continue;
                const int s_full_try = s_try - t_corr.crop_start;

                bool ok_sep = true;
                for (size_t k = 0; k < tks_full.size(); ++k) {
                    if (std::abs(tks_full[k] - s_full_try) < minSep) { ok_sep = false; break; }
                }
                if (!ok_sep) continue;

                std::vector<int> tks_full_try = tks_full;
                tks_full_try.push_back(s_full_try);
                std::vector<double> a_try;
                // 採用判定中はベースライン固定（b_fit）で、新規候補のみ時刻局所探索
                if (!FitAmplitudesFixedBaseline(x, t_full.T, tks_full_try, b_fit, a_try, ridge)) {
                    if (dbg) {
                        std::cerr << "[event " << event_id << "] fit failed (singular/ill-conditioned)\n";
                    }
                    continue;
                }

                if (a_try.empty()) continue;
                const double a_new = a_try.back();
                // 基準は相関最大点に固定（s_try ではなく s_center）
                if (!no_sign_cut && a_new * c_res[s_center] <= 0.0) {
                    continue;
                }
                if (tks_full_try.size() > 1 && std::fabs(a_new) < aMin) {
                    continue;
                }
                if (tks_full_try.size() > 1 && std::fabs(c_res[s_center]) < thr) {
                    continue;
                }

                const double rss_try = ComputeRSS(x, t_full.T, tks_full_try, b_fit, a_try);
                if (rss_prev > 0.0) {
                    const double improve = (rss_prev - rss_try) / rss_prev;
                    if (improve < rss_improve_min) {
                        continue;
                    }
                }

                if (!found_local || rss_try < local_best_rss) {
                    found_local = true;
                    local_best_rss = rss_try;
                    local_best_b = b_fit;
                    local_best_tks = tks_full_try;
                    local_best_a = a_try;
                }
            }

            if (!found_local) continue;

            if (!found_best || local_best_rss < best_rss_iter) {
                found_best = true;
                best_anchor = s_center;
                best_s_full = local_best_tks.back();
                best_rss_iter = local_best_rss;
                best_b_fit = local_best_b;
                best_tks_full = local_best_tks;
                best_a_fit = local_best_a;
                // iter=0 は相関順位を優先し、「最初に通った候補」を採用する
                if (kDefaultIter0RankFirst && iter == 0) break;
            }
        }

        if (found_best) {
            // 採用判定段階では「ベースライン固定 + 新規候補のみ局所調整」で更新
            // 採用後の全体リファイン（baseline + 全時刻）はループ終了後に1回だけ実施
            tks_anchor_full.push_back(best_s_full);
            tks_full.swap(best_tks_full);
            a_fit.swap(best_a_fit);
            b_fit = best_b_fit;
            rss_prev = best_rss_iter;
            added = true;
            if (dbg) {
                std::cerr << "[event " << event_id << "] accept(best-RSS) t_anchor=" << best_anchor
                          << " t_start_full(decision)=" << best_s_full
                          << " a=" << (a_fit.empty() ? 0.0 : a_fit.back())
                          << " rss=" << rss_prev << "\n";
            }
        }
        if (!added) break;

        ComputeResidual(x, t_full.T, tks_full, b_fit, a_fit, residual);
    }
    if (tks_full.empty()) return;

    // 採用完了後に1回だけ、ベースライン+全時刻を再調整
    {
        std::vector<int> tks_refined;
        std::vector<double> a_refined;
        double b_refined = 0.0;
        double rss_refined = 0.0;
        bool refined_ok = RefineTimesAroundAnchors(x, t_full.T, tks_anchor_full, minSep,
                                                   kDefaultPostAcceptRefineHalfWidth, ridge,
                                                   tks_refined, b_refined, a_refined, rss_refined);
        if (refined_ok) {
            tks_full.swap(tks_refined);
            a_fit.swap(a_refined);
            b_fit = b_refined;
            rss_prev = rss_refined;
        }
    }

    tks_crop.resize(tks_full.size());
    for (size_t k = 0; k < tks_full.size(); ++k) {
        tks_crop[k] = tks_full[k] + t_corr.crop_start;
    }

    const double baseline_fit = baseline + b_fit;
    int pulse_index = 0;
    for (size_t k = 0; k < tks_full.size(); ++k) {
        double a = a_fit[k];
        bool keep = true;
        if (!no_sign_cut) {
            double c_now = 0.0;
            if (ComputeCorrelationAtStart(x, t_corr.T, tks_crop[k], c_now)) {
                keep = (a * c_now > 0.0);
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

        if (dbg) {
            double cv = 0.0;
            ComputeCorrelationAtStart(x, t_corr.T, t_start, cv);
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
    int minSep = kDefaultMinSepSamples; // 既定はサンプル数で指定
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
    if (minSep < 1) minSep = 1;

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
