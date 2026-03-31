#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

// ============================================================
// PSsignal.cc
//
// 目的:
//  - event_baseline_noise.csv と pileup_waveforms.csv から
//    PS のパルスピーク時刻を抽出する
//
// 入力:
//  1) event_baseline_noise.csv
//     run,event_id,CH,baseline,noise
//     (旧形式: event_id,baseline,noise も許容)
//  2) pileup_waveforms.csv
//     1行1イベント: event_id, y0, y1, ..., y(N-1)
//
// 出力:
//  data/shapeddata/PSsignal.csv
//     run,event_id,pulse_index,peak_time
//
// 手順:
//  - baseline を差し引き、noise を用いて閾値を設定
//  - start閾値を下回ったらパルス開始
//  - end閾値を上回ったらパルス終了
//  - 開始〜終了の最小値の時刻をピーク時刻として出力
//  - ピークがなければ何も出力しない
// ============================================================

static const char* kDefaultBaselinePath = "data/shapeddata/event_baseline_noise.csv";
static const char* kDefaultWavePath = "data/shapeddata/pileup_waveforms.csv";
static const char* kDefaultOutPath = "data/shapeddata/PSsignal.csv";

// 閾値係数（noise に掛ける）
static const double kDefaultStartSigma = 9.0; // 立ち下がり開始判定
static const double kDefaultEndSigma = 9.0;   // 立ち下がり終了判定
// ピーク波高（|wave-baseline|）の下限
static const double kDefaultPeakHeightMin = 100.0;

struct BaselineNoise {
    long long event_id;
    int run_no;
    double baseline;
    double noise;
};

static std::string Trim(const std::string& s) {
    size_t a = 0;
    while (a < s.size() && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
    size_t b = s.size();
    while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1]))) --b;
    return s.substr(a, b - a);
}

static bool ParseLongLong(const std::string& s, long long& out) {
    char* endp = nullptr;
    errno = 0;
    long long v = std::strtoll(s.c_str(), &endp, 10);
    if (errno != 0 || endp == s.c_str() || *endp != '\0') return false;
    out = v;
    return true;
}

static bool ParseDouble(const std::string& s, double& out) {
    char* endp = nullptr;
    errno = 0;
    double v = std::strtod(s.c_str(), &endp);
    if (errno != 0 || endp == s.c_str() || *endp != '\0') return false;
    out = v;
    return true;
}

static std::vector<std::string> SplitCsv(const std::string& line) {
    std::vector<std::string> cols;
    std::stringstream ss(line);
    std::string token;
    while (std::getline(ss, token, ',')) {
        cols.push_back(Trim(token));
    }
    return cols;
}

static bool IsHeaderLine(const std::string& line) {
    if (line.find("run") != std::string::npos) return true;
    if (line.find("event") != std::string::npos) return true;
    if (line.find("baseline") != std::string::npos) return true;
    return false;
}

static bool LoadBaselineNoise(const std::string& path,
                              std::unordered_map<long long, BaselineNoise>& out_map,
                              int& run_no_global) {
    std::ifstream ifs(path.c_str());
    if (!ifs) return false;

    run_no_global = -1;
    std::string line;
    while (std::getline(ifs, line)) {
        line = Trim(line);
        if (line.empty()) continue;
        if (IsHeaderLine(line)) continue;

        std::vector<std::string> cols = SplitCsv(line);
        if (cols.size() < 3) continue;

        BaselineNoise bn;
        bn.run_no = -1;
        bn.event_id = -1;
        bn.baseline = 0.0;
        bn.noise = 0.0;

        // 新形式: run,event_id,CH,baseline,noise
        if (cols.size() >= 5) {
            long long run_ll = 0;
            if (!ParseLongLong(cols[0], run_ll)) continue;
            if (!ParseLongLong(cols[1], bn.event_id)) continue;
            if (!ParseDouble(cols[3], bn.baseline)) continue;
            if (!ParseDouble(cols[4], bn.noise)) continue;
            bn.run_no = static_cast<int>(run_ll);
        } else {
            // 旧形式: event_id,baseline,noise
            if (!ParseLongLong(cols[0], bn.event_id)) continue;
            if (!ParseDouble(cols[1], bn.baseline)) continue;
            if (!ParseDouble(cols[2], bn.noise)) continue;
        }

        if (run_no_global < 0 && bn.run_no >= 0) run_no_global = bn.run_no;
        out_map[bn.event_id] = bn;
    }
    return true;
}

static bool ReadWaveformsCsv(const std::string& path,
                             std::vector<long long>& event_ids,
                             std::vector<std::vector<double>>& waves) {
    std::ifstream ifs(path.c_str());
    if (!ifs) return false;

    std::string line;
    while (std::getline(ifs, line)) {
        line = Trim(line);
        if (line.empty()) continue;
        if (IsHeaderLine(line)) continue;

        std::vector<std::string> cols = SplitCsv(line);
        if (cols.size() < 2) continue;

        long long event_id = 0;
        if (!ParseLongLong(cols[0], event_id)) continue;

        std::vector<double> y;
        y.reserve(cols.size() - 1);
        bool ok = true;
        for (size_t i = 1; i < cols.size(); ++i) {
            double v = 0.0;
            if (!ParseDouble(cols[i], v)) { ok = false; break; }
            y.push_back(v);
        }
        if (!ok || y.empty()) continue;

        event_ids.push_back(event_id);
        waves.push_back(std::move(y));
    }
    return true;
}

static void DetectPulses(const std::vector<double>& y,
                         double baseline,
                         double noise,
                         double start_sigma,
                         double end_sigma,
                         std::vector<int>& peak_indices) {
    peak_indices.clear();
    if (y.empty()) return;
    if (noise <= 0.0) return;
    if (start_sigma <= 0.0 || end_sigma <= 0.0) return;

    double start_thr = -start_sigma * noise;
    double end_thr = -end_sigma * noise;

    bool in_pulse = false;
    int min_idx = -1;
    double min_val = 0.0;

    for (int i = 0; i < static_cast<int>(y.size()); ++i) {
        double v = y[i] - baseline;

        if (!in_pulse) {
            if (v <= start_thr) {
                in_pulse = true;
                min_val = v;
                min_idx = i;
            }
            continue;
        }

        if (v < min_val) {
            min_val = v;
            min_idx = i;
        }

        if (v >= end_thr) {
            if (min_idx >= 0) peak_indices.push_back(min_idx);
            in_pulse = false;
            min_idx = -1;
        }
    }

    // 終端まで戻らなかった場合は出力しない
}

static void Usage(const char* prog) {
    std::cerr
        << "Usage:\n"
        << "  " << prog << " [event_baseline_noise.csv pileup_waveforms.csv out.csv]\n"
        << "  Options:\n"
        << "    --start-sigma S      立ち下がり開始閾値 (default: 9.0)\n"
        << "    --end-sigma E        立ち下がり終了閾値 (default: 9.0)\n"
        << "    --peak-height-min H  波高下限 |wave-baseline| (default: 100)\n"
        << "\n";
}

int main(int argc, char** argv) {
    std::string baseline_path = kDefaultBaselinePath;
    std::string wave_path = kDefaultWavePath;
    std::string out_path = kDefaultOutPath;
    double start_sigma = kDefaultStartSigma;
    double end_sigma = kDefaultEndSigma;
    double peak_height_min = kDefaultPeakHeightMin;

    std::vector<std::string> args;
    for (int i = 1; i < argc; ++i) args.push_back(argv[i]);

    for (size_t i = 0; i < args.size(); ++i) {
        const std::string& a = args[i];
        if (a == "-h" || a == "--help") {
            Usage(argv[0]);
            return 0;
        } else if (a == "--start-sigma" && i + 1 < args.size()) {
            start_sigma = std::atof(args[++i].c_str());
        } else if (a == "--end-sigma" && i + 1 < args.size()) {
            end_sigma = std::atof(args[++i].c_str());
        } else if (a == "--peak-height-min" && i + 1 < args.size()) {
            peak_height_min = std::atof(args[++i].c_str());
        } else if (!a.empty() && a[0] == '-') {
            std::cerr << "Unknown option: " << a << "\n";
            Usage(argv[0]);
            return 1;
        } else {
            // positional
            if (baseline_path == kDefaultBaselinePath) baseline_path = a;
            else if (wave_path == kDefaultWavePath) wave_path = a;
            else if (out_path == kDefaultOutPath) out_path = a;
            else {
                std::cerr << "Too many positional arguments.\n";
                Usage(argv[0]);
                return 1;
            }
        }
    }

    if (end_sigma > start_sigma) {
        std::cerr << "warn: end_sigma > start_sigma. swapping to keep hysteresis.\n";
        std::swap(end_sigma, start_sigma);
    }

    std::unordered_map<long long, BaselineNoise> bnmap;
    int run_no_global = -1;
    if (!LoadBaselineNoise(baseline_path, bnmap, run_no_global)) {
        std::cerr << "error: failed to load baseline file: " << baseline_path << "\n";
        return 1;
    }

    std::vector<long long> event_ids;
    std::vector<std::vector<double>> waves;
    if (!ReadWaveformsCsv(wave_path, event_ids, waves)) {
        std::cerr << "error: failed to load waveform csv: " << wave_path << "\n";
        return 1;
    }

    std::error_code ec;
    std::filesystem::path out_parent = std::filesystem::path(out_path).parent_path();
    if (!out_parent.empty()) {
        std::filesystem::create_directories(out_parent, ec);
        if (ec) {
            std::cerr << "error: failed to create output dir: " << ec.message() << "\n";
            return 1;
        }
    }

    std::ofstream ofs(out_path.c_str());
    if (!ofs) {
        std::cerr << "error: failed to open output: " << out_path << "\n";
        return 1;
    }

    ofs << "run,event_id,pulse_index,peak_time\n";
    ofs << std::setprecision(12);

    for (size_t i = 0; i < waves.size(); ++i) {
        long long event_id = event_ids[i];
        const std::vector<double>& y = waves[i];

        auto it = bnmap.find(event_id);
        if (it == bnmap.end()) {
            std::cerr << "warn: baseline not found for event " << event_id << " (skip)\n";
            continue;
        }

        const BaselineNoise& bn = it->second;
        int run_no = (bn.run_no >= 0) ? bn.run_no : run_no_global;

        std::vector<int> peaks;
        DetectPulses(y, bn.baseline, bn.noise, start_sigma, end_sigma, peaks);
        if (peaks.empty()) continue;

        int pulse_index = 0;
        for (int t_peak : peaks) {
            if (t_peak < 0 || t_peak >= static_cast<int>(y.size())) continue;
            const double peak_height = std::fabs(y[static_cast<size_t>(t_peak)] - bn.baseline);
            if (peak_height <= peak_height_min) continue; // 既定では「100を超える」
            ofs << run_no << "," << event_id << "," << pulse_index << "," << t_peak << "\n";
            ++pulse_index;
        }
    }

    return 0;
}
