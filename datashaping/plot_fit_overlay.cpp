#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <iostream>
#include <cstdlib>

// 既定パラメータ
static const int kBaselineEstN = 50;
static const bool kDefaultRawTime = true;
// 既定パス
static const char* kDefaultTemplatePath = "data/shapeddata/template.json";
static const char* kDefaultBaselinePath = "data/shapeddata/event_baseline_noise.csv";
static const char* kDefaultWavePath = "data/shapeddata/pileup_waveforms.csv";
static const char* kDefaultPulsesPath = "data/shapeddata/pulses.csv";
static const char* kDefaultOutPath = "";
static const long long kDefaultEventId = 0;

// ============================================================
// plot_fit_overlay.cpp
//
// ビルド例:
//   g++ -O2 -std=c++17 -o plot_fit_overlay datashaping/plot_fit_overlay.cpp
//
// 入力:
//   1) template.json
//      keys: dt, polarity, template[], areaT, sumT2, t0_index, ...
//   2) event_baseline_noise.csv
//      run,event_id,CH,baseline,noise
//      (旧形式: event_id,baseline,sigma も許容)
//   3) pileup_waveforms.csv
//      1行1イベント: event_id, y0, y1, ..., y(N-1)
//   4) pulses.csv
//      run,event_id,pulse_index,area,t_fall_start_sample_raw,amplitude,baseline_fit
//
// 出力:
//   out.csv
//     time, y_raw, y_fit
//
// コマンドライン:
//   plot_fit_overlay [template.json event_baseline_noise.csv pileup_waveforms.csv pulses.csv event_id out.csv]
//   plot_fit_overlay ... --debug
//   plot_fit_overlay ... --plot-png out.png
//   plot_fit_overlay ... --template-time
//
//
// 注意:
//   - time = (sample_index - t0_index) * dt
//   - y_fit は baseline + Σ a_k T[i - t_k]
// ============================================================

struct TemplateData {
    double dt;
    double polarity;
    double areaT;
    double sumT2;
    int t0_index;
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
    std::stringstream ss(s.substr(p, end - p));
    ss >> out;
    return !ss.fail();
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
    return true;
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
        std::stringstream s0(col);
        s0 >> bn.event_id;
        if (s0.fail()) {
            // 新形式: run,event_id,CH,baseline,noise
            // 旧形式: event_id,baseline,sigma
            // ヘッダ or 非数値行はスキップ
            // new format: 先頭は run なので event_id は次の列
            if (!std::getline(ss, col, ',')) continue;
            std::stringstream s1(col);
            s1 >> bn.event_id;
            if (s1.fail()) continue;
            // CH
            if (!std::getline(ss, col, ',')) continue;
            // baseline
            if (!std::getline(ss, col, ',')) continue;
            std::stringstream s2(col);
            s2 >> bn.baseline;
            if (s2.fail()) continue;
            // noise
            if (!std::getline(ss, col, ',')) continue;
            std::stringstream s3(col);
            s3 >> bn.sigma;
            if (s3.fail()) continue;
            out.push_back(bn);
            continue;
        }
        // 旧形式: event_id,baseline,sigma
        if (!std::getline(ss, col, ',')) continue;
        std::stringstream s1(col);
        s1 >> bn.baseline;
        if (s1.fail()) continue;
        if (!std::getline(ss, col, ',')) continue;
        std::stringstream s2(col);
        s2 >> bn.sigma;
        if (s2.fail()) continue;
        out.push_back(bn);
    }
    return true;
}

static bool FindBaselineNoise(const std::vector<BaselineNoise>& v, long long event_id, double& baseline, double& sigma) {
    for (size_t i = 0; i < v.size(); ++i) {
        if (v[i].event_id == event_id) {
            baseline = v[i].baseline;
            sigma = v[i].sigma;
            return true;
        }
    }
    return false;
}

static void EstimateBaseline(const std::vector<double>& y, double& baseline) {
    const int n0 = (y.size() < static_cast<size_t>(kBaselineEstN))
                       ? static_cast<int>(y.size())
                       : kBaselineEstN;
    if (n0 <= 0) { baseline = 0.0; return; }
    double sum = 0.0;
    for (int i = 0; i < n0; ++i) sum += y[i];
    baseline = sum / n0;
}

static bool LoadWaveformForEvent(const std::string& path, long long event_id, std::vector<double>& y) {
    std::ifstream ifs(path.c_str());
    if (!ifs) return false;
    std::string line;
    while (std::getline(ifs, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string col;
        if (!std::getline(ss, col, ',')) continue;
        long long eid = 0;
        std::stringstream s0(col);
        s0 >> eid;
        if (eid != event_id) continue;
        y.clear();
        while (std::getline(ss, col, ',')) {
            double v = 0.0;
            std::stringstream s1(col);
            s1 >> v;
            y.push_back(v);
        }
        return !y.empty();
    }
    return false;
}

static bool LoadPulsesForEvent(const std::string& path, long long event_id, int t0_index,
                               std::vector<Pulse>& pulses,
                               double& baseline_fit, bool& has_baseline_fit) {
    std::ifstream ifs(path.c_str());
    if (!ifs) return false;
    pulses.clear();
    has_baseline_fit = false;
    baseline_fit = 0.0;
    std::string line;
    while (std::getline(ifs, line)) {
        if (line.empty()) continue;
        // header skip
        if (line.find("event_id") != std::string::npos) continue;
        std::stringstream ss(line);
        std::string col;
        std::vector<std::string> cols;
        while (std::getline(ss, col, ',')) cols.push_back(col);
        if (cols.size() < 5) continue;
        long long eid = 0;
        std::stringstream s0(cols[0]);
        s0 >> eid;
        bool has_run = !s0.fail();
        if (!has_run) {
            continue;
        }
        // 新形式: run,event_id,pulse_index,area,t_fall_start_sample_raw,amplitude,baseline_fit
        if (cols.size() >= 7) {
            std::stringstream s1(cols[1]);
            s1 >> eid;
        }
        if (eid != event_id) continue;
        int t_start_full = 0;
        double a = 0.0;
        if (cols.size() == 7) {
            // new format: run,event_id,pulse_index,area,t_fall_start_sample_raw,amplitude,baseline_fit
            double t_fall_start_sample_raw = 0.0;
            std::stringstream s2(cols[4]);
            s2 >> t_fall_start_sample_raw;
            t_start_full = static_cast<int>(std::floor(t_fall_start_sample_raw - t0_index + 0.5));
            std::stringstream s3(cols[5]);
            s3 >> a;
            if (!has_baseline_fit) {
                std::stringstream s4(cols[6]);
                s4 >> baseline_fit;
                if (!s4.fail()) has_baseline_fit = true;
            }
        } else if (cols.size() >= 10) {
            // old format: event_id,pulse_index,t_start_sample,t_start_sample_full,t_start_time,amplitude,...,baseline_fit
            std::stringstream s2(cols[3]);
            s2 >> t_start_full;
            std::stringstream s3(cols[5]);
            s3 >> a;
            if (!has_baseline_fit) {
                std::stringstream s4(cols[9]);
                s4 >> baseline_fit;
                if (!s4.fail()) has_baseline_fit = true;
            }
        } else if (cols.size() >= 9) {
            // mid format: event_id,pulse_index,t_start_sample,t_start_sample_full,t_start_time,amplitude,...
            std::stringstream s2(cols[3]);
            s2 >> t_start_full;
            std::stringstream s3(cols[5]);
            s3 >> a;
        } else {
            // old format: event_id,pulse_index,t_start_sample,t_start_time,amplitude,...
            std::stringstream s3(cols[4]);
            s3 >> a;
            t_start_full = 0;
        }
        Pulse p;
        p.t_start = t_start_full;
        p.amplitude = a;
        pulses.push_back(p);
    }
    // ファイルを正常に読めたら true（当該 event のパルス0件も許容）
    return true;
}

static void BuildFit(const std::vector<double>& y, const TemplateData& t, double baseline,
                     const std::vector<Pulse>& pulses, std::vector<double>& yfit) {
    const int N = static_cast<int>(y.size());
    const int L = static_cast<int>(t.T.size());
    yfit.assign(N, baseline);
    for (size_t k = 0; k < pulses.size(); ++k) {
        int tk = pulses[k].t_start;
        double a = pulses[k].amplitude;
        for (int i = 0; i < N; ++i) {
            int j = i - tk;
            if (j < 0 || j >= L) continue;
            yfit[i] += a * t.T[j];
        }
    }
}

int main(int argc, char** argv) {
    bool debug = false;
    bool plot_png = false;
    std::string plot_png_path;
    bool raw_time = kDefaultRawTime;
    std::vector<std::string> positional;
    positional.reserve(argc);
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--debug") {
            debug = true;
        } else if (a == "--plot-png") {
            plot_png = true;
            if (i + 1 < argc) {
                std::string next = argv[i + 1];
                if (next.rfind("--", 0) != 0) {
                    plot_png_path = argv[++i];
                }
            }
        } else if (a == "--template-time") {
            raw_time = false;
        } else if (a.rfind("--", 0) == 0) {
            // unknown option -> ignore
        } else {
            positional.push_back(a);
        }
    }

    std::string template_path = kDefaultTemplatePath;
    std::string baseline_path = kDefaultBaselinePath;
    std::string wave_path = kDefaultWavePath;
    std::string pulses_path = kDefaultPulsesPath;
    long long event_id = kDefaultEventId;
    std::string out_path = kDefaultOutPath;

    auto is_numeric = [](const std::string& s) -> bool {
        if (s.empty()) return false;
        if (!(std::isdigit(s[0]) || s[0] == '-' || s[0] == '+')) return false;
        for (size_t i = 1; i < s.size(); ++i) {
            if (!std::isdigit(s[i])) return false;
        }
        return true;
    };

    if (positional.size() == 1 && is_numeric(positional[0])) {
        std::stringstream ss(positional[0]);
        ss >> event_id;
    } else if ((positional.size() == 2 || positional.size() == 1) && is_numeric(positional[0])) {
        std::stringstream ss(positional[0]);
        ss >> event_id;
        if (positional.size() == 2) out_path = positional[1];
    } else {
        if (positional.size() >= 1) template_path = positional[0];
        if (positional.size() >= 2) baseline_path = positional[1];
        if (positional.size() >= 3) wave_path = positional[2];
        if (positional.size() >= 4) pulses_path = positional[3];
        if (positional.size() >= 5) {
            std::stringstream ss(positional[4]);
            ss >> event_id;
        }
        if (positional.size() >= 6) out_path = positional[5];
    }

    if (out_path.empty()) {
        out_path = "data/shapeddata/overlay_event" + std::to_string(event_id) + ".csv";
    }
    if (plot_png && plot_png_path.empty()) {
        plot_png_path = "data/shapeddata/overlay_event" + std::to_string(event_id) + ".png";
    }

    TemplateData t;
    if (!LoadTemplate(template_path, t)) {
        std::cerr << "error: failed to load template: " << template_path << "\n";
        return 1;
    }
    if (debug) {
        std::cerr << "info: template loaded. L=" << t.T.size()
                  << " dt=" << t.dt << " t0_index=" << t.t0_index << "\n";
    }

    std::vector<double> y;
    if (!LoadWaveformForEvent(wave_path, event_id, y)) {
        std::cerr << "error: event_id " << event_id << " not found in " << wave_path << "\n";
        return 1;
    }
    if (debug) {
        std::cerr << "info: waveform loaded. N=" << y.size() << "\n";
    }

    std::vector<BaselineNoise> bnlist;
    if (!LoadBaselineNoise(baseline_path, bnlist)) {
        std::cerr << "warn: failed to load baseline list: " << baseline_path << "\n";
    }

    double baseline = 0.0, sigma = 0.0;
    if (!FindBaselineNoise(bnlist, event_id, baseline, sigma)) {
        EstimateBaseline(y, baseline);
        if (debug) {
            std::cerr << "info: baseline estimated: " << baseline << "\n";
        }
    } else if (debug) {
        std::cerr << "info: baseline found: " << baseline << " sigma=" << sigma << "\n";
    }

    std::vector<Pulse> pulses;
    double baseline_fit = 0.0;
    bool has_baseline_fit = false;
    if (!LoadPulsesForEvent(pulses_path, event_id, t.t0_index, pulses, baseline_fit, has_baseline_fit)) {
        std::cerr << "error: failed to open pulses file: " << pulses_path << "\n";
        return 1;
    }
    if (debug) {
        std::cerr << "info: pulses loaded. M=" << pulses.size() << "\n";
        if (pulses.empty()) {
            std::cerr << "info: no pulse info for event " << event_id
                      << ". y_fit is baseline only.\n";
        }
    }
    if (has_baseline_fit) {
        baseline = baseline_fit;
        if (debug) {
            std::cerr << "info: baseline overridden by fit: " << baseline << "\n";
        }
    }

    std::vector<double> yfit;
    BuildFit(y, t, baseline, pulses, yfit);
    if (debug) {
        std::cerr << "info: fit built. yfit size=" << yfit.size() << "\n";
    }

    std::ofstream ofs(out_path.c_str());
    if (!ofs) {
        std::cerr << "error: failed to open output: " << out_path << "\n";
        return 1;
    }
    ofs << "time,y_raw,y_fit\n";
    for (int i = 0; i < static_cast<int>(y.size()); ++i) {
        double time = raw_time ? (i * t.dt) : ((i - t.t0_index) * t.dt);
        ofs << time << "," << y[i] << "," << yfit[i] << "\n";
    }
    ofs.close();

    if (plot_png) {
        // gnuplot が入っている場合のみPNGを出力
        std::string cmd = "gnuplot -e \"infile='" + out_path +
                          "'; outfile='" + plot_png_path +
                          "'\" datashaping/plot_fit_overlay.gp";
        int ret = std::system(cmd.c_str());
        if (ret != 0) {
            std::cerr << "error: gnuplot failed (ret=" << ret << ")\n";
            return 1;
        }
    }
    return 0;
}
