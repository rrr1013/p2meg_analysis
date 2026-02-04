// step1.cc
// Build dataset (run, event, CH, integral, peak_time) from 6 waveform txt files.
//
// Assumptions:
// - Each input file contains one number per line (one bin per line).
// - Each event consists of N bins (bins_per_event). Default N=500, can be overridden by --bins-per-event.
// - Waveform is a falling pulse (negative-going). We compute integral as sum(baseline - wf[t]).
// - peak_time is the bin index where wf is minimum (0-based by default; can output 1-based with --one-based-peak).
//
// Input files (example):
//   wave_PS_A_run12.txt
//   wave_PS_B_run12.txt
//   wave_NaI_A1_run12.txt
//   wave_NaI_A2_run12.txt
//   wave_NaI_B1_run12.txt
//   wave_NaI_B2_run12.txt
//
// Usage examples:
//   g++ -O2 -std=c++17 step1.cc -o step1
//   ./step1 --input-dir ./data
//   ./step1 --input-dir ./data --bins-per-event 1024 --baseline-bins 50
//   ./step1 --input-dir ./data --integ-start 0 --integ-end 400 --one-based-peak
//
// Output:
//   step1_runM.txt (tab-separated) in output-dir (or specified by --out)

#include <algorithm>
#include <cerrno>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <regex>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace fs = std::filesystem;

static const std::vector<std::string> CH_ORDER = {
    "PS_A", "PS_B", "NaI_A1", "NaI_A2", "NaI_B1", "NaI_B2"
};

struct Config {
    std::string input_dir = "../data/rawdata";
    std::string output_dir = "../data/shapeddata";
    std::string out_file;          // empty => auto step1_runM.txt in output-dir
    int bins_per_event = 500;      // DEFAULT = 500 (requested)
    int baseline_bins  = 50;
    int integ_start    = 0;
    int integ_end      = -1;       // -1 => to end
    double ps_integ_threshold = 100.0;   // <= threshold => zero integral, peak_time=bins_per_event
    double nai_integ_threshold = 2000.0; // <= threshold => zero integral, peak_time=bins_per_event
    bool one_based_peak = false;
    bool michel_mode = false;
};

static void usage(const char* prog) {
    std::cerr
        << "Usage:\n"
        << "  " << prog << " [options]\n\n"
        << "This program auto-detects run number from filenames in the input directory.\n"
        << "It expects exactly one consistent run for the 6 required channels.\n\n"
        << "Options:\n"
        << "  --input-dir DIR        input directory (default: ../data/rawdata)\n"
        << "  --output-dir DIR       output directory (default: ../data/shapeddata)\n"
        << "  --out FILE             output txt file (default: step1_runM.txt in output-dir)\n"
        << "  --bins-per-event N      bins per event (default: 500)\n"
        << "  --baseline-bins K       baseline bins from head (default: 50)\n"
        << "  --integ-start S         integration start bin (default: 0)\n"
        << "  --integ-end E           integration end bin (exclusive). -1 means end (default: -1)\n"
        << "  --ps-integ-threshold T  PS: if integral <= T, set integral=0 and peak_time=bins_per_event (default: 100)\n"
        << "  --nai-integ-threshold T NaI: if integral <= T, set integral=0 and peak_time=bins_per_event (default: 2000)\n"
        << "  --one-based-peak        output peak_time as 1-based (default: 0-based)\n"
        << "  --michel               Michel mode (omit peak_time column)\n"
        << "  -h, --help              show this help\n";
}

static std::vector<double> load_bins(const std::string& path) {
    std::ifstream ifs(path);
    if (!ifs) throw std::runtime_error("Failed to open input file: " + path);

    std::vector<double> bins;
    bins.reserve(1 << 20);

    std::string line;
    long long line_no = 0;
    while (std::getline(ifs, line)) {
        ++line_no;

        // trim spaces
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

static std::vector<std::vector<double>> split_events(const std::vector<double>& all_bins, int bins_per_event) {
    if (bins_per_event <= 0) throw std::runtime_error("bins_per_event must be > 0");

    const size_t n = all_bins.size();
    const size_t N = static_cast<size_t>(bins_per_event);
    if (n % N != 0) {
        throw std::runtime_error(
            "Total bins (" + std::to_string(n) + ") is not divisible by bins_per_event (" +
            std::to_string(bins_per_event) + "). Check N or missing/extra lines."
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

static double baseline_mean(const std::vector<double>& wf, int baseline_bins) {
    int k = std::max(1, baseline_bins);
    k = std::min(k, static_cast<int>(wf.size()));
    double sum = std::accumulate(wf.begin(), wf.begin() + k, 0.0);
    return sum / static_cast<double>(k);
}

struct Features {
    double integral;
    int peak_idx;  // index of minimum value in wf
};

static Features compute_features(const std::vector<double>& wf,
                                 int baseline_bins,
                                 int integ_start,
                                 int integ_end) {
    const int n = static_cast<int>(wf.size());
    if (n <= 0) throw std::runtime_error("Empty waveform");

    const double bl = baseline_mean(wf, baseline_bins);

    int s = std::max(0, integ_start);
    int e = (integ_end < 0) ? n : std::min(n, integ_end);
    if (s >= e) {
        throw std::runtime_error(
            "Invalid integration window: start=" + std::to_string(integ_start) +
            " end=" + std::to_string(integ_end) + " n=" + std::to_string(n)
        );
    }

    double integral = 0.0;
    for (int t = s; t < e; ++t) {
        integral += (bl - wf[t]); // falling pulse => positive area
    }

    double min_val = *std::min_element(wf.begin(), wf.end());
    long long sum_idx = 0;
    int count = 0;
    for (int i = 0; i < n; ++i) {
        if (wf[i] == min_val) {
            sum_idx += i;
            ++count;
        }
    }
    int peak_idx = (count > 0) ? static_cast<int>((sum_idx + count / 2) / count) : 0;
    return {integral, peak_idx};
}

// Scan input directory and find required files.
// Returns: run_number, and map ch->filepath
static std::pair<int, std::unordered_map<std::string, std::string>>
discover_files_and_run(const std::string& input_dir, bool michel_mode) {
    // Match:
    //  - default : wave_<CH>_run<M>.txt
    //  - michel  : Michel_<CH>_run<M>.txt (例: Michel_NaI_A1_run12.txt)
    // CH could contain letters/numbers; we only accept the required CH list.
    const std::string prefix = michel_mode ? "Michel_" : "wave_";
    std::regex re(("^" + prefix + "([A-Za-z0-9_]+)_run([0-9]+)\\.txt$"));

    std::unordered_map<std::string, std::string> ch_to_path;
    ch_to_path.reserve(CH_ORDER.size());

    int run_detected = -1;
    bool run_set = false;

    for (const auto& entry : fs::directory_iterator(input_dir)) {
        if (!entry.is_regular_file()) continue;

        const std::string fname = entry.path().filename().string();
        std::smatch m;
        if (!std::regex_match(fname, m, re)) continue;

        const std::string ch = m[1].str();
        const int run = std::stoi(m[2].str());

        // Only accept required channels
        if (std::find(CH_ORDER.begin(), CH_ORDER.end(), ch) == CH_ORDER.end()) continue;

        if (!run_set) {
            run_detected = run;
            run_set = true;
        } else if (run != run_detected) {
            // If multiple runs exist, this is ambiguous and dangerous
            throw std::runtime_error(
                "Multiple run numbers found in input-dir. Found run" + std::to_string(run_detected) +
                " and run" + std::to_string(run) +
                ". Please keep one run per directory or extend selection logic."
            );
        }

        // If duplicates, keep the first and warn
        if (ch_to_path.find(ch) == ch_to_path.end()) {
            ch_to_path[ch] = entry.path().string();
        }
    }

    if (!run_set) {
        throw std::runtime_error(
            "No matching files found in input-dir. Expected files like wave_PS_A_run12.txt, "
            "wave_NaI_A1_run12.txt (or Michel_NaI_A1_run12.txt in Michel mode), etc."
        );
    }

    // Ensure all required channels exist
    for (const auto& ch : CH_ORDER) {
        if (ch_to_path.find(ch) == ch_to_path.end()) {
            throw std::runtime_error("Missing required channel file for " + ch +
                                     " in input-dir: " + input_dir);
        }
    }

    return {run_detected, ch_to_path};
}

static Config parse_args(int argc, char** argv) {
    Config cfg;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];

        auto need = [&](const std::string& opt) -> std::string {
            if (i + 1 >= argc) throw std::runtime_error("Missing value for " + opt);
            return std::string(argv[++i]);
        };

        if (a == "--input-dir") cfg.input_dir = need(a);
        else if (a == "--output-dir") cfg.output_dir = need(a);
        else if (a == "--out") cfg.out_file = need(a);
        else if (a == "--bins-per-event") cfg.bins_per_event = std::stoi(need(a));
        else if (a == "--baseline-bins") cfg.baseline_bins = std::stoi(need(a));
        else if (a == "--integ-start") cfg.integ_start = std::stoi(need(a));
        else if (a == "--integ-end") cfg.integ_end = std::stoi(need(a));
        else if (a == "--ps-integ-threshold") cfg.ps_integ_threshold = std::stod(need(a));
        else if (a == "--nai-integ-threshold") cfg.nai_integ_threshold = std::stod(need(a));
        else if (a == "--one-based-peak") cfg.one_based_peak = true;
        else if (a == "--michel") cfg.michel_mode = true;
        else if (a == "-h" || a == "--help") { usage(argv[0]); std::exit(0); }
        else throw std::runtime_error("Unknown argument: " + a);
    }

    if (cfg.bins_per_event <= 0) throw std::runtime_error("--bins-per-event must be > 0");
    if (cfg.baseline_bins <= 0) throw std::runtime_error("--baseline-bins must be > 0");

    return cfg;
}

int main(int argc, char** argv) {
    try {
        Config cfg = parse_args(argc, argv);

        // 1) Discover run number and required files from filenames
        auto [run, ch_to_path] = discover_files_and_run(cfg.input_dir, cfg.michel_mode);

        // 2) Load and split per channel
        std::unordered_map<std::string, std::vector<std::vector<double>>> ch_events;
        ch_events.reserve(CH_ORDER.size());

        size_t nevt_ref = 0;
        bool nevt_set = false;

        for (const auto& ch : CH_ORDER) {
            const std::string& path = ch_to_path.at(ch);
            auto bins = load_bins(path);
            auto events = split_events(bins, cfg.bins_per_event);

            if (!nevt_set) {
                nevt_ref = events.size();
                nevt_set = true;
            } else if (events.size() != nevt_ref) {
                throw std::runtime_error(
                    "Event count mismatch across channels. " + ch + " has " +
                    std::to_string(events.size()) + " events, expected " +
                    std::to_string(nevt_ref) + ". Check data consistency."
                );
            }
            ch_events.emplace(ch, std::move(events));
        }

        // 3) Decide output file path
        std::string out_path;
        if (!cfg.out_file.empty()) {
            out_path = cfg.out_file;
        } else if (cfg.michel_mode) {
            out_path = (fs::path(cfg.output_dir) / ("step1_michel_run" + std::to_string(run) + ".txt")).string();
        } else {
            out_path = (fs::path(cfg.output_dir) / ("step1_run" + std::to_string(run) + ".txt")).string();
        }

        // 4) Write output
        std::ofstream ofs(out_path);
        if (!ofs) throw std::runtime_error("Failed to open output file: " + out_path);

        if (cfg.michel_mode) {
            ofs << "run\tevent\tCH\tintegral\n";
        } else {
            ofs << "run\tevent\tCH\tintegral\tpeak_time\n";
        }
        ofs << std::fixed << std::setprecision(6);

        for (size_t evt = 0; evt < nevt_ref; ++evt) {
            for (const auto& ch : CH_ORDER) {
                const auto& wf = ch_events.at(ch).at(evt);
                Features feat = compute_features(wf, cfg.baseline_bins, cfg.integ_start, cfg.integ_end);
                double integral = feat.integral;
                int peak_time = cfg.one_based_peak ? (feat.peak_idx + 1) : feat.peak_idx;
                const bool is_ps = (ch.rfind("PS", 0) == 0);
                const double threshold = is_ps ? cfg.ps_integ_threshold : cfg.nai_integ_threshold;
                if (integral <= threshold) {
                    integral = 0.0;
                    peak_time = cfg.bins_per_event;
                }

                if (cfg.michel_mode) {
                    ofs << run << "\t" << evt << "\t" << ch
                        << "\t" << integral << "\n";
                } else {
                    ofs << run << "\t" << evt << "\t" << ch
                        << "\t" << integral << "\t" << peak_time << "\n";
                }
            }
        }

        std::cout
            << "Run detected: " << run << "\n"
            << "bins_per_event: " << cfg.bins_per_event << " (default=500 unless overridden)\n"
            << "events: " << nevt_ref << ", channels: " << CH_ORDER.size() << "\n"
            << "wrote: " << out_path << "\n";

        return 0;

    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }
}

