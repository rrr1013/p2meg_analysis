// eventdisplay.cc
// Visualize per-event waveforms and energy (NaI) / hit (PS) with gnuplot.
//
// Input files (example):
//   wave_PS_A_run12.txt, wave_PS_B_run12.txt,
//   wave_NaI_A1_run12.txt, wave_NaI_A2_run12.txt,
//   wave_NaI_B1_run12.txt, wave_NaI_B2_run12.txt
//
// Build:
//   g++ -O2 -std=c++17 src/eventdisplay.cc -o eventdisplay
//
// Run:
//   ./eventdisplay --input-dir ./data/rawdata
//   ./eventdisplay --bins-per-event 500 --baseline-bins 50
//   ./eventdisplay --nai-a1-coeff 1.0 --nai-a1-offset 0.0 --nai-a2-coeff 1.0 --nai-a2-offset 0.0 \
//                 --nai-b1-coeff 1.0 --nai-b1-offset 0.0 --nai-b2-coeff 1.0 --nai-b2-offset 0.0
//
// Notes:
// - Energy per NaI channel = coeff * integral + offset (same coefficients as step2).
// - NaI integrals below --nai-integ-threshold are set to 0.
// - PS does not use energy conversion; only HIT/NO based on amplitude (baseline - min).
// - Events are displayed one by one; press Enter for next, or 'q' + Enter to quit.

#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace fs = std::filesystem;

static const std::vector<std::string> CH_ORDER = {
    "NaI_A1", "NaI_A2", "NaI_B1", "NaI_B2", "PS_A", "PS_B"
};

struct Config {
    std::string input_dir = "../data/rawdata";
    int bins_per_event = 500;
    int baseline_bins = 50;
    int integ_start = 0;
    int integ_end = -1; // -1 => to end
    bool michel_mode = false;
    int run_filter = -1;

    double nai_a1_coeff = 1.0;
    double nai_a1_offset = 0.0;
    double nai_a2_coeff = 1.0;
    double nai_a2_offset = 0.0;
    double nai_b1_coeff = 1.0;
    double nai_b1_offset = 0.0;
    double nai_b2_coeff = 1.0;
    double nai_b2_offset = 0.0;

    double ps_amp_threshold = 50.0; // PS: amplitude threshold (baseline - min), unit: ADC count
    double nai_integ_threshold = 500.0;

    int start_event = 0;
    int max_events = -1; // -1 => all
};

static void usage(const char* prog) {
    std::cerr
        << "Usage:\n"
        << "  " << prog << " [options]\n\n"
        << "Options:\n"
        << "  --input-dir DIR          input directory (default: ../data/rawdata)\n"
        << "  --bins-per-event N        bins per event (default: 500)\n"
        << "  --baseline-bins K         baseline bins (default: 50)\n"
        << "  --integ-start S           integration start bin (default: 0)\n"
        << "  --integ-end E             integration end bin (exclusive), -1 means end\n"
        << "  --michel                  read Michel_*.txt files (default: wave_*.txt)\n"
        << "  --run N                   select run number N (default: auto-detect)\n"
        << "  --nai-a1-coeff C          coeff for NaI_A1 (default: 1.0)\n"
        << "  --nai-a1-offset B         offset for NaI_A1 (default: 0.0)\n"
        << "  --nai-a2-coeff C          coeff for NaI_A2 (default: 1.0)\n"
        << "  --nai-a2-offset B         offset for NaI_A2 (default: 0.0)\n"
        << "  --nai-b1-coeff C          coeff for NaI_B1 (default: 1.0)\n"
        << "  --nai-b1-offset B         offset for NaI_B1 (default: 0.0)\n"
        << "  --nai-b2-coeff C          coeff for NaI_B2 (default: 1.0)\n"
        << "  --nai-b2-offset B         offset for NaI_B2 (default: 0.0)\n"
        << "  --ps-amp-threshold T      PS amplitude threshold (baseline - min, default: 100)\n"
        << "  --nai-integ-threshold T   NaI integral threshold (default: 2000)\n"
        << "  --start-event N           start event index (default: 0)\n"
        << "  --max-events N            max events to display (default: all)\n"
        << "  -h, --help                show this help\n";
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

static double baseline_mean(const std::vector<double>& wf, int baseline_bins) {
    int k = std::max(1, baseline_bins);
    k = std::min(k, static_cast<int>(wf.size()));
    double sum = std::accumulate(wf.begin(), wf.begin() + k, 0.0);
    return sum / static_cast<double>(k);
}

struct Features {
    double integral = 0.0;
    int peak_idx = 0;
    double amplitude = 0.0; // baseline - min (>=0)
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
        integral += (bl - wf[t]);
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
    double amplitude = bl - min_val;
    if (amplitude < 0.0) amplitude = 0.0;
    return {integral, peak_idx, amplitude};
}

static std::pair<int, std::unordered_map<std::string, std::string>>
discover_files_and_run(const std::string& input_dir, bool michel_mode, int run_filter) {
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

        if (std::find(CH_ORDER.begin(), CH_ORDER.end(), ch) == CH_ORDER.end()) continue;

        if (run_filter >= 0 && run != run_filter) {
            continue;
        }

        if (!run_set) {
            run_detected = run;
            run_set = true;
        } else if (run != run_detected) {
            throw std::runtime_error(
                "Multiple run numbers found in input-dir. Found run" + std::to_string(run_detected) +
                " and run" + std::to_string(run) + ". Use --run to select."
            );
        }

        if (ch_to_path.find(ch) == ch_to_path.end()) {
            ch_to_path[ch] = entry.path().string();
        }
    }

    if (!run_set) {
        throw std::runtime_error(
            "No matching files found in input-dir. Expected files like wave_PS_A_run12.txt, "
            "wave_NaI_A1_run12.txt, etc."
        );
    }

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
        else if (a == "--bins-per-event") cfg.bins_per_event = std::stoi(need(a));
        else if (a == "--baseline-bins") cfg.baseline_bins = std::stoi(need(a));
        else if (a == "--integ-start") cfg.integ_start = std::stoi(need(a));
        else if (a == "--integ-end") cfg.integ_end = std::stoi(need(a));
        else if (a == "--michel") cfg.michel_mode = true;
        else if (a == "--run") cfg.run_filter = std::stoi(need(a));
        else if (a == "--nai-a1-coeff") cfg.nai_a1_coeff = std::stod(need(a));
        else if (a == "--nai-a1-offset") cfg.nai_a1_offset = std::stod(need(a));
        else if (a == "--nai-a2-coeff") cfg.nai_a2_coeff = std::stod(need(a));
        else if (a == "--nai-a2-offset") cfg.nai_a2_offset = std::stod(need(a));
        else if (a == "--nai-b1-coeff") cfg.nai_b1_coeff = std::stod(need(a));
        else if (a == "--nai-b1-offset") cfg.nai_b1_offset = std::stod(need(a));
        else if (a == "--nai-b2-coeff") cfg.nai_b2_coeff = std::stod(need(a));
        else if (a == "--nai-b2-offset") cfg.nai_b2_offset = std::stod(need(a));
        else if (a == "--ps-amp-threshold") cfg.ps_amp_threshold = std::stod(need(a));
        else if (a == "--nai-integ-threshold") cfg.nai_integ_threshold = std::stod(need(a));
        else if (a == "--start-event") cfg.start_event = std::stoi(need(a));
        else if (a == "--max-events") cfg.max_events = std::stoi(need(a));
        else if (a == "-h" || a == "--help") { usage(argv[0]); std::exit(0); }
        else throw std::runtime_error("Unknown argument: " + a);
    }

    if (cfg.bins_per_event <= 0) throw std::runtime_error("--bins-per-event must be > 0");
    if (cfg.baseline_bins <= 0) throw std::runtime_error("--baseline-bins must be > 0");
    if (cfg.start_event < 0) throw std::runtime_error("--start-event must be >= 0");

    return cfg;
}

static double coeff_for_channel(const Config& cfg, const std::string& ch) {
    if (ch == "NaI_A1") return cfg.nai_a1_coeff;
    if (ch == "NaI_A2") return cfg.nai_a2_coeff;
    if (ch == "NaI_B1") return cfg.nai_b1_coeff;
    if (ch == "NaI_B2") return cfg.nai_b2_coeff;
    return 0.0;
}

static double offset_for_channel(const Config& cfg, const std::string& ch) {
    if (ch == "NaI_A1") return cfg.nai_a1_offset;
    if (ch == "NaI_A2") return cfg.nai_a2_offset;
    if (ch == "NaI_B1") return cfg.nai_b1_offset;
    if (ch == "NaI_B2") return cfg.nai_b2_offset;
    return 0.0;
}

static bool is_ps(const std::string& ch) {
    return ch.rfind("PS", 0) == 0;
}

static void gnuplot_event(FILE* gp,
                          int run,
                          int evt,
                          const std::unordered_map<std::string, std::vector<double>>& wf_map,
                          const std::unordered_map<std::string, double>& energy_map,
                          const std::unordered_map<std::string, double>& integral_map,
                          const std::unordered_map<std::string, bool>& ps_hit_map) {
    std::fprintf(gp, "reset\n");
    std::fprintf(gp, "set multiplot layout 3,2 rowsfirst title 'Run %d  Event %d'\n", run, evt);
    std::fprintf(gp, "set grid\n");

    for (const auto& ch : CH_ORDER) {
        std::string title;
        if (is_ps(ch)) {
            bool hit = ps_hit_map.at(ch);
            title = ch + std::string("  ") + (hit ? "HIT" : "NO");
        } else {
            double e = energy_map.at(ch);
            double integ = integral_map.at(ch);
            std::ostringstream oss;
            oss << ch << "  E=" << std::fixed << std::setprecision(1) << e
                << "  (I=" << std::fixed << std::setprecision(1) << integ << ")";
            title = oss.str();
        }

        std::fprintf(gp, "set title '%s'\n", title.c_str());
        std::fprintf(gp, "plot '-' using 1:2 with lines notitle\n");

        const auto& wf = wf_map.at(ch);
        for (size_t i = 0; i < wf.size(); ++i) {
            std::fprintf(gp, "%zu %f\n", i, wf[i]);
        }
        std::fprintf(gp, "e\n");
    }

    std::fprintf(gp, "unset multiplot\n");
    std::fflush(gp);
}

int main(int argc, char** argv) {
    try {
        Config cfg = parse_args(argc, argv);

        auto [run, ch_to_path] = discover_files_and_run(cfg.input_dir, cfg.michel_mode, cfg.run_filter);

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
                    std::to_string(nevt_ref) + "."
                );
            }
            ch_events.emplace(ch, std::move(events));
        }

        if (cfg.start_event >= static_cast<int>(nevt_ref)) {
            throw std::runtime_error("start_event is beyond total events: " + std::to_string(nevt_ref));
        }

        FILE* gp = popen("gnuplot -persist", "w");
        if (!gp) throw std::runtime_error("Failed to start gnuplot. Is it installed?");

        const int last_evt = (cfg.max_events < 0)
            ? static_cast<int>(nevt_ref)
            : std::min(static_cast<int>(nevt_ref), cfg.start_event + cfg.max_events);

        for (int evt = cfg.start_event; evt < last_evt; ++evt) {
            std::unordered_map<std::string, std::vector<double>> wf_map;
            std::unordered_map<std::string, double> integral_map;
            std::unordered_map<std::string, double> energy_map;
            std::unordered_map<std::string, bool> ps_hit_map;

            for (const auto& ch : CH_ORDER) {
                const auto& wf = ch_events.at(ch).at(static_cast<size_t>(evt));
                wf_map[ch] = wf;

                Features feat = compute_features(wf, cfg.baseline_bins, cfg.integ_start, cfg.integ_end);
                double integral = feat.integral;

                if (is_ps(ch)) {
                    bool hit = (feat.amplitude > cfg.ps_amp_threshold);
                    ps_hit_map[ch] = hit;
                    integral_map[ch] = integral;
                    energy_map[ch] = 0.0;
                } else {
                    if (integral <= cfg.nai_integ_threshold) {
                        integral = 0.0;
                    }
                    double energy = coeff_for_channel(cfg, ch) * integral + offset_for_channel(cfg, ch);
                    integral_map[ch] = integral;
                    energy_map[ch] = energy;
                }
            }

            gnuplot_event(gp, run, evt, wf_map, energy_map, integral_map, ps_hit_map);

            std::cout << "Event " << evt << " / " << nevt_ref - 1
                      << "  (Enter: next, q: quit): " << std::flush;
            std::string input;
            if (!std::getline(std::cin, input)) break;
            if (!input.empty() && (input[0] == 'q' || input[0] == 'Q')) break;
        }

        pclose(gp);
        return 0;

    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }
}
