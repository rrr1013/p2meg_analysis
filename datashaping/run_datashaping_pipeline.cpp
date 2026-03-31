#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

// ============================================================
// run_datashaping_pipeline.cpp
//
// 目的:
//  - data/rawdata/runN のようなフォルダを入力し、最終出力まで一括実行する
//  - 複数フォルダの同時指定に対応
//  - runごとに最終出力を分離
// ============================================================

namespace fs = std::filesystem;

struct Config {
    std::string out_base = "data/shapeddata";
    std::string bin_dir;
    int bins_per_event = 1000;
    bool cleanup_intermediate = false;
    bool debug = false;
    bool dry_run = false;
    bool continue_on_error = false;
    std::vector<std::string> run_dirs;
};

static std::string Quote(const std::string& s) {
    std::string out = "'";
    for (char c : s) {
        if (c == '\'') out += "'\\''";
        else out.push_back(c);
    }
    out.push_back('\'');
    return out;
}

static int RunCommand(const std::string& cmd, bool dry_run) {
    std::cout << "[cmd] " << cmd << "\n";
    if (dry_run) return 0;
    return std::system(cmd.c_str());
}

static void EnsureDir(const fs::path& p) {
    std::error_code ec;
    fs::create_directories(p, ec);
    if (ec) {
        throw std::runtime_error("Failed to create directory: " + p.string() +
                                 " (" + ec.message() + ")");
    }
}

static std::pair<int, std::map<std::string, std::string>> DiscoverWaveFiles(const std::string& dir) {
    std::map<std::string, std::string> ch_to_path;
    std::set<int> runs;

    const std::vector<std::string> channels = {
        "NaI_A1", "NaI_A2", "NaI_B1", "NaI_B2", "PS_A", "PS_B"
    };

    for (const auto& entry : fs::directory_iterator(dir)) {
        if (!entry.is_regular_file()) continue;

        const std::string name = entry.path().filename().string();
        if (name.rfind("wave_", 0) != 0) continue;
        if (name.size() < 13) continue;

        // wave_<CH>_run<NNN>.txt
        const size_t p_run = name.rfind("_run");
        if (p_run == std::string::npos) continue;
        const size_t p_dot = name.rfind(".txt");
        if (p_dot == std::string::npos || p_dot <= p_run + 4) continue;

        const std::string ch = name.substr(5, p_run - 5);
        bool known = false;
        for (const auto& c : channels) {
            if (c == ch) {
                known = true;
                break;
            }
        }
        if (!known) continue;

        const std::string run_s = name.substr(p_run + 4, p_dot - (p_run + 4));
        if (run_s.empty()) continue;
        int run = std::atoi(run_s.c_str());
        if (run <= 0) continue;

        runs.insert(run);
        if (ch_to_path.find(ch) == ch_to_path.end()) {
            ch_to_path[ch] = entry.path().string();
        }
    }

    if (runs.empty()) {
        throw std::runtime_error("No wave_*.txt found in directory: " + dir);
    }
    if (runs.size() != 1) {
        throw std::runtime_error("Multiple run numbers found in directory: " + dir);
    }

    for (const auto& ch : channels) {
        if (ch_to_path.find(ch) == ch_to_path.end()) {
            throw std::runtime_error("Missing channel file " + ch + " in directory: " + dir);
        }
    }

    return {*runs.begin(), ch_to_path};
}

static bool ExecuteOrThrow(const std::string& cmd,
                           bool dry_run,
                           bool continue_on_error,
                           const std::string& context) {
    const int rc = RunCommand(cmd, dry_run);
    if (rc == 0) return true;

    const std::string msg = "Command failed (" + context + "), rc=" + std::to_string(rc);
    if (continue_on_error) {
        std::cerr << "warn: " << msg << "\n";
        return false;
    }
    throw std::runtime_error(msg);
}

static void CleanupIntermediate(const std::vector<std::string>& files) {
    for (const auto& f : files) {
        std::error_code ec;
        fs::remove(f, ec);
    }
}

static Config ParseArgs(int argc, char** argv) {
    Config cfg;

    fs::path exe_path(argv[0]);
    if (!exe_path.parent_path().empty()) {
        cfg.bin_dir = exe_path.parent_path().string();
    } else {
        cfg.bin_dir = ".";
    }

    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto need = [&](const std::string& opt) -> std::string {
            if (i + 1 >= argc) throw std::runtime_error("Missing value for " + opt);
            return std::string(argv[++i]);
        };

        if (a == "--out-base") cfg.out_base = need(a);
        else if (a == "--bin-dir") cfg.bin_dir = need(a);
        else if (a == "--bins-per-event") cfg.bins_per_event = std::atoi(need(a).c_str());
        else if (a == "--cleanup-intermediate") cfg.cleanup_intermediate = true;
        else if (a == "--debug") cfg.debug = true;
        else if (a == "--dry-run") cfg.dry_run = true;
        else if (a == "--continue-on-error") cfg.continue_on_error = true;
        else if (a == "-h" || a == "--help") {
            std::cout
                << "Usage:\n"
                << "  " << argv[0] << " [options] <run_dir1> [run_dir2 ...]\n"
                << "Options:\n"
                << "  --out-base DIR             output base directory (default: data/shapeddata)\n"
                << "  --bin-dir DIR              tool binary directory (default: this executable dir)\n"
                << "  --bins-per-event N         bins per event (default: 1000)\n"
                << "  --cleanup-intermediate     remove intermediate files after each run\n"
                << "  --debug                    enable debug mode where supported\n"
                << "  --dry-run                  print commands only\n"
                << "  --continue-on-error        continue to next run on failure\n";
            std::exit(0);
        } else if (!a.empty() && a[0] == '-') {
            throw std::runtime_error("Unknown option: " + a);
        } else {
            cfg.run_dirs.push_back(a);
        }
    }

    if (cfg.run_dirs.empty()) {
        throw std::runtime_error("No run directories specified.");
    }
    if (cfg.bins_per_event <= 0) {
        throw std::runtime_error("--bins-per-event must be > 0");
    }

    return cfg;
}

int main(int argc, char** argv) {
    try {
        const Config cfg = ParseArgs(argc, argv);

        const std::string exe_wave = (fs::path(cfg.bin_dir) / "wave_to_pileup_csv").string();
        const std::string exe_base = (fs::path(cfg.bin_dir) / "wave_to_baseline_template").string();
        const std::string exe_pile = (fs::path(cfg.bin_dir) / "run_pileup_residual_standalone").string();
        const std::string exe_ps = (fs::path(cfg.bin_dir) / "run_pssignal_standalone").string();
        const std::string exe_merge = (fs::path(cfg.bin_dir) / "build_module_summary").string();

        std::set<int> processed_runs;

        for (size_t idir = 0; idir < cfg.run_dirs.size(); ++idir) {
            const std::string run_dir_in = cfg.run_dirs[idir];
            std::cout << "[run_dir] " << run_dir_in << "\n";

            bool run_ok = true;
            try {
                const auto discovered = DiscoverWaveFiles(run_dir_in);
                const int run = discovered.first;
                const std::map<std::string, std::string>& files = discovered.second;

                if (processed_runs.find(run) != processed_runs.end()) {
                    throw std::runtime_error("Duplicate run number across inputs: run" + std::to_string(run));
                }
                processed_runs.insert(run);

                const fs::path out_run = fs::path(cfg.out_base) / ("run" + std::to_string(run));
                EnsureDir(out_run);

                std::vector<std::string> intermediates;

                const std::vector<std::string> nai_ch = {"NaI_A1", "NaI_A2", "NaI_B1", "NaI_B2"};
                const std::vector<std::string> ps_ch = {"PS_A", "PS_B"};

                for (const auto& ch : nai_ch) {
                    const fs::path out_ch = out_run / ch;
                    EnsureDir(out_ch);

                    const std::string in_wave = files.at(ch);
                    const std::string out_pileup = (out_ch / "pileup_waveforms.csv").string();
                    const std::string out_base_csv = (out_ch / "event_baseline_noise.csv").string();
                    const std::string out_template = (out_ch / "template.json").string();
                    const std::string out_pulses = (out_ch / "pulses.csv").string();

                    std::string cmd;

                    cmd = Quote(exe_wave) + " " + Quote(in_wave) + " " + Quote(out_pileup) +
                          " --bins-per-event " + std::to_string(cfg.bins_per_event);
                    if (!ExecuteOrThrow(cmd, cfg.dry_run, cfg.continue_on_error, ch + ": wave_to_pileup_csv")) {
                        run_ok = false;
                        break;
                    }

                    cmd = Quote(exe_base) + " " + Quote(in_wave) +
                          " --out-baseline " + Quote(out_base_csv) +
                          " --out-template " + Quote(out_template) +
                          " --bins-per-event " + std::to_string(cfg.bins_per_event);
                    if (cfg.debug) cmd += " --debug";
                    if (!ExecuteOrThrow(cmd, cfg.dry_run, cfg.continue_on_error, ch + ": wave_to_baseline_template")) {
                        run_ok = false;
                        break;
                    }

                    cmd = Quote(exe_pile) + " " + Quote(out_template) + " " + Quote(out_base_csv) +
                          " " + Quote(out_pileup) + " --out " + Quote(out_pulses);
                    if (cfg.debug) cmd += " --debug";
                    if (!ExecuteOrThrow(cmd, cfg.dry_run, cfg.continue_on_error, ch + ": pileup_decompose")) {
                        run_ok = false;
                        break;
                    }

                    intermediates.push_back(out_pileup);
                    intermediates.push_back(out_base_csv);
                    intermediates.push_back(out_template);
                    intermediates.push_back(out_pulses);
                }
                if (!run_ok) continue;

                for (const auto& ch : ps_ch) {
                    const fs::path out_ch = out_run / ch;
                    EnsureDir(out_ch);

                    const std::string in_wave = files.at(ch);
                    const std::string out_pileup = (out_ch / "pileup_waveforms.csv").string();
                    const std::string out_base_csv = (out_ch / "event_baseline_noise.csv").string();
                    const std::string out_pssignal = (out_ch / "PSsignal.csv").string();

                    std::string cmd;

                    cmd = Quote(exe_wave) + " " + Quote(in_wave) + " " + Quote(out_pileup) +
                          " --bins-per-event " + std::to_string(cfg.bins_per_event);
                    if (!ExecuteOrThrow(cmd, cfg.dry_run, cfg.continue_on_error, ch + ": wave_to_pileup_csv")) {
                        run_ok = false;
                        break;
                    }

                    cmd = Quote(exe_base) + " " + Quote(in_wave) +
                          " --out-baseline " + Quote(out_base_csv) +
                          " --bins-per-event " + std::to_string(cfg.bins_per_event) +
                          " --baseline-only";
                    if (cfg.debug) cmd += " --debug";
                    if (!ExecuteOrThrow(cmd, cfg.dry_run, cfg.continue_on_error, ch + ": baseline_only")) {
                        run_ok = false;
                        break;
                    }

                    cmd = Quote(exe_ps) + " " + Quote(out_base_csv) + " " + Quote(out_pileup) +
                          " " + Quote(out_pssignal);
                    if (!ExecuteOrThrow(cmd, cfg.dry_run, cfg.continue_on_error, ch + ": PSsignal")) {
                        run_ok = false;
                        break;
                    }

                    intermediates.push_back(out_pileup);
                    intermediates.push_back(out_base_csv);
                    intermediates.push_back(out_pssignal);
                }
                if (!run_ok) continue;

                const std::string out_final =
                    (out_run / ("final_module_pulses_run" + std::to_string(run) + ".txt")).string();
                const std::string out_unmatched =
                    (out_run / ("unmatched_ps_run" + std::to_string(run) + ".csv")).string();

                std::string cmd = Quote(exe_merge) +
                                  " --nai-a1 " + Quote((out_run / "NaI_A1" / "pulses.csv").string()) +
                                  " --nai-a2 " + Quote((out_run / "NaI_A2" / "pulses.csv").string()) +
                                  " --nai-b1 " + Quote((out_run / "NaI_B1" / "pulses.csv").string()) +
                                  " --nai-b2 " + Quote((out_run / "NaI_B2" / "pulses.csv").string()) +
                                  " --ps-a " + Quote((out_run / "PS_A" / "PSsignal.csv").string()) +
                                  " --ps-b " + Quote((out_run / "PS_B" / "PSsignal.csv").string()) +
                                  " --out-main " + Quote(out_final) +
                                  " --out-unmatched " + Quote(out_unmatched);
                if (cfg.debug) cmd += " --debug";

                if (!ExecuteOrThrow(cmd, cfg.dry_run, cfg.continue_on_error, "build_module_summary")) {
                    run_ok = false;
                }

                if (run_ok && cfg.cleanup_intermediate) {
                    CleanupIntermediate(intermediates);
                }

                if (run_ok) {
                    std::cout << "[done] run" << run
                              << " final=" << out_final
                              << " unmatched=" << out_unmatched << "\n";
                }
            } catch (const std::exception& e) {
                run_ok = false;
                if (cfg.continue_on_error) {
                    std::cerr << "warn: " << e.what() << "\n";
                } else {
                    throw;
                }
            }

            if (!run_ok && !cfg.continue_on_error) {
                return 1;
            }
        }

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
}
