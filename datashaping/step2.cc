// step2.cc
// Convert channel-based dataset into module-based physics dataset.
//
// Input (per line) from step1:
//   run  event  CH  integral  peak_time
//
// Output (per line, per module):
//   run  event  module  energy  peak_time_real  particle
//
// Rules:
// - module 1: PS_A, NaI_A1, NaI_A2
// - module 2: PS_B, NaI_B1, NaI_B2
// - energy = cA * NaI_A_integral + cB * NaI_B_integral
// - peak_time_real = min(NaI_A_peak, NaI_B_peak) * time_bin
// - particle = "positron" if PS_integral > 0 else "gamma"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

// =====================
// Config
// =====================
struct Config {
    double nai1_coeff = 1.0;   // coefficient for first NaI
    double nai2_coeff = 1.0;   // coefficient for second NaI
    double time_bin   = 4e-9;  // [s]
    std::string input_dir = "../data/shapeddata";
    std::string output_dir = "../data/shapeddata";
    std::string input_file;
    std::string output_file;
};

// =====================
// Data structures
// =====================
struct ChannelData {
    double integral = 0.0;
    int peak_time = -1;
};

struct EventData {
    int run = -1;
    int event = -1;
    std::map<std::string, ChannelData> ch;
};

// =====================
// Utilities
// =====================
static void usage(const char* prog) {
    std::cerr
        << "Usage:\n"
        << "  " << prog << " [options]\n\n"
        << "Input/Output:\n"
        << "  --input DIR/FILE        input file (default: auto-detect step1_runM.txt in --input-dir)\n"
        << "  --input-dir DIR         input directory (default: ../data/shapeddata)\n"
        << "  --output-dir DIR        output directory (default: ../data/shapeddata)\n"
        << "  --out FILE              output file (default: step2_runM.txt in --output-dir)\n\n"
        << "Options:\n"
        << "  --nai1-coeff C   coefficient for first NaI (default: 1.0)\n"
        << "  --nai2-coeff C   coefficient for second NaI (default: 1.0)\n"
        << "  --time-bin T     time per bin [s] (default: 4e-9)\n";
}

static std::pair<int, std::string> discover_step1_file(const std::string& input_dir) {
    std::regex re(R"(step1_run([0-9]+)\.txt$)");
    int run_detected = -1;
    bool run_set = false;
    std::string found_path;

    for (const auto& entry : fs::directory_iterator(input_dir)) {
        if (!entry.is_regular_file()) continue;
        const std::string fname = entry.path().filename().string();
        std::smatch m;
        if (!std::regex_match(fname, m, re)) continue;

        const int run = std::stoi(m[1].str());
        if (!run_set) {
            run_detected = run;
            run_set = true;
            found_path = entry.path().string();
        } else if (run != run_detected) {
            throw std::runtime_error(
                "Multiple step1_runM.txt files found in input-dir. Found run" +
                std::to_string(run_detected) + " and run" + std::to_string(run) +
                ". Please keep one run per directory or specify --input."
            );
        }
    }

    if (!run_set) {
        throw std::runtime_error(
            "No step1_runM.txt found in input-dir. Use --input or check --input-dir."
        );
    }

    return {run_detected, found_path};
}

static Config parse_args(int argc, char** argv) {
    Config cfg;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto need = [&](const std::string& opt) -> std::string {
            if (i + 1 >= argc) throw std::runtime_error("Missing value for " + opt);
            return std::string(argv[++i]);
        };

        if (a == "--input") cfg.input_file = need(a);
        else if (a == "--input-dir") cfg.input_dir = need(a);
        else if (a == "--output-dir") cfg.output_dir = need(a);
        else if (a == "--out") cfg.output_file = need(a);
        else if (a == "--nai1-coeff") cfg.nai1_coeff = std::stod(need(a));
        else if (a == "--nai2-coeff") cfg.nai2_coeff = std::stod(need(a));
        else if (a == "--time-bin") cfg.time_bin = std::stod(need(a));
        else if (a == "-h" || a == "--help") { usage(argv[0]); std::exit(0); }
        else if (!a.empty() && a[0] != '-') {
            if (cfg.input_file.empty()) cfg.input_file = a;
            else if (cfg.output_file.empty()) cfg.output_file = a;
            else throw std::runtime_error("Too many positional arguments: " + a);
        } else {
            throw std::runtime_error("Unknown option: " + a);
        }
    }

    return cfg;
}

// =====================
// Main
// =====================
int main(int argc, char** argv) {
    try {
        Config cfg = parse_args(argc, argv);

        int run = -1;
        if (cfg.input_file.empty()) {
            auto [detected_run, path] = discover_step1_file(cfg.input_dir);
            run = detected_run;
            cfg.input_file = path;
        }
        if (cfg.output_file.empty()) {
            std::regex re(R"(step1_run([0-9]+)\.txt$)");
            std::smatch m;
            if (std::regex_search(cfg.input_file, m, re)) {
                run = std::stoi(m[1].str());
            } else if (run < 0) {
                auto [detected_run, _path] = discover_step1_file(cfg.input_dir);
                run = detected_run;
            }
            if (run < 0) {
                throw std::runtime_error("Cannot infer run number from input file. Specify --out.");
            }
            cfg.output_file = (fs::path(cfg.output_dir) /
                               ("step2_run" + std::to_string(run) + ".txt")).string();
        }

        std::ifstream ifs(cfg.input_file);
        if (!ifs) throw std::runtime_error("Cannot open input file: " + cfg.input_file);

        // Group by (run, event)
        std::map<std::pair<int,int>, EventData> events;

        std::string line;
        while (std::getline(ifs, line)) {
            if (line.empty() || line[0] == '#') continue;

            std::istringstream iss(line);
            int run, event, peak;
            std::string ch;
            double integral;

            if (!(iss >> run >> event >> ch >> integral >> peak)) {
                // skip header or malformed line
                continue;
            }

            auto& ev = events[{run, event}];
            ev.run = run;
            ev.event = event;
            ev.ch[ch] = {integral, peak};
        }

        std::ofstream ofs(cfg.output_file);
        if (!ofs) throw std::runtime_error("Cannot open output file: " + cfg.output_file);

        ofs << "run\tevent\tmodule\tenergy\tpeak_time_real\tparticle\n";

        // Process events
        for (const auto& [key, ev] : events) {
            // -------- module 1 --------
            {
                const auto& ps   = ev.ch.at("PS_A");
                const auto& nai1 = ev.ch.at("NaI_A1");
                const auto& nai2 = ev.ch.at("NaI_A2");

                double energy =
                    cfg.nai1_coeff * nai1.integral +
                    cfg.nai2_coeff * nai2.integral;

                int peak_bin = std::min(nai1.peak_time, nai2.peak_time);
                double peak_time_real = peak_bin * cfg.time_bin;

                std::string particle =
                    (ps.integral > 0.0) ? "positron" : "gamma";

                ofs << ev.run << "\t"
                    << ev.event << "\t"
                    << 1 << "\t"
                    << std::fixed << std::setprecision(6) << energy << "\t"
                    << std::scientific << std::setprecision(2) << peak_time_real << "\t"
                    << particle << "\n";
            }

            // -------- module 2 --------
            {
                const auto& ps   = ev.ch.at("PS_B");
                const auto& nai1 = ev.ch.at("NaI_B1");
                const auto& nai2 = ev.ch.at("NaI_B2");

                double energy =
                    cfg.nai1_coeff * nai1.integral +
                    cfg.nai2_coeff * nai2.integral;

                int peak_bin = std::min(nai1.peak_time, nai2.peak_time);
                double peak_time_real = peak_bin * cfg.time_bin;

                std::string particle =
                    (ps.integral > 0.0) ? "positron" : "gamma";

                ofs << ev.run << "\t"
                    << ev.event << "\t"
                    << 2 << "\t"
                    << std::fixed << std::setprecision(6) << energy << "\t"
                    << std::scientific << std::setprecision(2) << peak_time_real << "\t"
                    << particle << "\n";
            }
        }

        std::cout << "Processed events: " << events.size() << "\n";
        std::cout << "Output written to: " << cfg.output_file << "\n";
        return 0;

    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }
}
