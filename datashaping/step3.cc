// step3.cc
// Read module-based dataset and split into 3 files by (module1 particle, module2 particle) combination
// within the same (run, event).
//
// Input per line (tab/space separated, header allowed):
//   run  event  module  energy  peak_time_real  particle
//
// Output: 3 files
//   out_prefix_pp.txt  : (positron, positron)
//   out_prefix_pg.txt  : (positron, gamma) OR (gamma, positron)  <-- treated as the same class
//   out_prefix_gg.txt  : (gamma, gamma)
//
// Each output file contains the ORIGINAL TWO LINES (module 1 and module 2) for that (run, event),
// i.e., it outputs the same columns:
//   run  event  module  energy  peak_time_real  particle
//
// Notes / assumptions:
// - For each (run,event), there should be exactly two lines: module 1 and module 2.
// - If duplicates or missing modules are found, the event is skipped with a warning to stderr.
//
// Build:
//   g++ -O2 -std=c++17 step3.cc -o step3
//
// Run:
//   ./step3
//   ./step3 input.txt out_prefix
//   -> writes step3_runM_pp.txt, step3_runM_pg.txt, step3_runM_gg.txt

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

namespace fs = std::filesystem;

struct Config {
    std::string input_dir = "../data/shapeddata";
    std::string output_dir = "../data/shapeddata";
    std::string input_file;
    std::string out_prefix;
};

struct LineRec {
    int run = -1;
    int event = -1;
    int module = -1;            // 1 or 2
    double energy = 0.0;
    double peak_time_real = 0.0;
    std::string particle;       // "positron" or "gamma"
    std::string raw_line;       // normalized output line (tab-separated)
};

static void usage(const char* prog) {
    std::cerr
        << "Usage:\n"
        << "  " << prog << " [options]\n\n"
        << "Input/Output:\n"
        << "  --input DIR/FILE        input file (default: auto-detect step2_runM.txt in --input-dir)\n"
        << "  --input-dir DIR         input directory (default: ../data/shapeddata)\n"
        << "  --output-dir DIR        output directory (default: ../data/shapeddata)\n"
        << "  --out-prefix PREFIX     output prefix (default: step3_runM in --output-dir)\n\n"
        << "Outputs:\n"
        << "  out_prefix_pp.txt  (positron,positron)\n"
        << "  out_prefix_pg.txt  (positron,gamma) or (gamma,positron)\n"
        << "  out_prefix_gg.txt  (gamma,gamma)\n";
}

static std::pair<int, std::string> discover_step2_file(const std::string& input_dir) {
    std::regex re(R"(step2_run([0-9]+)\.txt$)");
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
                "Multiple step2_runM.txt files found in input-dir. Found run" +
                std::to_string(run_detected) + " and run" + std::to_string(run) +
                ". Please keep one run per directory or specify --input."
            );
        }
    }

    if (!run_set) {
        throw std::runtime_error(
            "No step2_runM.txt found in input-dir. Use --input or check --input-dir."
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
        else if (a == "--out-prefix") cfg.out_prefix = need(a);
        else if (a == "-h" || a == "--help") { usage(argv[0]); std::exit(0); }
        else if (!a.empty() && a[0] != '-') {
            if (cfg.input_file.empty()) cfg.input_file = a;
            else if (cfg.out_prefix.empty()) cfg.out_prefix = a;
            else throw std::runtime_error("Too many positional arguments: " + a);
        } else {
            throw std::runtime_error("Unknown option: " + a);
        }
    }

    return cfg;
}

static bool parse_line(const std::string& line, LineRec& rec) {
    if (line.empty()) return false;
    if (line[0] == '#') return false;

    std::istringstream iss(line);

    int run, event, module;
    double energy, t;
    std::string particle;

    if (!(iss >> run >> event >> module >> energy >> t >> particle)) {
        return false; // header or malformed
    }

    // Normalize particle string (just in case)
    if (particle != "positron" && particle != "gamma") {
        return false;
    }

    rec.run = run;
    rec.event = event;
    rec.module = module;
    rec.energy = energy;
    rec.peak_time_real = t;
    rec.particle = particle;

    // Create a normalized, tab-separated line for output
    std::ostringstream oss;
    oss << rec.run << "\t"
        << rec.event << "\t"
        << rec.module << "\t"
        << std::fixed << std::setprecision(6) << rec.energy << "\t"
        << std::scientific << std::setprecision(2) << rec.peak_time_real << "\t"
        << rec.particle;
    rec.raw_line = oss.str();

    return true;
}

struct EventPair {
    bool has1 = false;
    bool has2 = false;
    LineRec m1;
    LineRec m2;
};

static std::string classify_pair(const EventPair& ep) {
    // Determine combination using module1 and module2 particles
    // Return: "pp", "pg", "gg"
    const std::string& a = ep.m1.particle;
    const std::string& b = ep.m2.particle;

    if (a == "positron" && b == "positron") return "pp";
    if (a == "gamma"    && b == "gamma")    return "gg";
    // mixed
    return "pg";
}

int main(int argc, char** argv) {
    try {
        Config cfg = parse_args(argc, argv);

        int run = -1;
        if (cfg.input_file.empty()) {
            auto [detected_run, path] = discover_step2_file(cfg.input_dir);
            run = detected_run;
            cfg.input_file = path;
        }
        if (cfg.out_prefix.empty()) {
            std::regex re(R"(step2_run([0-9]+)\.txt$)");
            std::smatch m;
            if (std::regex_search(cfg.input_file, m, re)) {
                run = std::stoi(m[1].str());
            } else if (run < 0) {
                auto [detected_run, _path] = discover_step2_file(cfg.input_dir);
                run = detected_run;
            }
            if (run < 0) {
                throw std::runtime_error("Cannot infer run number from input file. Specify --out-prefix.");
            }
            cfg.out_prefix = (fs::path(cfg.output_dir) /
                              ("step3_run" + std::to_string(run))).string();
        }

        std::ifstream ifs(cfg.input_file);
        if (!ifs) throw std::runtime_error("Cannot open input file: " + cfg.input_file);

        std::ofstream ofs_pp(cfg.out_prefix + "_pp.txt");
        std::ofstream ofs_pg(cfg.out_prefix + "_pg.txt");
        std::ofstream ofs_gg(cfg.out_prefix + "_gg.txt");
        if (!ofs_pp || !ofs_pg || !ofs_gg) {
            throw std::runtime_error("Cannot open one of output files with prefix: " + cfg.out_prefix);
        }

        // Headers
        const std::string header = "run\tevent\tmodule\tenergy\tpeak_time_real\tparticle\n";
        ofs_pp << header;
        ofs_pg << header;
        ofs_gg << header;

        // Collect per (run,event)
        std::map<std::pair<int,int>, EventPair> mp;

        std::string line;
        long long line_no = 0;
        while (std::getline(ifs, line)) {
            ++line_no;
            LineRec rec;
            if (!parse_line(line, rec)) continue;

            if (rec.module != 1 && rec.module != 2) {
                std::cerr << "WARN: line " << line_no << " has unexpected module=" << rec.module
                          << " (skipped)\n";
                continue;
            }

            auto& ep = mp[{rec.run, rec.event}];

            if (rec.module == 1) {
                if (ep.has1) {
                    std::cerr << "WARN: duplicate module 1 for run=" << rec.run
                              << " event=" << rec.event << " (keeping first, skipping this line)\n";
                    continue;
                }
                ep.has1 = true;
                ep.m1 = rec;
            } else {
                if (ep.has2) {
                    std::cerr << "WARN: duplicate module 2 for run=" << rec.run
                              << " event=" << rec.event << " (keeping first, skipping this line)\n";
                    continue;
                }
                ep.has2 = true;
                ep.m2 = rec;
            }
        }

        // Write out by class
        size_t written_pp = 0, written_pg = 0, written_gg = 0, skipped = 0;

        for (const auto& [key, ep] : mp) {
            if (!ep.has1 || !ep.has2) {
                ++skipped;
                std::cerr << "WARN: missing module for run=" << key.first << " event=" << key.second
                          << " (skipped)\n";
                continue;
            }

            const std::string cls = classify_pair(ep);

            if (cls == "pp") {
                ofs_pp << ep.m1.raw_line << "\n";
                ofs_pp << ep.m2.raw_line << "\n";
                ++written_pp;
            } else if (cls == "gg") {
                ofs_gg << ep.m1.raw_line << "\n";
                ofs_gg << ep.m2.raw_line << "\n";
                ++written_gg;
            } else { // "pg"
                ofs_pg << ep.m1.raw_line << "\n";
                ofs_pg << ep.m2.raw_line << "\n";
                ++written_pg;
            }
        }

        std::cout << "Done.\n"
                  << "Events written: pp=" << written_pp
                  << ", pg=" << written_pg
                  << ", gg=" << written_gg
                  << "\nSkipped (incomplete): " << skipped << "\n"
                  << "Output files:\n"
                  << "  " << cfg.out_prefix << "_pp.txt\n"
                  << "  " << cfg.out_prefix << "_pg.txt\n"
                  << "  " << cfg.out_prefix << "_gg.txt\n";

        return 0;

    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }
}
