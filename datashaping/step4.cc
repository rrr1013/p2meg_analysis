// step4.cc
// From (positron,gamma) class file (module-based dataset), produce:
//   (E_e, E_g, t, phi_e, phi_g)
//
// Input per line (tab/space separated, header allowed):
//   run  event  module  energy  peak_time_real  particle
//
// Assumptions:
// - The input file contains ONLY events that are (positron,gamma) pairs,
//   and therefore for each (run,event) there should be exactly two lines:
//   one with particle=positron and one with particle=gamma.
// - t is taken as (positron peak_time_real) - (gamma peak_time_real).
//
// Phi rules (UPDATED per request):
// - One of (phi_e, phi_g) is 0. The other is a user-provided constant angle.
// - If module 1 is positron:
//     phi_e = 0
//     phi_g = PHI_NONZERO   (provided by user)
// - If module 1 is gamma:
//     phi_g = 0
//     phi_e = PHI_NONZERO   (provided by user)
//
// Output per event (one line):
//   E_e  E_g  t  phi_e  phi_g
//
// Build:
//   g++ -O2 -std=c++17 step4.cc -o step4
//
// Run:
//   ./step4 --phi-nonzero 0.34906585
//   ./step4 pg_input.txt eg_features.txt --phi-nonzero 0.34906585
//   (example phi=20 deg in rad)

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
    double phi_nonzero = 0.0; // fallback (e.g., radians)
    bool   phi_set = false;
    std::string input_dir = "../data/shapeddata";
    std::string output_dir = "../data/shapeddata";
    std::string input_file;
    std::string output_file;
};

static const std::map<int, double> RUN_PHI = {
    // Defaults: run1..run9 -> 20,40,60,...,180 deg (radians)
    {1, 0.34906585},  // 20 deg
    {2, 0.69813170},  // 40 deg
    {3, 1.04719755},  // 60 deg
    {4, 1.39626340},  // 80 deg
    {5, 1.74532925},  // 100 deg
    {6, 2.09439510},  // 120 deg
    {7, 2.44346095},  // 140 deg
    {8, 2.79252680},  // 160 deg
    {9, 3.14159265},  // 180 deg
};

struct Rec {
    int run = -1;
    int event = -1;
    int module = -1;             // 1 or 2
    double energy = 0.0;
    double t = 0.0;              // peak_time_real
    std::string particle;        // "positron" or "gamma"
};

static void usage(const char* prog) {
    std::cerr
        << "Usage:\n"
        << "  " << prog << " [options] [--phi-nonzero PHI]\n\n"
        << "Input/Output:\n"
        << "  --input DIR/FILE        input file (default: auto-detect step3_runM_pg.txt in --input-dir)\n"
        << "  --input-dir DIR         input directory (default: ../data/shapeddata)\n"
        << "  --output-dir DIR        output directory (default: ../data/shapeddata)\n"
        << "  --out FILE              output file (default: step4_runM_eg.txt in --output-dir)\n\n"
        << "Input must be (positron,gamma) events with two lines per (run,event).\n"
        << "PHI is the constant angle assigned to the non-zero side (recommend radians).\n"
        << "If RUN_PHI is set in code, it overrides --phi-nonzero per run.\n";
}

static std::vector<std::pair<int, std::string>> discover_step3_pg_files(const std::string& input_dir) {
    std::regex re(R"(step3_run([0-9]+)_pg\.txt$)");
    std::map<int, std::string> run_to_path;

    for (const auto& entry : fs::directory_iterator(input_dir)) {
        if (!entry.is_regular_file()) continue;
        const std::string fname = entry.path().filename().string();
        std::smatch m;
        if (!std::regex_match(fname, m, re)) continue;

        const int run = std::stoi(m[1].str());
        if (run_to_path.find(run) == run_to_path.end()) {
            run_to_path[run] = entry.path().string();
        }
    }

    if (run_to_path.empty()) {
        throw std::runtime_error(
            "No step3_runM_pg.txt found in input-dir. Use --input or check --input-dir."
        );
    }

    std::vector<std::pair<int, std::string>> out;
    out.reserve(run_to_path.size());
    for (const auto& [run, path] : run_to_path) {
        out.emplace_back(run, path);
    }
    return out;
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
        else if (a == "--phi-nonzero") {
            cfg.phi_nonzero = std::stod(need(a));
            cfg.phi_set = true;
        } else if (a == "-h" || a == "--help") {
            usage(argv[0]);
            std::exit(0);
        } else if (!a.empty() && a[0] != '-') {
            if (cfg.input_file.empty()) cfg.input_file = a;
            else if (cfg.output_file.empty()) cfg.output_file = a;
            else throw std::runtime_error("Too many positional arguments: " + a);
        } else {
            throw std::runtime_error("Unknown option: " + a);
        }
    }

    return cfg;
}

static bool parse_line(const std::string& line, Rec& r) {
    if (line.empty()) return false;
    if (line[0] == '#') return false;

    std::istringstream iss(line);
    int run, event, module;
    double energy, t;
    std::string particle;

    if (!(iss >> run >> event >> module >> energy >> t >> particle)) {
        return false; // header or malformed
    }
    if (particle != "positron" && particle != "gamma") return false;
    if (module != 1 && module != 2) return false;

    r.run = run;
    r.event = event;
    r.module = module;
    r.energy = energy;
    r.t = t;
    r.particle = particle;
    return true;
}

struct Pair {
    bool has_pos = false;
    bool has_gam = false;
    Rec pos;
    Rec gam;
};

static void process_file(const std::string& input_file,
                         std::ofstream& ofs,
                         double phi_nonzero,
                         int expected_run) {
    std::ifstream ifs(input_file);
    if (!ifs) throw std::runtime_error("Cannot open input file: " + input_file);

    // group by (run,event)
    std::map<std::pair<int,int>, Pair> mp;

    std::string line;
    long long line_no = 0;
    while (std::getline(ifs, line)) {
        ++line_no;
        Rec r;
        if (!parse_line(line, r)) continue;

        auto& p = mp[{r.run, r.event}];

        if (r.particle == "positron") {
            if (p.has_pos) {
                std::cerr << "WARN: duplicate positron line for run=" << r.run
                          << " event=" << r.event << " (skipping line " << line_no << ")\n";
                continue;
            }
            p.has_pos = true;
            p.pos = r;
        } else { // gamma
            if (p.has_gam) {
                std::cerr << "WARN: duplicate gamma line for run=" << r.run
                          << " event=" << r.event << " (skipping line " << line_no << ")\n";
                continue;
            }
            p.has_gam = true;
            p.gam = r;
        }
    }

    size_t written = 0, skipped = 0;

    for (const auto& [key, p] : mp) {
        if (!p.has_pos || !p.has_gam) {
            ++skipped;
            std::cerr << "WARN: incomplete pair for run=" << key.first
                      << " event=" << key.second << " (skipped)\n";
            continue;
        }

        const double E_e = p.pos.energy;
        const double E_g = p.gam.energy;

            // t: positron peak_time_real minus gamma peak_time_real
            const double t = p.pos.t - p.gam.t;

        // Determine particle on module 1
        std::string m1_particle;
        if (p.pos.module == 1) m1_particle = "positron";
        else if (p.gam.module == 1) m1_particle = "gamma";
        else {
            ++skipped;
            std::cerr << "WARN: no module 1 line for run=" << key.first
                      << " event=" << key.second << " (skipped)\n";
            continue;
        }

        double phi_e = 0.0;
        double phi_g = 0.0;

        if (m1_particle == "positron") {
            // module1=positron => phi_e=0, phi_g=given
            phi_e = 0.0;
            phi_g = phi_nonzero;
        } else {
            // module1=gamma => phi_g=0, phi_e=given
            phi_g = 0.0;
            phi_e = phi_nonzero;
        }

        ofs << std::fixed << std::setprecision(6) << E_e << "\t"
            << std::fixed << std::setprecision(6) << E_g << "\t"
            << std::scientific << std::setprecision(2) << t << "\t"
            << std::fixed << std::setprecision(6) << phi_e << "\t"
            << std::fixed << std::setprecision(6) << phi_g << "\n";
        ++written;
    }

    std::cout << "Processed: " << input_file << "\n"
              << "Events written: " << written << "\n"
              << "Events skipped (incomplete/invalid): " << skipped << "\n";
}

int main(int argc, char** argv) {
    try {
        Config cfg = parse_args(argc, argv);

        if (cfg.input_file.empty()) {
            const auto files = discover_step3_pg_files(cfg.input_dir);
            if (cfg.output_file.empty()) {
                cfg.output_file = (fs::path(cfg.output_dir) / "step4_eg_all.txt").string();
            }
            std::ofstream ofs(cfg.output_file);
            if (!ofs) throw std::runtime_error("Cannot open output file: " + cfg.output_file);
            ofs << "E_e\tE_g\tt\tphi_e\tphi_g\n";
            ofs << std::fixed << std::setprecision(6);
            for (const auto& [run, path] : files) {
                auto it = RUN_PHI.find(run);
                if (it == RUN_PHI.end() && !cfg.phi_set) {
                    throw std::runtime_error(
                        "No phi defined for run " + std::to_string(run) +
                        ". Set RUN_PHI in code or pass --phi-nonzero."
                    );
                }
                const double phi = (it != RUN_PHI.end()) ? it->second : cfg.phi_nonzero;
                process_file(path, ofs, phi, run);
            }
            std::cout << "Output: " << cfg.output_file << "\n";
        } else {
            if (cfg.output_file.empty()) {
                std::regex re(R"(step3_run([0-9]+)_pg\.txt$)");
                std::smatch m;
                int run = -1;
                if (std::regex_search(cfg.input_file, m, re)) {
                    run = std::stoi(m[1].str());
                }
                if (run < 0) {
                    throw std::runtime_error("Cannot infer run number from input file. Specify --out.");
                }
                cfg.output_file = (fs::path(cfg.output_dir) /
                                   ("step4_run" + std::to_string(run) + "_eg.txt")).string();
            }
            std::ofstream ofs(cfg.output_file);
            if (!ofs) throw std::runtime_error("Cannot open output file: " + cfg.output_file);
            ofs << "E_e\tE_g\tt\tphi_e\tphi_g\n";
            ofs << std::fixed << std::setprecision(6);
            int run = -1;
            std::regex re(R"(step3_run([0-9]+)_pg\.txt$)");
            std::smatch m;
            if (std::regex_search(cfg.input_file, m, re)) {
                run = std::stoi(m[1].str());
            }
            auto it = (run >= 0) ? RUN_PHI.find(run) : RUN_PHI.end();
            if (it == RUN_PHI.end() && !cfg.phi_set) {
                throw std::runtime_error(
                    "No phi defined for this run. Set RUN_PHI in code or pass --phi-nonzero."
                );
            }
            const double phi = (it != RUN_PHI.end()) ? it->second : cfg.phi_nonzero;
            process_file(cfg.input_file, ofs, phi, run);
            std::cout << "Output: " << cfg.output_file << "\n";
        }

        return 0;

    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }
}
