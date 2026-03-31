#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace fs = std::filesystem;

struct Config {
    std::string input_dir = "data/shapeddata";
    std::string output_dir = "data/shapeddata";
    std::string input_file;
    std::string out_prefix;
    std::string out_eg;
    int run = -1;
    bool eg_all_patterns = false;
    double ab_window_ns = 100.0;
};

struct Rec {
    int run = -1;
    long long event = -1;
    std::string module;   // A or B
    double energy_like = 0.0;   // 互換列名 integral の実体（エネルギー相当量）
    // final_module_pulses_runN の t_fall_start_ns:
    //  - 元の t_fall_start_sample_raw（テンプレ基準）へ CHオフセットを加え、
    //    4 ns/sample で ns化した時刻
    double t_fall_start_ns = 0.0;
    std::string particle; // positron or gamma
    std::string nai_ch;
    std::string ps_ch;
};

// サンプリング周期 [ns/sample]
static const double kSamplePeriodNs = 4.0;
// ns -> 実時間 [s]
static const double kNsToSec = 1.0e-9;
static const double kPi = 3.14159265358979323846;

// 方位角設定:
//  - module B は常に +pi/2
//  - module A は run で切り替え
//  - この run 区分は測定配置（A側検出器の設置角）に対応
//      run 1..7   : -pi/2
//      run 8..14  : -pi/6
//      run 15..21 : +pi/6
static const double kPhiModuleB = 0.5 * kPi;

static bool PhiModuleAFromRun(int run, double& phi_a) {
    if (run >= 1 && run <= 7) {
        phi_a = -0.5 * kPi;
        return true;
    }
    if (run >= 8 && run <= 14) {
        phi_a = -kPi / 6.0;
        return true;
    }
    if (run >= 15 && run <= 21) {
        phi_a = kPi / 6.0;
        return true;
    }
    return false;
}

static void Usage(const char* prog) {
    std::cerr
        << "Usage:\n"
        << "  " << prog << " [options]\n\n"
        << "Options:\n"
        << "  --input FILE       input file (default: auto-detect final_module_pulses_runM.txt)\n"
        << "  --input-dir DIR    input directory (default: data/shapeddata)\n"
        << "  --output-dir DIR   output directory (default: data/shapeddata)\n"
        << "  --out-prefix PFX   step3出力prefix (default: step3f_runM in output-dir)\n"
        << "  --out-eg FILE      step4出力 (default: step4f_runM_eg.txt in output-dir)\n"
        << "  --run N            run number\n"
        << "  --ab-window-ns X   A-Bペア許容時間差[ns] (default: 100)\n"
        << "  --eg-all-patterns  step4(eg)はA-Bの全組合せ(|dt|窓内)で出力\n";
}

static std::pair<int, std::string> DiscoverInput(const std::string& input_dir, int run_filter) {
    std::regex re(R"(final_module_pulses_run([0-9]+)\.txt$)");
    int run = -1;
    bool set = false;
    std::string path;

    for (const auto& ent : fs::recursive_directory_iterator(input_dir)) {
        if (!ent.is_regular_file()) continue;
        const std::string fname = ent.path().filename().string();
        std::smatch m;
        if (!std::regex_match(fname, m, re)) continue;
        const int r = std::stoi(m[1].str());
        if (run_filter >= 0 && r != run_filter) continue;
        if (!set) {
            run = r;
            set = true;
            path = ent.path().string();
        } else if (r != run) {
            throw std::runtime_error("Multiple run files found. Use --run or --input.");
        }
    }

    if (!set) {
        if (run_filter >= 0) {
            throw std::runtime_error("No final_module_pulses_run" + std::to_string(run_filter) + ".txt found.");
        }
        throw std::runtime_error("No final_module_pulses_runM.txt found.");
    }
    return {run, path};
}

static int InferRunFromPath(const std::string& path) {
    std::regex re(R"(final_module_pulses_run([0-9]+)\.txt$)");
    std::smatch m;
    if (std::regex_search(path, m, re)) return std::stoi(m[1].str());
    return -1;
}

static Config ParseArgs(int argc, char** argv) {
    Config cfg;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto need = [&](const std::string& opt) {
            if (i + 1 >= argc) throw std::runtime_error("Missing value for " + opt);
            return std::string(argv[++i]);
        };

        if (a == "--input") cfg.input_file = need(a);
        else if (a == "--input-dir") cfg.input_dir = need(a);
        else if (a == "--output-dir") cfg.output_dir = need(a);
        else if (a == "--out-prefix") cfg.out_prefix = need(a);
        else if (a == "--out-eg") cfg.out_eg = need(a);
        else if (a == "--run") cfg.run = std::stoi(need(a));
        else if (a == "--ab-window-ns") cfg.ab_window_ns = std::stod(need(a));
        else if (a == "--eg-all-patterns") cfg.eg_all_patterns = true;
        else if (a == "-h" || a == "--help") { Usage(argv[0]); std::exit(0); }
        else if (!a.empty() && a[0] != '-') {
            if (cfg.input_file.empty()) cfg.input_file = a;
            else throw std::runtime_error("Too many positional arguments: " + a);
        } else {
            throw std::runtime_error("Unknown option: " + a);
        }
    }
    return cfg;
}

static bool ParseLine(const std::string& line, Rec& r) {
    if (line.empty()) return false;
    if (line.find("run event module") != std::string::npos) return false;

    std::istringstream iss(line);
    if (!(iss >> r.run >> r.event >> r.module >> r.energy_like >> r.t_fall_start_ns >> r.particle >> r.nai_ch >> r.ps_ch)) {
        return false;
    }
    if (r.module != "A" && r.module != "B") return false;
    if (r.particle != "positron" && r.particle != "gamma") return false;
    return true;
}

struct MatchScore {
    int pairs = -1;
    double sum_abs_dt = std::numeric_limits<double>::infinity();
};

static bool BetterScore(const MatchScore& a, const MatchScore& b) {
    if (a.pairs != b.pairs) return a.pairs > b.pairs;
    return a.sum_abs_dt + 1.0e-12 < b.sum_abs_dt;
}

// 時刻1次元に対する1対1全体最適化:
//  1) ペア数最大化
//  2) 同数なら Σ|Δt| 最小化
// A/B を昇順ソート済みとして DP で厳密に解く。
static std::vector<std::pair<int,int>> OptimizeABPairsGlobal(const std::vector<Rec>& va,
                                                             const std::vector<Rec>& vb,
                                                             double window_ns) {
    const int na = static_cast<int>(va.size());
    const int nb = static_cast<int>(vb.size());
    const int ncol = nb + 1;
    const int nstate = (na + 1) * (nb + 1);
    auto idx = [ncol](int i, int j) { return i * ncol + j; };

    std::vector<MatchScore> dp(static_cast<size_t>(nstate));
    std::vector<int> prev_i(static_cast<size_t>(nstate), -1);
    std::vector<int> prev_j(static_cast<size_t>(nstate), -1);
    std::vector<char> prev_act(static_cast<size_t>(nstate), 0); // 'A':skip A, 'B':skip B, 'M':match

    dp[0] = MatchScore{0, 0.0};
    prev_act[0] = 'S';

    for (int i = 0; i <= na; ++i) {
        for (int j = 0; j <= nb; ++j) {
            const int cur = idx(i, j);
            if (dp[cur].pairs < 0) continue;

            // A側を未使用にする
            if (i < na) {
                const int nxt = idx(i + 1, j);
                const MatchScore cand = dp[cur];
                if (BetterScore(cand, dp[nxt])) {
                    dp[nxt] = cand;
                    prev_i[nxt] = i;
                    prev_j[nxt] = j;
                    prev_act[nxt] = 'A';
                }
            }
            // B側を未使用にする
            if (j < nb) {
                const int nxt = idx(i, j + 1);
                const MatchScore cand = dp[cur];
                if (BetterScore(cand, dp[nxt])) {
                    dp[nxt] = cand;
                    prev_i[nxt] = i;
                    prev_j[nxt] = j;
                    prev_act[nxt] = 'B';
                }
            }
            // A_i と B_j をマッチ
            if (i < na && j < nb) {
                const double dt = std::fabs(va[static_cast<size_t>(i)].t_fall_start_ns -
                                            vb[static_cast<size_t>(j)].t_fall_start_ns);
                if (dt <= window_ns) {
                    const int nxt = idx(i + 1, j + 1);
                    const MatchScore cand{dp[cur].pairs + 1, dp[cur].sum_abs_dt + dt};
                    if (BetterScore(cand, dp[nxt])) {
                        dp[nxt] = cand;
                        prev_i[nxt] = i;
                        prev_j[nxt] = j;
                        prev_act[nxt] = 'M';
                    }
                }
            }
        }
    }

    std::vector<std::pair<int,int>> pairs;
    int i = na;
    int j = nb;
    while (!(i == 0 && j == 0)) {
        const int cur = idx(i, j);
        const char act = prev_act[cur];
        const int pi = prev_i[cur];
        const int pj = prev_j[cur];
        if (act == 'M') {
            pairs.emplace_back(i - 1, j - 1);
        }
        if (pi < 0 || pj < 0) break;
        i = pi;
        j = pj;
    }
    std::reverse(pairs.begin(), pairs.end());
    return pairs;
}

static void WriteStep3Header(std::ofstream& ofs) {
    // 互換性のため列名 integral を維持（値は energy_like）
    ofs << "run\tevent\tpair_id\tmodule\tintegral\tt_fall_start_ns\tt_fall_start_sec\tparticle\tnai_ch\tps_ch\n";
}

static void WriteStep3Rec(std::ofstream& ofs, int run, long long event, int pair_id, const Rec& r) {
    ofs << run << "\t" << event << "\t" << pair_id << "\t"
        << r.module << "\t" << std::abs(r.energy_like) << "\t" << r.t_fall_start_ns << "\t"
        << (r.t_fall_start_ns * kNsToSec) << "\t"
        << r.particle << "\t" << r.nai_ch << "\t" << r.ps_ch << "\n";
}

static bool WriteEgIfOppositePair(std::ofstream& ofs_eg,
                                  const Rec& a,
                                  const Rec& b,
                                  double phi_module_a) {
    const Rec* pos = nullptr;
    const Rec* gam = nullptr;
    if (a.particle == "positron" && b.particle == "gamma") {
        pos = &a;
        gam = &b;
    } else if (a.particle == "gamma" && b.particle == "positron") {
        pos = &b;
        gam = &a;
    }
    if (!pos || !gam) return false;

    const double t = (pos->t_fall_start_ns - gam->t_fall_start_ns) * kNsToSec;
    const double phi_e = (pos->module == "A") ? phi_module_a : kPhiModuleB;
    const double phi_g = (gam->module == "A") ? phi_module_a : kPhiModuleB;

    ofs_eg << std::fixed << std::setprecision(6) << std::abs(pos->energy_like) << "\t"
           << std::fixed << std::setprecision(6) << std::abs(gam->energy_like) << "\t"
           << std::scientific << std::setprecision(6) << t << "\t"
           << std::fixed << std::setprecision(6) << phi_e << "\t"
           << std::fixed << std::setprecision(6) << phi_g << "\n";
    return true;
}

int main(int argc, char** argv) {
    try {
        Config cfg = ParseArgs(argc, argv);
        if (cfg.ab_window_ns <= 0.0) {
            throw std::runtime_error("ab-window-ns must be > 0.");
        }

        int run = -1;
        if (cfg.input_file.empty()) {
            auto [r, p] = DiscoverInput(cfg.input_dir, cfg.run);
            run = r;
            cfg.input_file = p;
        }
        if (run < 0) run = (cfg.run >= 0) ? cfg.run : InferRunFromPath(cfg.input_file);
        if (run < 0) throw std::runtime_error("Cannot infer run number.");

        if (cfg.out_prefix.empty()) {
            cfg.out_prefix = (fs::path(cfg.output_dir) / ("step3f_run" + std::to_string(run))).string();
        }
        if (cfg.out_eg.empty()) {
            cfg.out_eg = (fs::path(cfg.output_dir) / ("step4f_run" + std::to_string(run) + "_eg.txt")).string();
        }

        double phi_module_a = 0.0;
        if (!PhiModuleAFromRun(run, phi_module_a)) {
            throw std::runtime_error("No phi setting for run" + std::to_string(run) +
                                     ". Supported runs: 1..21");
        }

        std::ifstream ifs(cfg.input_file);
        if (!ifs) throw std::runtime_error("Cannot open input: " + cfg.input_file);

        std::ofstream ofs_pp(cfg.out_prefix + "_pp.txt");
        std::ofstream ofs_pg(cfg.out_prefix + "_pg.txt");
        std::ofstream ofs_gg(cfg.out_prefix + "_gg.txt");
        std::ofstream ofs_eg(cfg.out_eg);
        if (!ofs_pp || !ofs_pg || !ofs_gg || !ofs_eg) throw std::runtime_error("Cannot open output files.");

        WriteStep3Header(ofs_pp);
        WriteStep3Header(ofs_pg);
        WriteStep3Header(ofs_gg);
        ofs_eg << "E_e\tE_g\tt\tphi_e\tphi_g\n";

        std::map<std::pair<int,long long>, std::vector<Rec>> events;
        std::string line;
        while (std::getline(ifs, line)) {
            Rec r;
            if (!ParseLine(line, r)) continue;
            events[{r.run, r.event}].push_back(r);
        }

        size_t n_pair = 0;
        size_t n_skip = 0;
        size_t n_eg = 0;

        for (auto& kv : events) {
            const int r = kv.first.first;
            const long long ev = kv.first.second;
            std::vector<Rec> va;
            std::vector<Rec> vb;
            for (const auto& rec : kv.second) {
                if (rec.module == "A") va.push_back(rec);
                else if (rec.module == "B") vb.push_back(rec);
            }
            std::sort(va.begin(), va.end(), [](const Rec& l, const Rec& r2){ return l.t_fall_start_ns < r2.t_fall_start_ns; });
            std::sort(vb.begin(), vb.end(), [](const Rec& l, const Rec& r2){ return l.t_fall_start_ns < r2.t_fall_start_ns; });

            // 1対1全体最適化:
            //  1) ペア数最大化
            //  2) 同数なら Σ|Δt| 最小化
            const auto pairs = OptimizeABPairsGlobal(va, vb, cfg.ab_window_ns);
            int pair_id = 0;
            for (const auto& ij : pairs) {
                const int ia = ij.first;
                const int ib = ij.second;
                const int pid = pair_id++;
                const Rec& a = va[static_cast<size_t>(ia)];
                const Rec& b = vb[static_cast<size_t>(ib)];

                if (a.particle == "positron" && b.particle == "positron") {
                    WriteStep3Rec(ofs_pp, r, ev, pid, a);
                    WriteStep3Rec(ofs_pp, r, ev, pid, b);
                } else if (a.particle == "gamma" && b.particle == "gamma") {
                    WriteStep3Rec(ofs_gg, r, ev, pid, a);
                    WriteStep3Rec(ofs_gg, r, ev, pid, b);
                } else {
                    WriteStep3Rec(ofs_pg, r, ev, pid, a);
                    WriteStep3Rec(ofs_pg, r, ev, pid, b);

                    if (!cfg.eg_all_patterns) {
                        if (WriteEgIfOppositePair(ofs_eg, a, b, phi_module_a)) {
                            ++n_eg;
                        }
                    }
                }
                ++n_pair;
            }

            // step4(eg)のみ全組合せで作るモード:
            // A/B の全候補から |dt| 窓内かつ異種粒子(positron,gamma)を全て採用
            if (cfg.eg_all_patterns) {
                for (const auto& a : va) {
                    for (const auto& b : vb) {
                        const double dt = std::fabs(a.t_fall_start_ns - b.t_fall_start_ns);
                        if (dt > cfg.ab_window_ns) continue;
                        if (WriteEgIfOppositePair(ofs_eg, a, b, phi_module_a)) {
                            ++n_eg;
                        }
                    }
                }
            }
            n_skip += (va.size() + vb.size() - 2 * pairs.size());
        }

        std::cout << "Input: " << cfg.input_file << "\n";
        std::cout << "Step3 output: " << cfg.out_prefix << "_{pp,pg,gg}.txt\n";
        std::cout << "Step4 output: " << cfg.out_eg << "\n";
        std::cout << "AB pair window: " << cfg.ab_window_ns << " ns (" << (cfg.ab_window_ns / kSamplePeriodNs)
                  << " samples)\n";
        std::cout << "EG matching mode: " << (cfg.eg_all_patterns ? "all-patterns" : "from-1to1-pairs") << "\n";
        std::cout << "Pairs written: " << n_pair << "\n";
        std::cout << "Unpaired pulses skipped: " << n_skip << "\n";
        std::cout << "EG events written: " << n_eg << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
}
