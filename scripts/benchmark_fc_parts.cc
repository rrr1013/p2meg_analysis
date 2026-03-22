// scripts/benchmark_fc_parts.cc
//
// 使い方:
//   ./build/benchmark_fc_parts
//       <datafile> <N_mu_eff_nom> <N_mu_eff_sigma> <BR_test> <n_toys> [seed]
//
// 目的:
//   FC の 1 点評価で、toy 生成と q_toy 計算（free fit + fixed fit）の
//   どちらが時間を支配しているかを切り分ける。

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "p2meg/ACCGridPdf.h"
#include "p2meg/AnalysisWindow.h"
#include "p2meg/Constants.h"
#include "p2meg/DetectorResolution.h"
#include "p2meg/Event.h"
#include "p2meg/MathUtils.h"
#include "p2meg/NLLFit.h"
#include "p2meg/PdfWrappers.h"
#include "p2meg/RMDGridPdf.h"
#include "p2meg/UpperLimit.h"

static const char* kDefaultRmdRoot = "data/pdf_cache/rmd_grid.root";
static const char* kDefaultRmdKey  = "rmd_grid";
static const char* kDefaultAccRoot = "data/pdf_cache/acc_grid.root";
static const char* kDefaultAccKey  = "acc_grid";

static bool ParseDoubleStrict(const char* s, double& out)
{
    if (!s) return false;
    char* end = nullptr;
    out = std::strtod(s, &end);
    if (end == s || *end != '\0') return false;
    return std::isfinite(out);
}

static bool ParseIntStrict(const char* s, int& out)
{
    if (!s) return false;
    char* end = nullptr;
    const long v = std::strtol(s, &end, 10);
    if (end == s || *end != '\0') return false;
    out = static_cast<int>(v);
    return true;
}

static bool ParseULLStrict(const char* s, unsigned long long& out)
{
    if (!s) return false;
    char* end = nullptr;
    out = std::strtoull(s, &end, 10);
    if (end == s || *end != '\0') return false;
    return true;
}

static bool ParseDoublesFromLine(const std::string& line, std::vector<double>& out)
{
    out.clear();
    std::istringstream iss(line);

    std::string first;
    if (!(iss >> first)) return false;
    if (!first.empty() && first[0] == '#') return false;

    char* endptr = nullptr;
    const double v0 = std::strtod(first.c_str(), &endptr);
    if (endptr == first.c_str() || *endptr != '\0') return false;

    out.push_back(v0);
    double v = 0.0;
    while (iss >> v) out.push_back(v);
    return true;
}

static bool ThetaFromPhiPair(const Event& ev, double& theta_out)
{
    int idx_e = -1;
    int idx_g = -1;
    if (!Detector_IsAllowedPhiPairValue(ev.phi_detector_e, ev.phi_detector_g,
                                        detres, idx_e, idx_g)) {
        return false;
    }

    const int N_phi_e = Math_GetNPhiE(detres);
    const int N_phi_g = Math_GetNPhiG(detres);
    const double phi_e =
        Detector_PhiGridPoint(idx_e, detres.phi_e_min, detres.phi_e_max, N_phi_e);
    const double phi_g =
        Detector_PhiGridPoint(idx_g, detres.phi_g_min, detres.phi_g_max, N_phi_g);
    theta_out = std::fabs(phi_e - phi_g);
    return std::isfinite(theta_out);
}

static bool IsInsideAnalysisWindow(const Event& ev)
{
    if (!std::isfinite(ev.Ee) || !std::isfinite(ev.Eg) ||
        !std::isfinite(ev.t) || !std::isfinite(ev.phi_detector_e) ||
        !std::isfinite(ev.phi_detector_g)) {
        return false;
    }

    double theta = 0.0;
    if (!ThetaFromPhiPair(ev, theta)) return false;

    if (ev.Ee < analysis_window.Ee_min || ev.Ee > analysis_window.Ee_max) return false;
    if (ev.Eg < analysis_window.Eg_min || ev.Eg > analysis_window.Eg_max) return false;
    if (ev.t  < analysis_window.t_min  || ev.t  > analysis_window.t_max ) return false;
    if (theta < analysis_window.theta_min || theta > analysis_window.theta_max) return false;
    return true;
}

static bool LoadEventsFromDat(const char* filepath,
                              std::vector<Event>& events)
{
    events.clear();

    std::ifstream fin(filepath);
    if (!fin) {
        std::cerr << "[benchmark_fc_parts] cannot open: " << filepath << "\n";
        return false;
    }

    std::string line;
    std::vector<double> cols;
    while (std::getline(fin, line)) {
        if (!ParseDoublesFromLine(line, cols)) continue;
        if (cols.size() != 5) continue;

        Event ev{};
        ev.Ee = cols[0];
        ev.Eg = cols[1];
        ev.t = cols[2];
        ev.phi_detector_e = cols[3];
        ev.phi_detector_g = cols[4];

        if (!IsInsideAnalysisWindow(ev)) continue;
        events.push_back(ev);
    }

    return !events.empty();
}

int main(int argc, char** argv)
{
    if (argc != 6 && argc != 7) {
        std::cerr
            << "Usage: " << argv[0]
            << " <datafile> <N_mu_eff_nom> <N_mu_eff_sigma> <BR_test> <n_toys> [seed]\n";
        return 1;
    }

    const char* datafile = argv[1];
    double N_mu_eff_nom = 0.0;
    double N_mu_eff_sigma = 0.0;
    double BR_test = 0.0;
    int n_toys = 0;
    unsigned long long seed = 123456789ULL;

    if (!ParseDoubleStrict(argv[2], N_mu_eff_nom) ||
        !ParseDoubleStrict(argv[3], N_mu_eff_sigma) ||
        !ParseDoubleStrict(argv[4], BR_test) ||
        !ParseIntStrict(argv[5], n_toys)) {
        std::cerr << "[benchmark_fc_parts] failed to parse arguments\n";
        return 1;
    }
    if (argc == 7 && !ParseULLStrict(argv[6], seed)) {
        std::cerr << "[benchmark_fc_parts] failed to parse seed\n";
        return 1;
    }

    std::vector<Event> events;
    if (!LoadEventsFromDat(datafile, events)) {
        std::cerr << "[benchmark_fc_parts] no valid events in analysis window\n";
        return 2;
    }

    if (!RMDGridPdf_Load(kDefaultRmdRoot, kDefaultRmdKey)) {
        std::cerr << "[benchmark_fc_parts] RMDGridPdf_Load failed\n";
        return 3;
    }
    if (!ACCGridPdf_Load(kDefaultAccRoot, kDefaultAccKey)) {
        std::cerr << "[benchmark_fc_parts] ACCGridPdf_Load failed\n";
        return 3;
    }

    static SignalPdfContext sigctx{
        analysis_window,
        detres,
        kMassesPDG
    };

    std::vector<PdfComponent> components;
    components.push_back(MakeSignalComponent(&sigctx));
    components.push_back(MakeRMDComponent());
    components.push_back(MakeACCComponent());

    const double N0 = static_cast<double>(events.size());
    FitConfig free_cfg;
    free_cfg.start_yields = {N0 / 3.0, N0 / 3.0, N0 / 3.0};
    free_cfg.max_calls = 20000;
    free_cfg.tol = 1e-3;

    FitConfig prof_cfg = free_cfg;

    NormalizationUncertaintyConfig norm_cfg;
    norm_cfg.N_mu_eff_nom = N_mu_eff_nom;
    norm_cfg.N_mu_eff_sigma = N_mu_eff_sigma;

    UpperLimitBRPointConfig point_cfg;
    point_cfg.BR_test = BR_test;
    point_cfg.n_toys = 1;
    point_cfg.cl = 0.90;
    point_cfg.free_fit_cfg = free_cfg;
    point_cfg.prof_fit_cfg = prof_cfg;
    point_cfg.toy_cfg.seed = seed;
    point_cfg.toy_cfg.pmax_scan_trials = 20000;
    point_cfg.toy_cfg.pmax_safety = 5.0;
    point_cfg.toy_cfg.pmax_update = 1.2;
    point_cfg.toy_cfg.event_pool_size_per_component = 20000;
    point_cfg.norm_cfg = norm_cfg;

    UpperLimitBRPointResult obs = EvaluateUpperLimitBRPoint(events, components, point_cfg);
    if (!std::isfinite(obs.q_obs) || obs.fit_prof_obs.yields_hat.size() != components.size()) {
        std::cerr << "[benchmark_fc_parts] failed to build observed point\n";
        return 4;
    }

    std::vector<double> prof_mean_yields = obs.fit_prof_obs.yields_hat;
    std::vector<double> pmax_cache(components.size(), 0.0);

    double total_gen_ms = 0.0;
    double total_q_ms = 0.0;
    long long total_events = 0;
    int n_valid = 0;

    for (int itoy = 0; itoy < n_toys; ++itoy) {
        std::mt19937_64 rng_norm(seed + 0xB0000000ULL + static_cast<unsigned long long>(itoy));
        double N_mu_eff_toy = N_mu_eff_nom;
        if (!std::isfinite(N_mu_eff_sigma) || N_mu_eff_sigma < 0.0) continue;
        if (N_mu_eff_sigma > 0.0) {
            std::normal_distribution<double> gaus(N_mu_eff_nom, N_mu_eff_sigma);
            for (int itry = 0; itry < 10000; ++itry) {
                const double v = gaus(rng_norm);
                if (std::isfinite(v) && v > 0.0) {
                    N_mu_eff_toy = v;
                    break;
                }
            }
        }

        std::vector<double> mean_yields = prof_mean_yields;
        mean_yields[0] = BR_test * N_mu_eff_toy;
        if (!std::isfinite(mean_yields[0]) || mean_yields[0] < 0.0) continue;

        std::vector<Event> toy_events;
        std::vector<double> toy_generated_yields;

        const auto t0 = std::chrono::steady_clock::now();
        if (!GenerateToyDatasetFromModel(components, mean_yields, point_cfg.toy_cfg,
                                         static_cast<unsigned long long>(itoy),
                                         toy_events, &pmax_cache, &toy_generated_yields)) {
            continue;
        }
        const auto t1 = std::chrono::steady_clock::now();
        if (toy_events.empty()) continue;

        FitConfig toy_free_cfg = free_cfg;
        if (toy_free_cfg.start_yields.size() == components.size()) {
            toy_free_cfg.start_yields = toy_generated_yields;
            const double nsum =
                std::accumulate(toy_free_cfg.start_yields.begin(),
                                toy_free_cfg.start_yields.end(), 0.0);
            if (!(nsum > 0.0)) {
                for (double& v : toy_free_cfg.start_yields) v = 1.0;
            }
        }

        FitConfig toy_prof_cfg = prof_cfg;
        if (toy_prof_cfg.start_yields.size() == components.size()) {
            toy_prof_cfg.start_yields = toy_generated_yields;
            toy_prof_cfg.start_yields[0] = obs.N_sig_test_nominal;
        }

        FitResult fit_free_toy;
        FitResult fit_prof_toy;
        const auto t2 = std::chrono::steady_clock::now();
        const double q_toy = EvaluateProfileLikelihoodQ(toy_events, components,
                                                        toy_free_cfg, toy_prof_cfg,
                                                        obs.N_sig_test_nominal,
                                                        fit_free_toy, fit_prof_toy);
        const auto t3 = std::chrono::steady_clock::now();
        if (!std::isfinite(q_toy)) continue;

        total_gen_ms +=
            std::chrono::duration<double, std::milli>(t1 - t0).count();
        total_q_ms +=
            std::chrono::duration<double, std::milli>(t3 - t2).count();
        total_events += static_cast<long long>(toy_events.size());
        ++n_valid;
    }

    std::cout << std::setprecision(10);
    std::cout << "================ benchmark_fc_parts ================\n";
    std::cout << "datafile             = " << datafile << "\n";
    std::cout << "events_used          = " << events.size() << "\n";
    std::cout << "BR_test              = " << BR_test << "\n";
    std::cout << "N_sig_test_nominal   = " << obs.N_sig_test_nominal << "\n";
    std::cout << "q_obs                = " << obs.q_obs << "\n";
    std::cout << "n_toys_requested     = " << n_toys << "\n";
    std::cout << "n_toys_valid         = " << n_valid << "\n";
    std::cout << "toy_mean_events      = "
              << ((n_valid > 0) ? static_cast<double>(total_events) / n_valid : 0.0) << "\n";
    std::cout << "toy_gen_total_ms     = " << total_gen_ms << "\n";
    std::cout << "toy_q_total_ms       = " << total_q_ms << "\n";
    std::cout << "toy_gen_mean_ms      = "
              << ((n_valid > 0) ? total_gen_ms / n_valid : 0.0) << "\n";
    std::cout << "toy_q_mean_ms        = "
              << ((n_valid > 0) ? total_q_ms / n_valid : 0.0) << "\n";
    std::cout << "q_over_gen_ratio     = "
              << ((total_gen_ms > 0.0) ? total_q_ms / total_gen_ms : 0.0) << "\n";
    std::cout << "====================================================\n";

    return 0;
}
