// scripts/diagnose_fc_qtoy_steps.cc
//
// 使い方:
//   ./build/diagnose_fc_qtoy_steps
//       <datafile> <N_mu_eff_nom> <N_mu_eff_sigma> <BR_test> <n_toys> [seed]
//
// 目的:
//   図8の p-value の段構造が何に由来するかを調べるため、
//   ある BR 点で FC toy の q_toy 分布を生成事象数ごとに要約する。

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
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

struct QSummary {
    int n_total;
    int n_ge;
    std::vector<double> q_values;
};

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
                              std::vector<Event>& events,
                              long& n_skipped_non5,
                              long& n_skipped_outside)
{
    events.clear();
    n_skipped_non5 = 0;
    n_skipped_outside = 0;

    std::ifstream fin(filepath);
    if (!fin) return false;

    std::string line;
    std::vector<double> cols;
    while (std::getline(fin, line)) {
        if (!ParseDoublesFromLine(line, cols)) continue;
        if (cols.size() != 5) {
            ++n_skipped_non5;
            continue;
        }

        Event ev{};
        ev.Ee = cols[0];
        ev.Eg = cols[1];
        ev.t = cols[2];
        ev.phi_detector_e = cols[3];
        ev.phi_detector_g = cols[4];

        if (!IsInsideAnalysisWindow(ev)) {
            ++n_skipped_outside;
            continue;
        }
        events.push_back(ev);
    }

    return !events.empty();
}

static double Quantile(std::vector<double> values, double prob)
{
    if (values.empty()) return 0.0;
    std::sort(values.begin(), values.end());
    if (values.size() == 1U) return values.front();
    const double x = prob * static_cast<double>(values.size() - 1U);
    const std::size_t i0 = static_cast<std::size_t>(std::floor(x));
    const std::size_t i1 = static_cast<std::size_t>(std::ceil(x));
    const double t = x - static_cast<double>(i0);
    return (1.0 - t) * values[i0] + t * values[i1];
}

static unsigned long long MixLocalSeed(unsigned long long seed, unsigned long long salt)
{
    unsigned long long z = seed + 0x9e3779b97f4a7c15ULL + salt;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    z = z ^ (z >> 31);
    return z;
}

static bool SamplePositiveGaussianLocal(double mean,
                                        double sigma,
                                        std::mt19937_64& rng,
                                        double& sampled_out)
{
    sampled_out = 0.0;
    if (!std::isfinite(mean) || !(mean > 0.0)) return false;
    if (!std::isfinite(sigma) || sigma <= 0.0) {
        sampled_out = mean;
        return true;
    }

    std::normal_distribution<double> gaus(mean, sigma);
    for (int itry = 0; itry < 10000; ++itry) {
        const double v = gaus(rng);
        if (std::isfinite(v) && v > 0.0) {
            sampled_out = v;
            return true;
        }
    }

    sampled_out = mean;
    return true;
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
        std::cerr << "[diagnose_fc_qtoy_steps] failed to parse arguments\n";
        return 1;
    }
    if (argc == 7 && !ParseULLStrict(argv[6], seed)) {
        std::cerr << "[diagnose_fc_qtoy_steps] failed to parse seed\n";
        return 1;
    }

    std::vector<Event> data_events;
    long n_skipped_non5 = 0;
    long n_skipped_outside = 0;
    if (!LoadEventsFromDat(datafile, data_events, n_skipped_non5, n_skipped_outside)) {
        std::cerr << "[diagnose_fc_qtoy_steps] failed to load events from " << datafile << "\n";
        return 2;
    }

    if (!RMDGridPdf_Load(kDefaultRmdRoot, kDefaultRmdKey)) {
        std::cerr << "[diagnose_fc_qtoy_steps] RMDGridPdf_Load failed\n";
        return 3;
    }
    if (!ACCGridPdf_Load(kDefaultAccRoot, kDefaultAccKey)) {
        std::cerr << "[diagnose_fc_qtoy_steps] ACCGridPdf_Load failed\n";
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

    const double N0 = static_cast<double>(data_events.size());
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

    UpperLimitBRPointResult obs = EvaluateUpperLimitBRPoint(data_events, components, point_cfg);
    if (!std::isfinite(obs.q_obs) || obs.fit_prof_obs.yields_hat.size() != components.size()) {
        std::cerr << "[diagnose_fc_qtoy_steps] failed to build observed point\n";
        return 4;
    }

    std::vector<double> prof_mean_yields = obs.fit_prof_obs.yields_hat;
    std::vector<double> pmax_cache(components.size(), 0.0);

    std::vector<double> all_q;
    all_q.reserve(static_cast<std::size_t>(n_toys));
    std::map<int, QSummary> by_nsig;
    std::map<int, QSummary> by_nprompt;
    int n_valid = 0;
    int n_ge_total = 0;

    for (int itoy = 0; itoy < n_toys; ++itoy) {
        std::mt19937_64 rng_norm(
            MixLocalSeed(seed, static_cast<unsigned long long>(0xB0000000ULL + static_cast<unsigned long long>(itoy))));

        double N_mu_eff_toy = N_mu_eff_nom;
        if (!SamplePositiveGaussianLocal(N_mu_eff_nom, N_mu_eff_sigma, rng_norm, N_mu_eff_toy)) {
            continue;
        }

        std::vector<double> mean_yields = prof_mean_yields;
        mean_yields[0] = BR_test * N_mu_eff_toy;
        if (!std::isfinite(mean_yields[0]) || mean_yields[0] < 0.0) continue;

        std::vector<Event> toy_events;
        std::vector<double> toy_generated_yields;
        if (!GenerateToyDatasetFromModel(components, mean_yields, point_cfg.toy_cfg,
                                         static_cast<unsigned long long>(itoy),
                                         toy_events, &pmax_cache, &toy_generated_yields)) {
            continue;
        }
        if (toy_events.empty()) continue;
        if (toy_generated_yields.size() != components.size()) continue;

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
        const double q_toy = EvaluateProfileLikelihoodQ(toy_events, components,
                                                        toy_free_cfg, toy_prof_cfg,
                                                        obs.N_sig_test_nominal,
                                                        fit_free_toy, fit_prof_toy);
        if (!std::isfinite(q_toy)) continue;

        ++n_valid;
        all_q.push_back(q_toy);
        if (q_toy >= obs.q_obs) ++n_ge_total;

        const int n_sig_gen = static_cast<int>(std::llround(toy_generated_yields[0]));
        const int n_rmd_gen = static_cast<int>(std::llround(toy_generated_yields[1]));
        QSummary& s_sig = by_nsig[n_sig_gen];
        ++s_sig.n_total;
        s_sig.q_values.push_back(q_toy);
        if (q_toy >= obs.q_obs) ++s_sig.n_ge;

        QSummary& s_prompt = by_nprompt[n_sig_gen + n_rmd_gen];
        ++s_prompt.n_total;
        s_prompt.q_values.push_back(q_toy);
        if (q_toy >= obs.q_obs) ++s_prompt.n_ge;
    }

    std::cout << std::setprecision(10);
    std::cout << "================ diagnose_fc_qtoy_steps ================\n";
    std::cout << "datafile              = " << datafile << "\n";
    std::cout << "events_used           = " << data_events.size() << "\n";
    std::cout << "BR_test               = " << BR_test << "\n";
    std::cout << "N_sig_test_nominal    = " << obs.N_sig_test_nominal << "\n";
    std::cout << "q_obs                 = " << obs.q_obs << "\n";
    std::cout << "fit_prof_mean_sig     = " << obs.fit_prof_obs.yields_hat[0] << "\n";
    std::cout << "fit_prof_mean_rmd     = " << obs.fit_prof_obs.yields_hat[1] << "\n";
    std::cout << "fit_prof_mean_acc     = " << obs.fit_prof_obs.yields_hat[2] << "\n";
    std::cout << "n_toys_requested      = " << n_toys << "\n";
    std::cout << "n_toys_valid          = " << n_valid << "\n";
    std::cout << "global_p_value        = "
              << ((n_valid > 0) ? static_cast<double>(n_ge_total) / static_cast<double>(n_valid) : 0.0) << "\n";
    std::cout << "q_toy_q50             = " << Quantile(all_q, 0.50) << "\n";
    std::cout << "q_toy_q90             = " << Quantile(all_q, 0.90) << "\n";
    std::cout << "q_toy_q95             = " << Quantile(all_q, 0.95) << "\n";
    std::cout << "========================================================\n";

    std::cout << "================ by generated N_sig ====================\n";
    std::cout << "# N_sig_gen  n_toys  frac  frac_ge_qobs  q50  q90  q95\n";
    for (const auto& kv : by_nsig) {
        const int n_sig_gen = kv.first;
        const QSummary& s = kv.second;
        if (s.n_total <= 0) continue;
        std::cout << n_sig_gen << " "
                  << s.n_total << " "
                  << static_cast<double>(s.n_total) / static_cast<double>(n_valid) << " "
                  << static_cast<double>(s.n_ge) / static_cast<double>(s.n_total) << " "
                  << Quantile(s.q_values, 0.50) << " "
                  << Quantile(s.q_values, 0.90) << " "
                  << Quantile(s.q_values, 0.95) << "\n";
    }
    std::cout << "========================================================\n";

    std::cout << "============= by generated (N_sig + N_rmd) =============\n";
    std::cout << "# N_prompt_gen  n_toys  frac  frac_ge_qobs  q50  q90  q95\n";
    for (const auto& kv : by_nprompt) {
        const int n_prompt_gen = kv.first;
        const QSummary& s = kv.second;
        if (s.n_total <= 0) continue;
        std::cout << n_prompt_gen << " "
                  << s.n_total << " "
                  << static_cast<double>(s.n_total) / static_cast<double>(n_valid) << " "
                  << static_cast<double>(s.n_ge) / static_cast<double>(s.n_total) << " "
                  << Quantile(s.q_values, 0.50) << " "
                  << Quantile(s.q_values, 0.90) << " "
                  << Quantile(s.q_values, 0.95) << "\n";
    }
    std::cout << "========================================================\n";

    return 0;
}
