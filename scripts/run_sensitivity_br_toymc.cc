// scripts/run_sensitivity_br_toymc.cc
//
// 使い方:
//   ./build/run_sensitivity_br_toymc
//       <datafile> <N_mu_eff_nom> <N_mu_eff_sigma>
//       <br_min> <br_max> <br_step>
//       <n_calib_toys> <n_sens_toys>
//       [seed] [n_validation_toys] [n_validation_inner]
//
// 目的:
//   - 既存の BR upper-limit toy MC を流用し、
//     background-only pseudo-experiment で得られる BR_90 分布の中央値
//     を感度として評価する。
//   - 計算時間短縮のため、FC belt の toy q 分布は
//     「N_sig=0 profile fit から得た nominal background point」で 1 回だけ較正し、
//     外側の sensitivity toy ではその belt を再利用する。
//   - この近似の健全性は、少数 toy に対して
//     既存 EvaluateUpperLimitBRScan（各 dataset ごとに内側 toy を回す）
//     と比較して確認する。

#include <algorithm>
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

struct CalibrationPoint {
    double BR_test;               // 分岐比仮説
    double N_sig_test_nominal;    // 公称 N_mu_eff に対応する signal yield
    std::vector<double> q_sorted; // q_toy を昇順に並べたもの
    int n_toys_valid;             // 有効 toy 数
};

struct LimitSummary {
    double BR_90;
    double N_sig_90_nominal;
    FitResult fit_free;
    std::vector<double> q_obs_scan;
    std::vector<double> p_value_scan;
    std::vector<int> accepted_scan;
    bool ok;
};

static constexpr int kLocalExactHalfWindow = 2;

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
    if (!fin) {
        std::cerr << "[run_sensitivity_br_toymc] cannot open: " << filepath << "\n";
        return false;
    }

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

static std::vector<double> BuildScanPoints(double x_min, double x_max, double x_step)
{
    std::vector<double> out;
    if (!std::isfinite(x_min) || !std::isfinite(x_max) || !std::isfinite(x_step)) return out;
    if (!(x_step > 0.0) || x_max < x_min) return out;

    for (double x = x_min; x <= x_max + 0.5 * x_step; x += x_step) {
        out.push_back(x);
    }
    return out;
}

static bool SamplePositiveGaussian(double mean,
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

static unsigned long long MixSeed(unsigned long long seed, unsigned long long salt)
{
    unsigned long long z = seed + 0x9e3779b97f4a7c15ULL + salt;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    z = z ^ (z >> 31);
    return z;
}

static double MedianFromSorted(const std::vector<double>& sorted)
{
    if (sorted.empty()) return 0.0;
    const std::size_t n = sorted.size();
    if ((n % 2U) == 1U) return sorted[n / 2U];
    return 0.5 * (sorted[n / 2U - 1U] + sorted[n / 2U]);
}

static double QuantileFromSorted(const std::vector<double>& sorted, double prob)
{
    if (sorted.empty()) return 0.0;
    if (!(prob >= 0.0) || !(prob <= 1.0) || !std::isfinite(prob)) return 0.0;
    if (sorted.size() == 1U) return sorted.front();

    const double x = prob * static_cast<double>(sorted.size() - 1U);
    const std::size_t i0 = static_cast<std::size_t>(std::floor(x));
    const std::size_t i1 = static_cast<std::size_t>(std::ceil(x));
    const double t = x - static_cast<double>(i0);
    return (1.0 - t) * sorted[i0] + t * sorted[i1];
}

static double PValueFromSortedQ(const std::vector<double>& q_sorted, double q_obs)
{
    if (q_sorted.empty() || !std::isfinite(q_obs)) return 0.0;
    const auto it = std::lower_bound(q_sorted.begin(), q_sorted.end(), q_obs);
    const std::size_t n_ge =
        static_cast<std::size_t>(q_sorted.end() - it);
    return static_cast<double>(n_ge) / static_cast<double>(q_sorted.size());
}

static bool BuildCalibrationDistributions(
    const std::vector<PdfComponent>& components,
    const std::vector<double>& br_scan,
    const std::vector<double>& n_sig_scan,
    const std::vector<std::vector<double>>& calib_mean_yields_by_point,
    const FitConfig& free_cfg_template,
    const FitConfig& prof_cfg_template,
    const ToyGeneratorConfig& toy_cfg,
    const NormalizationUncertaintyConfig& norm_cfg,
    int n_calib_toys,
    std::vector<double>& io_pmax_cache,
    std::vector<CalibrationPoint>& out_points)
{
    out_points.clear();
    out_points.reserve(br_scan.size());

    for (std::size_t ip = 0; ip < br_scan.size(); ++ip) {
        CalibrationPoint point;
        point.BR_test = br_scan[ip];
        point.N_sig_test_nominal = n_sig_scan[ip];
        point.q_sorted.clear();
        point.n_toys_valid = 0;

        for (int itoy = 0; itoy < n_calib_toys; ++itoy) {
            std::mt19937_64 rng_norm(
                MixSeed(toy_cfg.seed,
                        static_cast<unsigned long long>(0xC0000000ULL +
                        1000000ULL * static_cast<unsigned long long>(ip) +
                        static_cast<unsigned long long>(itoy))));

            double N_mu_eff_toy = norm_cfg.N_mu_eff_nom;
            if (!SamplePositiveGaussian(norm_cfg.N_mu_eff_nom,
                                        norm_cfg.N_mu_eff_sigma,
                                        rng_norm,
                                        N_mu_eff_toy)) {
                continue;
            }

            if (ip >= calib_mean_yields_by_point.size()) continue;
            std::vector<double> mean_yields = calib_mean_yields_by_point[ip];
            if (mean_yields.size() != components.size()) continue;
            mean_yields[0] = point.BR_test * N_mu_eff_toy;
            if (!std::isfinite(mean_yields[0]) || mean_yields[0] < 0.0) continue;

            std::vector<Event> toy_events;
            std::vector<double> generated_yields;
            if (!GenerateToyDatasetFromModel(components, mean_yields, toy_cfg,
                                             static_cast<unsigned long long>(
                                                 1000000ULL * static_cast<unsigned long long>(ip) +
                                                 static_cast<unsigned long long>(itoy)),
                                             toy_events, &io_pmax_cache,
                                             &generated_yields)) {
                continue;
            }
            if (toy_events.empty()) continue;

            FitConfig free_cfg = free_cfg_template;
            if (free_cfg.start_yields.size() == components.size()) {
                free_cfg.start_yields = generated_yields;
                const double nsum = std::accumulate(free_cfg.start_yields.begin(),
                                                    free_cfg.start_yields.end(), 0.0);
                if (!(nsum > 0.0)) {
                    for (double& v : free_cfg.start_yields) v = 1.0;
                }
            }

            FitConfig prof_cfg = prof_cfg_template;
            if (prof_cfg.start_yields.size() == components.size()) {
                prof_cfg.start_yields = generated_yields;
                prof_cfg.start_yields[0] = point.N_sig_test_nominal;
            }

            FitResult fit_free_toy;
            FitResult fit_prof_toy;
            const double q_toy = EvaluateProfileLikelihoodQ(
                toy_events, components, free_cfg, prof_cfg,
                point.N_sig_test_nominal, fit_free_toy, fit_prof_toy);
            if (!std::isfinite(q_toy)) continue;

            point.q_sorted.push_back(q_toy);
            ++point.n_toys_valid;
        }

        std::sort(point.q_sorted.begin(), point.q_sorted.end());
        out_points.push_back(point);
    }

    return true;
}

static LimitSummary EvaluateLimitFromCalibrationQScan(
    const ProfileLikelihoodQScanResult& qscan,
    const std::vector<CalibrationPoint>& calib_points,
    double cl,
    bool enforce_monotonic)
{
    LimitSummary out;
    out.BR_90 = 0.0;
    out.N_sig_90_nominal = 0.0;
    out.fit_free.status = 1;
    out.fit_free.nll_min = 0.0;
    out.q_obs_scan.clear();
    out.p_value_scan.clear();
    out.accepted_scan.clear();
    out.ok = false;

    out.fit_free = qscan.fit_free;
    if (qscan.fit_free.status != 0 || qscan.points.size() != calib_points.size()) {
        return out;
    }

    const double acceptance_threshold =
        (std::isfinite(cl) && cl > 0.0 && cl < 1.0) ? (1.0 - cl) : 0.1;

    out.q_obs_scan.reserve(calib_points.size());
    out.p_value_scan.reserve(calib_points.size());
    out.accepted_scan.reserve(calib_points.size());

    bool seen_reject = false;
    for (std::size_t i = 0; i < calib_points.size(); ++i) {
        const double q_obs = qscan.points[i].q_value;
        const double p_value = PValueFromSortedQ(calib_points[i].q_sorted, q_obs);
        int accepted = (p_value >= acceptance_threshold) ? 1 : 0;
        if (enforce_monotonic) {
            if (seen_reject) accepted = 0;
            if (!accepted) seen_reject = true;
        }

        out.q_obs_scan.push_back(q_obs);
        out.p_value_scan.push_back(p_value);
        out.accepted_scan.push_back(accepted);

        if (accepted) {
            out.BR_90 = calib_points[i].BR_test;
            out.N_sig_90_nominal = calib_points[i].N_sig_test_nominal;
        }
    }

    out.ok = true;
    return out;
}

static LimitSummary EvaluateFastLimitWithCalibration(
    const std::vector<Event>& events,
    const std::vector<PdfComponent>& components,
    const std::vector<CalibrationPoint>& calib_points,
    const std::vector<double>& n_sig_scan,
    const FitConfig& free_cfg,
    const FitConfig& prof_cfg,
    double cl)
{
    const ProfileLikelihoodQScanResult qscan =
        EvaluateProfileLikelihoodQScan(events, components, free_cfg, prof_cfg, n_sig_scan);
    return EvaluateLimitFromCalibrationQScan(qscan, calib_points, cl, true);
}

static bool BuildMeanYieldsByPointFromQScan(
    const ProfileLikelihoodQScanResult& qscan,
    const std::vector<double>& n_sig_scan,
    const std::vector<double>& fallback_yields,
    std::vector<std::vector<double>>& out_mean_yields)
{
    out_mean_yields.clear();
    if (qscan.points.size() != n_sig_scan.size()) return false;

    out_mean_yields.reserve(n_sig_scan.size());
    for (std::size_t i = 0; i < n_sig_scan.size(); ++i) {
        std::vector<double> mean_yields = fallback_yields;
        const FitResult& fit_prof = qscan.points[i].fit_prof;
        if (fit_prof.status == 0 && fit_prof.yields_hat.size() == fallback_yields.size()) {
            mean_yields = fit_prof.yields_hat;
        }
        if (mean_yields.size() != fallback_yields.size()) return false;
        if (!mean_yields.empty()) mean_yields[0] = n_sig_scan[i];
        out_mean_yields.push_back(mean_yields);
    }
    return true;
}

static int LastAcceptedIndex(const std::vector<int>& accepted_scan)
{
    int idx = -1;
    for (std::size_t i = 0; i < accepted_scan.size(); ++i) {
        if (accepted_scan[i]) idx = static_cast<int>(i);
    }
    return idx;
}

static LimitSummary EvaluateLocallyExactLimit(
    const std::vector<Event>& events,
    const std::vector<PdfComponent>& components,
    const std::vector<double>& br_scan,
    const std::vector<double>& n_sig_scan,
    const std::vector<CalibrationPoint>& calib_points_global,
    const FitConfig& free_cfg,
    const FitConfig& prof_cfg,
    const ToyGeneratorConfig& toy_cfg,
    const NormalizationUncertaintyConfig& norm_cfg,
    int n_local_exact_toys,
    int local_exact_half_window,
    double cl)
{
    LimitSummary out;
    out.BR_90 = 0.0;
    out.N_sig_90_nominal = 0.0;
    out.ok = false;

    const ProfileLikelihoodQScanResult qscan =
        EvaluateProfileLikelihoodQScan(events, components, free_cfg, prof_cfg, n_sig_scan);
    if (qscan.fit_free.status != 0 || qscan.points.size() != br_scan.size()) {
        return out;
    }

    const LimitSummary fast_raw =
        EvaluateLimitFromCalibrationQScan(qscan, calib_points_global, cl, false);
    if (!fast_raw.ok) return out;

    out = fast_raw;

    const int last_acc = LastAcceptedIndex(fast_raw.accepted_scan);
    if (last_acc < 0 || n_local_exact_toys <= 0 || local_exact_half_window < 0) {
        out.ok = true;
        return out;
    }

    const int i_min = std::max(0, last_acc - local_exact_half_window);
    const int i_max = std::min(static_cast<int>(br_scan.size()) - 1,
                               last_acc + local_exact_half_window + 1);

    std::vector<std::vector<double>> mean_yields_by_point;
    std::vector<double> fallback_yields = qscan.fit_free.yields_hat;
    if (fallback_yields.size() != components.size()) {
        fallback_yields.assign(components.size(), 0.0);
    }
    if (!BuildMeanYieldsByPointFromQScan(qscan, n_sig_scan, fallback_yields,
                                         mean_yields_by_point)) {
        return out;
    }

    std::vector<double> br_subset;
    std::vector<double> nsig_subset;
    std::vector<std::vector<double>> mean_subset;
    br_subset.reserve(static_cast<std::size_t>(i_max - i_min + 1));
    nsig_subset.reserve(static_cast<std::size_t>(i_max - i_min + 1));
    mean_subset.reserve(static_cast<std::size_t>(i_max - i_min + 1));
    for (int i = i_min; i <= i_max; ++i) {
        br_subset.push_back(br_scan[static_cast<std::size_t>(i)]);
        nsig_subset.push_back(n_sig_scan[static_cast<std::size_t>(i)]);
        mean_subset.push_back(mean_yields_by_point[static_cast<std::size_t>(i)]);
    }

    std::vector<double> pmax_cache_local(components.size(), 0.0);
    std::vector<CalibrationPoint> calib_subset;
    if (!BuildCalibrationDistributions(components, br_subset, nsig_subset,
                                       mean_subset, free_cfg, prof_cfg,
                                       toy_cfg, norm_cfg, n_local_exact_toys,
                                       pmax_cache_local, calib_subset)) {
        return out;
    }

    const ProfileLikelihoodQScanResult qscan_subset =
        EvaluateProfileLikelihoodQScan(events, components, free_cfg, prof_cfg, nsig_subset);
    const LimitSummary local_exact =
        EvaluateLimitFromCalibrationQScan(qscan_subset, calib_subset, cl, false);
    if (!local_exact.ok) {
        out.ok = true;
        return out;
    }

    for (int i = i_min; i <= i_max; ++i) {
        const std::size_t j = static_cast<std::size_t>(i - i_min);
        out.q_obs_scan[static_cast<std::size_t>(i)] = local_exact.q_obs_scan[j];
        out.p_value_scan[static_cast<std::size_t>(i)] = local_exact.p_value_scan[j];
        out.accepted_scan[static_cast<std::size_t>(i)] = local_exact.accepted_scan[j];
    }

    out.BR_90 = 0.0;
    out.N_sig_90_nominal = 0.0;
    for (std::size_t i = 0; i < br_scan.size(); ++i) {
        if (!out.accepted_scan[i]) continue;
        out.BR_90 = br_scan[i];
        out.N_sig_90_nominal = n_sig_scan[i];
    }

    out.ok = true;
    return out;
}

static FitConfig BuildDefaultFitConfig(double N0)
{
    FitConfig cfg;
    cfg.start_yields = {N0 / 3.0, N0 / 3.0, N0 / 3.0};
    cfg.max_calls = 20000;
    cfg.tol = 1e-3;
    return cfg;
}

int main(int argc, char** argv)
{
    if (argc != 9 && argc != 10 && argc != 11 && argc != 12) {
        std::cerr
            << "Usage: " << argv[0]
            << " <datafile> <N_mu_eff_nom> <N_mu_eff_sigma>"
            << " <br_min> <br_max> <br_step>"
            << " <n_calib_toys> <n_sens_toys>"
            << " [seed] [n_validation_toys] [n_validation_inner]\n";
        return 1;
    }

    const char* datafile = argv[1];

    double N_mu_eff_nom = 0.0;
    double N_mu_eff_sigma = 0.0;
    double br_min = 0.0;
    double br_max = 0.0;
    double br_step = 0.0;
    int n_calib_toys = 0;
    int n_sens_toys = 0;
    unsigned long long seed = 123456789ULL;
    int n_validation_toys = 0;
    int n_validation_inner = 0;

    if (!ParseDoubleStrict(argv[2], N_mu_eff_nom) ||
        !ParseDoubleStrict(argv[3], N_mu_eff_sigma) ||
        !ParseDoubleStrict(argv[4], br_min) ||
        !ParseDoubleStrict(argv[5], br_max) ||
        !ParseDoubleStrict(argv[6], br_step) ||
        !ParseIntStrict(argv[7], n_calib_toys) ||
        !ParseIntStrict(argv[8], n_sens_toys)) {
        std::cerr << "[run_sensitivity_br_toymc] failed to parse required arguments\n";
        return 1;
    }
    if (argc >= 10 && !ParseULLStrict(argv[9], seed)) {
        std::cerr << "[run_sensitivity_br_toymc] failed to parse seed\n";
        return 1;
    }
    if (argc >= 11 && !ParseIntStrict(argv[10], n_validation_toys)) {
        std::cerr << "[run_sensitivity_br_toymc] failed to parse n_validation_toys\n";
        return 1;
    }
    if (argc >= 12 && !ParseIntStrict(argv[11], n_validation_inner)) {
        std::cerr << "[run_sensitivity_br_toymc] failed to parse n_validation_inner\n";
        return 1;
    }

    const std::vector<double> br_scan = BuildScanPoints(br_min, br_max, br_step);
    if (br_scan.empty()) {
        std::cerr << "[run_sensitivity_br_toymc] invalid BR scan range\n";
        return 1;
    }
    if (!(N_mu_eff_nom > 0.0) || !std::isfinite(N_mu_eff_nom)) {
        std::cerr << "[run_sensitivity_br_toymc] N_mu_eff_nom must be positive\n";
        return 1;
    }
    if (!(N_mu_eff_sigma >= 0.0) || !std::isfinite(N_mu_eff_sigma)) {
        std::cerr << "[run_sensitivity_br_toymc] N_mu_eff_sigma must be non-negative\n";
        return 1;
    }
    if (n_calib_toys <= 0 || n_sens_toys <= 0) {
        std::cerr << "[run_sensitivity_br_toymc] toy counts must be positive\n";
        return 1;
    }

    std::vector<Event> data_events;
    long n_skipped_non5 = 0;
    long n_skipped_outside = 0;
    if (!LoadEventsFromDat(datafile, data_events, n_skipped_non5, n_skipped_outside)) {
        std::cerr << "[run_sensitivity_br_toymc] no valid events in analysis window: "
                  << datafile << "\n";
        return 2;
    }

    if (!RMDGridPdf_Load(kDefaultRmdRoot, kDefaultRmdKey)) {
        std::cerr << "[run_sensitivity_br_toymc] RMDGridPdf_Load failed\n";
        return 3;
    }
    if (!ACCGridPdf_Load(kDefaultAccRoot, kDefaultAccKey)) {
        std::cerr << "[run_sensitivity_br_toymc] ACCGridPdf_Load failed\n";
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
    const FitConfig free_cfg = BuildDefaultFitConfig(N0);
    FitConfig prof_cfg = free_cfg;

    const FitResult fit_data_free = FitNLL(data_events, components, free_cfg);

    FitConfig bg_cfg = free_cfg;
    if (bg_cfg.start_yields.size() >= 1U) bg_cfg.start_yields[0] = 0.0;
    const FitResult fit_data_bg0 =
        FitNLLFixedSignal(data_events, components, bg_cfg, 0.0);
    if (fit_data_bg0.status != 0 || fit_data_bg0.yields_hat.size() != components.size()) {
        std::cerr << "[run_sensitivity_br_toymc] background-only fit failed\n";
        return 4;
    }

    std::vector<double> n_sig_scan;
    n_sig_scan.reserve(br_scan.size());
    for (double br : br_scan) n_sig_scan.push_back(br * N_mu_eff_nom);

    const ProfileLikelihoodQScanResult qscan_data_nominal =
        EvaluateProfileLikelihoodQScan(data_events, components, free_cfg, prof_cfg, n_sig_scan);
    std::vector<std::vector<double>> calib_mean_yields_by_point;
    calib_mean_yields_by_point.reserve(br_scan.size());
    for (std::size_t i = 0; i < br_scan.size(); ++i) {
        std::vector<double> mean_yields = fit_data_bg0.yields_hat;
        if (i < qscan_data_nominal.points.size()) {
            const FitResult& fit_prof = qscan_data_nominal.points[i].fit_prof;
            if (fit_prof.status == 0 && fit_prof.yields_hat.size() == components.size()) {
                mean_yields = fit_prof.yields_hat;
            }
        }
        if (mean_yields.size() == components.size()) {
            mean_yields[0] = n_sig_scan[i];
        }
        calib_mean_yields_by_point.push_back(mean_yields);
    }

    ToyGeneratorConfig toy_cfg;
    toy_cfg.seed = seed;
    toy_cfg.pmax_scan_trials = 20000;
    toy_cfg.pmax_safety = 5.0;
    toy_cfg.pmax_update = 1.2;
    toy_cfg.event_pool_size_per_component = 20000;

    NormalizationUncertaintyConfig norm_cfg;
    norm_cfg.N_mu_eff_nom = N_mu_eff_nom;
    norm_cfg.N_mu_eff_sigma = N_mu_eff_sigma;

    std::vector<double> pmax_cache(components.size(), 0.0);
    std::vector<CalibrationPoint> calib_points;

    const auto t_calib_begin = std::chrono::steady_clock::now();
    if (!BuildCalibrationDistributions(components, br_scan, n_sig_scan,
                                       calib_mean_yields_by_point,
                                       free_cfg, prof_cfg,
                                       toy_cfg, norm_cfg,
                                       n_calib_toys,
                                       pmax_cache, calib_points)) {
        std::cerr << "[run_sensitivity_br_toymc] calibration failed\n";
        return 5;
    }
    const auto t_calib_end = std::chrono::steady_clock::now();

    const LimitSummary obs_fast = EvaluateFastLimitWithCalibration(
        data_events, components, calib_points, n_sig_scan, free_cfg, prof_cfg, 0.90);
    const LimitSummary obs_local_exact = EvaluateLocallyExactLimit(
        data_events, components, br_scan, n_sig_scan,
        calib_points, free_cfg, prof_cfg,
        toy_cfg, norm_cfg,
        n_calib_toys, kLocalExactHalfWindow, 0.90);

    std::vector<double> sens_br_values;
    std::vector<double> sens_nsig_values;
    sens_br_values.reserve(static_cast<std::size_t>(n_sens_toys));
    sens_nsig_values.reserve(static_cast<std::size_t>(n_sens_toys));

    int n_sens_valid = 0;
    for (int itoy = 0; itoy < n_sens_toys; ++itoy) {
        std::vector<Event> toy_events;
        std::vector<double> toy_generated_yields;
        if (!GenerateToyDatasetFromModel(components, fit_data_bg0.yields_hat, toy_cfg,
                                         static_cast<unsigned long long>(0xD0000000ULL +
                                         static_cast<unsigned long long>(itoy)),
                                         toy_events, &pmax_cache, &toy_generated_yields)) {
            continue;
        }
        if (toy_events.empty()) continue;

        const LimitSummary lim = EvaluateLocallyExactLimit(
            toy_events, components, br_scan, n_sig_scan,
            calib_points, free_cfg, prof_cfg,
            toy_cfg, norm_cfg,
            n_calib_toys, kLocalExactHalfWindow, 0.90);
        if (!lim.ok) continue;

        sens_br_values.push_back(lim.BR_90);
        sens_nsig_values.push_back(lim.N_sig_90_nominal);
        ++n_sens_valid;
    }

    std::sort(sens_br_values.begin(), sens_br_values.end());
    std::sort(sens_nsig_values.begin(), sens_nsig_values.end());

    const double br_sens_median = MedianFromSorted(sens_br_values);
    const double nsig_sens_median = MedianFromSorted(sens_nsig_values);

    const auto t_sens_end = std::chrono::steady_clock::now();

    double validation_fast_median = 0.0;
    double validation_exact_median = 0.0;
    double validation_mean_abs_diff = 0.0;
    int validation_valid = 0;
    double validation_fast_seconds = 0.0;
    double validation_exact_seconds = 0.0;

    if (n_validation_toys > 0 && n_validation_inner > 0) {
        std::vector<double> fast_vals;
        std::vector<double> exact_vals;
        fast_vals.reserve(static_cast<std::size_t>(n_validation_toys));
        exact_vals.reserve(static_cast<std::size_t>(n_validation_toys));

        for (int itoy = 0; itoy < n_validation_toys; ++itoy) {
            std::vector<Event> toy_events;
            std::vector<double> toy_generated_yields;
            if (!GenerateToyDatasetFromModel(components, fit_data_bg0.yields_hat, toy_cfg,
                                             static_cast<unsigned long long>(0xE0000000ULL +
                                             static_cast<unsigned long long>(itoy)),
                                             toy_events, &pmax_cache, &toy_generated_yields)) {
                continue;
            }
            if (toy_events.empty()) continue;

            const auto t_fast_begin = std::chrono::steady_clock::now();
            const LimitSummary lim_fast = EvaluateLocallyExactLimit(
                toy_events, components, br_scan, n_sig_scan,
                calib_points, free_cfg, prof_cfg,
                toy_cfg, norm_cfg,
                n_calib_toys, kLocalExactHalfWindow, 0.90);
            const auto t_fast_end = std::chrono::steady_clock::now();
            validation_fast_seconds +=
                std::chrono::duration<double>(t_fast_end - t_fast_begin).count();
            if (!lim_fast.ok) continue;

            UpperLimitBRScanConfig exact_cfg;
            exact_cfg.BR_scan = br_scan;
            exact_cfg.n_toys_per_point = n_validation_inner;
            exact_cfg.cl = 0.90;
            exact_cfg.free_fit_cfg = free_cfg;
            exact_cfg.prof_fit_cfg = prof_cfg;
            exact_cfg.toy_cfg = toy_cfg;
            exact_cfg.toy_cfg.seed =
                MixSeed(seed, static_cast<unsigned long long>(0xF0000000ULL + itoy));
            exact_cfg.norm_cfg = norm_cfg;

            const auto t_exact_begin = std::chrono::steady_clock::now();
            const UpperLimitBRScanResult lim_exact =
                EvaluateUpperLimitBRScan(toy_events, components, exact_cfg);
            const auto t_exact_end = std::chrono::steady_clock::now();
            validation_exact_seconds +=
                std::chrono::duration<double>(t_exact_end - t_exact_begin).count();

            fast_vals.push_back(lim_fast.BR_90);
            exact_vals.push_back(lim_exact.BR_90);
        }

        validation_valid = static_cast<int>(std::min(fast_vals.size(), exact_vals.size()));
        if (validation_valid > 0) {
            for (int i = 0; i < validation_valid; ++i) {
                validation_mean_abs_diff += std::fabs(fast_vals[i] - exact_vals[i]);
            }
            validation_mean_abs_diff /= static_cast<double>(validation_valid);

            std::sort(fast_vals.begin(), fast_vals.end());
            std::sort(exact_vals.begin(), exact_vals.end());
            validation_fast_median = MedianFromSorted(fast_vals);
            validation_exact_median = MedianFromSorted(exact_vals);
        }
    }

    std::cout << std::setprecision(10);
    std::cout << "==================== Input ====================\n";
    std::cout << "datafile                   = " << datafile << "\n";
    std::cout << "events_used                = " << data_events.size() << "\n";
    std::cout << "skipped_non5               = " << n_skipped_non5 << "\n";
    std::cout << "skipped_outside_win        = " << n_skipped_outside << "\n";
    std::cout << "N_mu_eff_nom               = " << N_mu_eff_nom << "\n";
    std::cout << "N_mu_eff_sigma             = " << N_mu_eff_sigma << "\n";
    std::cout << "br_scan_points             = " << br_scan.size() << "\n";
    std::cout << "n_calib_toys               = " << n_calib_toys << "\n";
    std::cout << "n_sens_toys                = " << n_sens_toys << "\n";
    std::cout << "seed                       = " << seed << "\n";
    std::cout << "n_validation_toys          = " << n_validation_toys << "\n";
    std::cout << "n_validation_inner         = " << n_validation_inner << "\n";
    std::cout << "local_exact_half_window    = " << kLocalExactHalfWindow << "\n";
    std::cout << "local_exact_inner_toys     = " << n_calib_toys << "\n";
    std::cout << "===============================================\n";

    std::cout << "==================== Fits =====================\n";
    std::cout << "data_free_status           = " << fit_data_free.status << "\n";
    std::cout << "data_free_nll             = " << fit_data_free.nll_min << "\n";
    if (fit_data_free.yields_hat.size() >= 3U) {
        std::cout << "data_free_N_sig_hat       = " << fit_data_free.yields_hat[0] << "\n";
        std::cout << "data_free_N_rmd_hat       = " << fit_data_free.yields_hat[1] << "\n";
        std::cout << "data_free_N_acc_hat       = " << fit_data_free.yields_hat[2] << "\n";
    }
    std::cout << "data_bg0_status            = " << fit_data_bg0.status << "\n";
    std::cout << "data_bg0_nll              = " << fit_data_bg0.nll_min << "\n";
    if (fit_data_bg0.yields_hat.size() >= 3U) {
        std::cout << "data_bg0_N_sig_fixed      = " << fit_data_bg0.yields_hat[0] << "\n";
        std::cout << "data_bg0_N_rmd_hat        = " << fit_data_bg0.yields_hat[1] << "\n";
        std::cout << "data_bg0_N_acc_hat        = " << fit_data_bg0.yields_hat[2] << "\n";
    }
    std::cout << "===============================================\n";

    std::cout << "================== Calibration =================\n";
    std::cout << "# BR  N_sig_nom  n_valid  q50  q90  q95\n";
    for (const auto& p : calib_points) {
        std::cout << p.BR_test << " "
                  << p.N_sig_test_nominal << " "
                  << p.n_toys_valid << " "
                  << QuantileFromSorted(p.q_sorted, 0.50) << " "
                  << QuantileFromSorted(p.q_sorted, 0.90) << " "
                  << QuantileFromSorted(p.q_sorted, 0.95) << "\n";
    }
    std::cout << "calibration_seconds        = "
              << std::chrono::duration<double>(t_calib_end - t_calib_begin).count() << "\n";
    std::cout << "===============================================\n";

    std::cout << "=================== Observed ===================\n";
    std::cout << "observed_fast_ok           = " << (obs_fast.ok ? 1 : 0) << "\n";
    std::cout << "observed_fast_BR90         = " << obs_fast.BR_90 << "\n";
    std::cout << "observed_fast_Nsig90_nom   = " << obs_fast.N_sig_90_nominal << "\n";
    std::cout << "observed_local_exact_ok    = " << (obs_local_exact.ok ? 1 : 0) << "\n";
    std::cout << "observed_local_exact_BR90  = " << obs_local_exact.BR_90 << "\n";
    std::cout << "observed_local_exact_Nsig90_nom = " << obs_local_exact.N_sig_90_nominal << "\n";
    std::cout << "# BR  q_obs  p_value  accepted\n";
    for (std::size_t i = 0; i < br_scan.size(); ++i) {
        const double q_obs =
            (i < obs_local_exact.q_obs_scan.size()) ? obs_local_exact.q_obs_scan[i] : 0.0;
        const double p_value =
            (i < obs_local_exact.p_value_scan.size()) ? obs_local_exact.p_value_scan[i] : 0.0;
        const int accepted =
            (i < obs_local_exact.accepted_scan.size()) ? obs_local_exact.accepted_scan[i] : 0;
        std::cout << br_scan[i] << " "
                  << q_obs << " "
                  << p_value << " "
                  << accepted << "\n";
    }
    std::cout << "===============================================\n";

    std::cout << "================= Sensitivity ==================\n";
    std::cout << "sensitivity_valid_toys     = " << n_sens_valid << "\n";
    std::cout << "sensitivity_BR50           = " << br_sens_median << "\n";
    std::cout << "sensitivity_Nsig50_nom     = " << nsig_sens_median << "\n";
    std::cout << "sensitivity_BR16           = " << QuantileFromSorted(sens_br_values, 0.16) << "\n";
    std::cout << "sensitivity_BR84           = " << QuantileFromSorted(sens_br_values, 0.84) << "\n";
    std::cout << "sensitivity_BR05           = " << QuantileFromSorted(sens_br_values, 0.05) << "\n";
    std::cout << "sensitivity_BR95           = " << QuantileFromSorted(sens_br_values, 0.95) << "\n";
    std::cout << "sensitivity_seconds_total  = "
              << std::chrono::duration<double>(t_sens_end - t_calib_end).count() << "\n";
    std::cout << "===============================================\n";

    std::cout << "=============== Sensitivity Toys ==============\n";
    std::cout << "# toy_index  BR90\n";
    for (std::size_t i = 0; i < sens_br_values.size(); ++i) {
        std::cout << i << " "
                  << sens_br_values[i] << "\n";
    }
    std::cout << "===============================================\n";

    std::cout << "================== Validation ==================\n";
    std::cout << "validation_valid           = " << validation_valid << "\n";
    std::cout << "validation_fast_median     = " << validation_fast_median << "\n";
    std::cout << "validation_exact_median    = " << validation_exact_median << "\n";
    std::cout << "validation_mean_abs_diff   = " << validation_mean_abs_diff << "\n";
    std::cout << "validation_fast_seconds    = " << validation_fast_seconds << "\n";
    std::cout << "validation_exact_seconds   = " << validation_exact_seconds << "\n";
    std::cout << "===============================================\n";

    std::cout << "note: sensitivity now uses a local exact refinement around the FC boundary.\n";
    std::cout << "      First a global calibration belt is built once, then each outer toy\n";
    std::cout << "      re-calibrates only the boundary-neighbour BR points with toy-specific\n";
    std::cout << "      nested FC toys. validation compares this local-exact method against\n";
    std::cout << "      the slower full nested toy construction.\n";

    return 0;
}
