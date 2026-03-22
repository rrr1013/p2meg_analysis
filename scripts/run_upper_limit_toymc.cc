// scripts/run_upper_limit_toymc.cc
//
// 使い方:
//   ./build/run_upper_limit_toymc <datafile> <N_mu_eff> <s_min> <s_max> <s_step> <n_toys> [seed]
//
// 出力:
//  - 実データの best fit yields
//  - 各 s 点の q_obs(s), toy p-value, 受容判定
//  - 90% C.L. upper limit N_sig^90
//  - BR_90 = N_sig^90 / N_mu_eff
//
// 注意:
//  - 既存の final likelihood / PDF / Event をそのまま使う
//  - toy 生成は PDF から直接行う
//  - sideband / normalisation の追加 nuisance は今回入れず、
//    現在の ConstraintNLL(yields) のみを制約として使う

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
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
                              std::vector<Event>& events,
                              long& n_skipped_non5,
                              long& n_skipped_outside)
{
    events.clear();
    n_skipped_non5 = 0;
    n_skipped_outside = 0;

    std::ifstream fin(filepath);
    if (!fin) {
        std::cerr << "[run_upper_limit_toymc] cannot open: " << filepath << "\n";
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

static std::vector<double> BuildScanPoints(double s_min, double s_max, double s_step)
{
    std::vector<double> out;
    if (!std::isfinite(s_min) || !std::isfinite(s_max) || !std::isfinite(s_step)) return out;
    if (!(s_step > 0.0) || s_max < s_min) return out;

    for (double s = s_min; s <= s_max + 0.5 * s_step; s += s_step) {
        out.push_back(s);
    }
    return out;
}

int main(int argc, char** argv)
{
    if (argc != 7 && argc != 8) {
        std::cerr
            << "Usage: " << argv[0]
            << " <datafile> <N_mu_eff> <s_min> <s_max> <s_step> <n_toys> [seed]\n";
        return 1;
    }

    const char* datafile = argv[1];

    double N_mu_eff = 0.0;
    double s_min = 0.0;
    double s_max = 0.0;
    double s_step = 0.0;
    int n_toys = 0;
    unsigned long long seed = 123456789ULL;

    if (!ParseDoubleStrict(argv[2], N_mu_eff) ||
        !ParseDoubleStrict(argv[3], s_min) ||
        !ParseDoubleStrict(argv[4], s_max) ||
        !ParseDoubleStrict(argv[5], s_step) ||
        !ParseIntStrict(argv[6], n_toys)) {
        std::cerr << "[run_upper_limit_toymc] failed to parse arguments\n";
        return 1;
    }
    if (argc == 8 && !ParseULLStrict(argv[7], seed)) {
        std::cerr << "[run_upper_limit_toymc] failed to parse seed\n";
        return 1;
    }

    std::vector<double> scan_points = BuildScanPoints(s_min, s_max, s_step);
    if (scan_points.empty()) {
        std::cerr << "[run_upper_limit_toymc] invalid scan range\n";
        return 1;
    }
    if (n_toys <= 0) {
        std::cerr << "[run_upper_limit_toymc] n_toys must be positive\n";
        return 1;
    }

    std::vector<Event> events;
    long n_skipped_non5 = 0;
    long n_skipped_outside = 0;
    if (!LoadEventsFromDat(datafile, events, n_skipped_non5, n_skipped_outside)) {
        std::cerr << "[run_upper_limit_toymc] no valid events in analysis window: "
                  << datafile << "\n";
        return 2;
    }

    if (!RMDGridPdf_Load(kDefaultRmdRoot, kDefaultRmdKey)) {
        std::cerr << "[run_upper_limit_toymc] RMDGridPdf_Load failed\n";
        return 3;
    }
    if (!ACCGridPdf_Load(kDefaultAccRoot, kDefaultAccKey)) {
        std::cerr << "[run_upper_limit_toymc] ACCGridPdf_Load failed\n";
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

    ToyGeneratorConfig toy_cfg;
    toy_cfg.seed = seed;
    toy_cfg.pmax_scan_trials = 20000;
    toy_cfg.pmax_safety = 5.0;
    toy_cfg.pmax_update = 1.2;
    toy_cfg.event_pool_size_per_component = 20000;

    UpperLimitScanConfig scan_cfg;
    scan_cfg.N_sig_scan = scan_points;
    scan_cfg.n_toys_per_point = n_toys;
    scan_cfg.cl = 0.90;
    scan_cfg.N_mu_eff = N_mu_eff;
    scan_cfg.free_fit_cfg = free_cfg;
    scan_cfg.prof_fit_cfg = prof_cfg;
    scan_cfg.toy_cfg = toy_cfg;

    const UpperLimitScanResult res =
        EvaluateUpperLimitScan(events, components, scan_cfg);

    std::cout << std::setprecision(10);
    std::cout << "==================== Input ====================\n";
    std::cout << "datafile             = " << datafile << "\n";
    std::cout << "events_used          = " << events.size() << "\n";
    std::cout << "skipped_non5         = " << n_skipped_non5 << "\n";
    std::cout << "skipped_outside_win  = " << n_skipped_outside << "\n";
    std::cout << "N_mu_eff             = " << N_mu_eff << "\n";
    std::cout << "scan_points          = " << scan_points.size() << "\n";
    std::cout << "n_toys_per_point     = " << n_toys << "\n";
    std::cout << "seed                 = " << seed << "\n";
    std::cout << "===============================================\n";

    std::cout << "==================== Best Fit =================\n";
    std::cout << "status   = " << res.fit_free_obs.status << "\n";
    std::cout << "nll_min  = " << res.fit_free_obs.nll_min << "\n";
    if (res.fit_free_obs.yields_hat.size() >= 3) {
        std::cout << "N_sig_hat = " << res.fit_free_obs.yields_hat[0] << "\n";
        std::cout << "N_rmd_hat = " << res.fit_free_obs.yields_hat[1] << "\n";
        std::cout << "N_acc_hat = " << res.fit_free_obs.yields_hat[2] << "\n";
    }
    std::cout << "===============================================\n";

    std::cout << "==================== Scan =====================\n";
    std::cout << "# s  q_obs  p_value  accepted  valid_toys  N_rmd_prof  N_acc_prof\n";
    for (const auto& p : res.points) {
        double n_rmd_prof = 0.0;
        double n_acc_prof = 0.0;
        if (p.fit_prof_obs.yields_hat.size() >= 3) {
            n_rmd_prof = p.fit_prof_obs.yields_hat[1];
            n_acc_prof = p.fit_prof_obs.yields_hat[2];
        }
        std::cout << p.N_sig_test << " "
                  << p.q_obs << " "
                  << p.p_value << " "
                  << (p.accepted ? 1 : 0) << " "
                  << p.n_toys_valid << " "
                  << n_rmd_prof << " "
                  << n_acc_prof << "\n";
    }
    std::cout << "===============================================\n";

    std::cout << "==================== Limit ====================\n";
    std::cout << "N_sig_90 = " << res.N_sig_90 << "\n";
    std::cout << "BR_90    = " << res.BR_90 << "\n";
    std::cout << "===============================================\n";

    std::cout << "note: current upper limit includes only the existing ConstraintNLL(yields)\n";
    std::cout << "      and does not add new sideband / normalisation nuisance parameters.\n";

    return 0;
}
