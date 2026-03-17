// scripts/run_upper_limit_br_toymc.cc
//
// 使い方:
//   ./build/run_upper_limit_br_toymc
//       <datafile> <N_mu_eff_nom> <N_mu_eff_sigma>
//       <br_min> <br_max> <br_step> <n_toys> [seed]
//
// 出力:
//  - 実データの best fit yields
//  - 各 BR 点の q_obs(BR), toy p-value, 受容判定
//  - 90% C.L. upper limit BR_90
//  - 公称 N_mu_eff に対応する N_sig_90 = BR_90 * N_mu_eff_nom
//
// 注意:
//  - toy 生成時のみ N_mu_eff を Gaussian で揺らして
//    normalisation uncertainty を入れる
//  - likelihood 自体には normalisation nuisance を明示的に追加していない
//  - それ以外の制約は現在の ConstraintNLL(yields) のみを使う

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
        std::cerr << "[run_upper_limit_br_toymc] cannot open: " << filepath << "\n";
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

int main(int argc, char** argv)
{
    if (argc != 8 && argc != 9) {
        std::cerr
            << "Usage: " << argv[0]
            << " <datafile> <N_mu_eff_nom> <N_mu_eff_sigma>"
            << " <br_min> <br_max> <br_step> <n_toys> [seed]\n";
        return 1;
    }

    const char* datafile = argv[1];

    double N_mu_eff_nom = 0.0;
    double N_mu_eff_sigma = 0.0;
    double br_min = 0.0;
    double br_max = 0.0;
    double br_step = 0.0;
    int n_toys = 0;
    unsigned long long seed = 123456789ULL;

    if (!ParseDoubleStrict(argv[2], N_mu_eff_nom) ||
        !ParseDoubleStrict(argv[3], N_mu_eff_sigma) ||
        !ParseDoubleStrict(argv[4], br_min) ||
        !ParseDoubleStrict(argv[5], br_max) ||
        !ParseDoubleStrict(argv[6], br_step) ||
        !ParseIntStrict(argv[7], n_toys)) {
        std::cerr << "[run_upper_limit_br_toymc] failed to parse arguments\n";
        return 1;
    }
    if (argc == 9 && !ParseULLStrict(argv[8], seed)) {
        std::cerr << "[run_upper_limit_br_toymc] failed to parse seed\n";
        return 1;
    }

    std::vector<double> scan_points = BuildScanPoints(br_min, br_max, br_step);
    if (scan_points.empty()) {
        std::cerr << "[run_upper_limit_br_toymc] invalid scan range\n";
        return 1;
    }
    if (!(N_mu_eff_nom > 0.0) || !std::isfinite(N_mu_eff_nom)) {
        std::cerr << "[run_upper_limit_br_toymc] N_mu_eff_nom must be positive\n";
        return 1;
    }
    if (!(N_mu_eff_sigma >= 0.0) || !std::isfinite(N_mu_eff_sigma)) {
        std::cerr << "[run_upper_limit_br_toymc] N_mu_eff_sigma must be non-negative\n";
        return 1;
    }
    if (n_toys <= 0) {
        std::cerr << "[run_upper_limit_br_toymc] n_toys must be positive\n";
        return 1;
    }

    std::vector<Event> events;
    long n_skipped_non5 = 0;
    long n_skipped_outside = 0;
    if (!LoadEventsFromDat(datafile, events, n_skipped_non5, n_skipped_outside)) {
        std::cerr << "[run_upper_limit_br_toymc] no valid events in analysis window: "
                  << datafile << "\n";
        return 2;
    }

    if (!RMDGridPdf_Load(kDefaultRmdRoot, kDefaultRmdKey)) {
        std::cerr << "[run_upper_limit_br_toymc] RMDGridPdf_Load failed\n";
        return 3;
    }
    if (!ACCGridPdf_Load(kDefaultAccRoot, kDefaultAccKey)) {
        std::cerr << "[run_upper_limit_br_toymc] ACCGridPdf_Load failed\n";
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

    UpperLimitBRScanConfig scan_cfg;
    scan_cfg.BR_scan = scan_points;
    scan_cfg.n_toys_per_point = n_toys;
    scan_cfg.cl = 0.90;
    scan_cfg.free_fit_cfg = free_cfg;
    scan_cfg.prof_fit_cfg = prof_cfg;
    scan_cfg.toy_cfg = toy_cfg;
    scan_cfg.norm_cfg.N_mu_eff_nom = N_mu_eff_nom;
    scan_cfg.norm_cfg.N_mu_eff_sigma = N_mu_eff_sigma;

    const UpperLimitBRScanResult res =
        EvaluateUpperLimitBRScan(events, components, scan_cfg);

    std::cout << std::setprecision(10);
    std::cout << "==================== Input ====================\n";
    std::cout << "datafile             = " << datafile << "\n";
    std::cout << "events_used          = " << events.size() << "\n";
    std::cout << "skipped_non5         = " << n_skipped_non5 << "\n";
    std::cout << "skipped_outside_win  = " << n_skipped_outside << "\n";
    std::cout << "N_mu_eff_nom         = " << N_mu_eff_nom << "\n";
    std::cout << "N_mu_eff_sigma       = " << N_mu_eff_sigma << "\n";
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
    std::cout << "# BR  N_sig_nom  q_obs  p_value  accepted  valid_toys  N_rmd_prof  N_acc_prof\n";
    for (const auto& p : res.points) {
        double n_rmd_prof = 0.0;
        double n_acc_prof = 0.0;
        if (p.fit_prof_obs.yields_hat.size() >= 3) {
            n_rmd_prof = p.fit_prof_obs.yields_hat[1];
            n_acc_prof = p.fit_prof_obs.yields_hat[2];
        }
        std::cout << p.BR_test << " "
                  << p.N_sig_test_nominal << " "
                  << p.q_obs << " "
                  << p.p_value << " "
                  << (p.accepted ? 1 : 0) << " "
                  << p.n_toys_valid << " "
                  << n_rmd_prof << " "
                  << n_acc_prof << "\n";
    }
    std::cout << "===============================================\n";

    std::cout << "==================== Limit ====================\n";
    std::cout << "BR_90             = " << res.BR_90 << "\n";
    std::cout << "N_sig_90_nominal  = " << res.N_sig_90_nominal << "\n";
    std::cout << "===============================================\n";

    std::cout << "note: current BR upper limit randomises only N_mu_eff in toy generation\n";
    std::cout << "      and otherwise uses the existing ConstraintNLL(yields).\n";

    return 0;
}
