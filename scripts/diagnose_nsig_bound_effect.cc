// scripts/diagnose_nsig_bound_effect.cc
//
// 目的:
//  - N_sig >= 0 制約あり/なしで free fit の挙動を比較する
//  - 開始値依存性を確認する
//  - 代表的な fixed N_sig に対して q_obs が定義できるか確認する
//  - 背景を固定したときの NLL(N_sig) を調べ、負の N_sig を好むかを見る

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
#include "p2meg/Likelihood.h"
#include "p2meg/MathUtils.h"
#include "p2meg/NLLFit.h"
#include "p2meg/PdfWrappers.h"
#include "p2meg/RMDGridPdf.h"
#include "p2meg/UpperLimit.h"

static const char* kDefaultRmdRoot = "data/pdf_cache/rmd_grid.root";
static const char* kDefaultRmdKey  = "rmd_grid";
static const char* kDefaultAccRoot = "data/pdf_cache/acc_grid.root";
static const char* kDefaultAccKey  = "acc_grid";

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

static void SetYieldLowerBoundsEnabled(bool enable)
{
    if (enable) {
        unsetenv("P2MEG_DISABLE_YIELD_BOUNDS");
    } else {
        setenv("P2MEG_DISABLE_YIELD_BOUNDS", "1", 1);
    }
}

static const char* BoundModeLabel(bool enable)
{
    return enable ? "bounded" : "unbounded";
}

static void PrintFitLine(const char* tag, const FitResult& fit)
{
    std::cout << std::setw(14) << tag
              << " status=" << std::setw(2) << fit.status
              << " nll=" << std::setw(14) << fit.nll_min;
    if (fit.yields_hat.size() >= 3) {
        std::cout << " Nsig=" << std::setw(12) << fit.yields_hat[0]
                  << " Nrmd=" << std::setw(12) << fit.yields_hat[1]
                  << " Nacc=" << std::setw(12) << fit.yields_hat[2];
    }
    std::cout << "\n";
}

static bool EvaluateEventIntensity(const Event& ev,
                                   const std::vector<PdfComponent>& components,
                                   const std::vector<double>& yields,
                                   double& pi_out,
                                   double& sig_pdf_out)
{
    pi_out = 0.0;
    sig_pdf_out = 0.0;
    if (components.size() != yields.size()) return false;
    if (components.empty()) return false;

    for (std::size_t k = 0; k < components.size(); ++k) {
        const auto& comp = components[k];
        double pk = comp.eval ? comp.eval(ev, comp.ctx) : 0.0;
        if (!std::isfinite(pk) || pk < 0.0) pk = 0.0;
        if (k == 0) sig_pdf_out = pk;
        pi_out += yields[k] * pk;
    }

    return std::isfinite(pi_out) && std::isfinite(sig_pdf_out);
}

static void PrintFixedBackgroundScan(const std::vector<Event>& events,
                                     const std::vector<PdfComponent>& components,
                                     const FitResult& ref_fit,
                                     double n_sig_min,
                                     double n_sig_max,
                                     double n_sig_step)
{
    if (ref_fit.yields_hat.size() < 3) return;

    const double N_rmd = ref_fit.yields_hat[1];
    const double N_acc = ref_fit.yields_hat[2];

    std::cout << "# fixed-background NLL scan"
              << " (N_rmd=" << N_rmd
              << ", N_acc=" << N_acc << ")\n";
    std::cout << "# N_sig  NLL  min_pi  dNLL_dNsig_at_point_like\n";

    for (double N_sig = n_sig_min; N_sig <= n_sig_max + 0.5 * n_sig_step; N_sig += n_sig_step) {
        std::vector<double> yields{N_sig, N_rmd, N_acc};
        const double nll = NLL(events, components, yields);

        double min_pi = 1e300;
        double score_sum = 0.0;
        bool ok = true;
        for (const auto& ev : events) {
            double pi = 0.0;
            double p_sig = 0.0;
            if (!EvaluateEventIntensity(ev, components, yields, pi, p_sig)) {
                ok = false;
                break;
            }
            if (!(pi > 0.0) || !std::isfinite(pi)) {
                ok = false;
                break;
            }
            if (pi < min_pi) min_pi = pi;
            score_sum += p_sig / pi;
        }

        double deriv_like = std::numeric_limits<double>::quiet_NaN();
        if (ok) {
            // dNLL/dNsig = 1 - Σ_i p_sig(x_i)/pi_i + N_sig/sigma^2
            // ここで sigma=1000 は現状の弱い Gaussian constraint。
            deriv_like = 1.0 - score_sum + (N_sig / (1000.0 * 1000.0));
        }

        std::cout << N_sig << " "
                  << nll << " "
                  << (ok ? min_pi : std::numeric_limits<double>::quiet_NaN()) << " "
                  << deriv_like << "\n";
    }
}

static void RunMode(const char* datafile,
                    const std::vector<Event>& events,
                    const std::vector<PdfComponent>& components,
                    const std::vector<std::pair<std::string, std::vector<double>>>& starts,
                    const std::vector<double>& fixed_nsig_tests)
{
    const double N0 = static_cast<double>(events.size());

    FitResult best_fit;
    best_fit.status = 999;
    best_fit.nll_min = std::numeric_limits<double>::quiet_NaN();

    std::cout << "================ mode: "
              << BoundModeLabel(std::getenv("P2MEG_DISABLE_YIELD_BOUNDS") == nullptr)
              << " ================\n";

    for (const auto& item : starts) {
        FitConfig cfg;
        cfg.start_yields = item.second;
        cfg.max_calls = 20000;
        cfg.tol = 1e-3;

        const FitResult fit = FitNLL(events, components, cfg);
        PrintFitLine(item.first.c_str(), fit);

        if (std::isfinite(fit.nll_min) &&
            (!std::isfinite(best_fit.nll_min) || fit.nll_min < best_fit.nll_min)) {
            best_fit = fit;
        }
    }

    FitConfig deep_cfg;
    deep_cfg.start_yields = {N0 / 3.0, N0 / 3.0, N0 / 3.0};
    deep_cfg.max_calls = 200000;
    deep_cfg.tol = 1e-4;
    const FitResult deep_fit = FitNLL(events, components, deep_cfg);
    PrintFitLine("equal_deep", deep_fit);
    if (std::isfinite(deep_fit.nll_min) &&
        (!std::isfinite(best_fit.nll_min) || deep_fit.nll_min < best_fit.nll_min)) {
        best_fit = deep_fit;
    }

    if (best_fit.yields_hat.size() >= 3) {
        std::cout << "# best finite fit selected for diagnostics\n";
        PrintFitLine("best", best_fit);
    }

    for (double fixed_nsig : fixed_nsig_tests) {
        FitConfig free_cfg;
        free_cfg.start_yields = {N0 / 3.0, N0 / 3.0, N0 / 3.0};
        free_cfg.max_calls = 20000;
        free_cfg.tol = 1e-3;

        FitConfig prof_cfg = free_cfg;

        FitResult fit_free;
        FitResult fit_prof;
        const double q_obs = EvaluateProfileLikelihoodQ(events, components,
                                                        free_cfg, prof_cfg,
                                                        fixed_nsig,
                                                        fit_free, fit_prof);
        std::cout << "# q_obs test for fixed N_sig = " << fixed_nsig << "\n";
        PrintFitLine("free", fit_free);
        PrintFitLine("prof", fit_prof);
        std::cout << "q_obs = " << q_obs << "\n";
    }

    if (best_fit.yields_hat.size() >= 3) {
        PrintFixedBackgroundScan(events, components, best_fit, -20.0, 5.0, 1.0);
    }

    std::cout << "===============================================\n";
    (void)datafile;
}

int main(int argc, char** argv)
{
    if (argc != 2 && argc != 3) {
        std::cerr << "Usage: " << argv[0]
                  << " <datafile> [fixed_nsig]\n";
        return 1;
    }

    const char* datafile = argv[1];
    const double fixed_nsig_user = (argc == 3) ? std::strtod(argv[2], nullptr) : 1.77664;

    std::vector<Event> events;
    long n_skipped_non5 = 0;
    long n_skipped_outside = 0;
    if (!LoadEventsFromDat(datafile, events, n_skipped_non5, n_skipped_outside)) {
        std::cerr << "[diagnose_nsig_bound_effect] no valid events in analysis window: "
                  << datafile << "\n";
        return 2;
    }

    if (!RMDGridPdf_Load(kDefaultRmdRoot, kDefaultRmdKey)) {
        std::cerr << "[diagnose_nsig_bound_effect] RMDGridPdf_Load failed\n";
        return 3;
    }
    if (!ACCGridPdf_Load(kDefaultAccRoot, kDefaultAccKey)) {
        std::cerr << "[diagnose_nsig_bound_effect] ACCGridPdf_Load failed\n";
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
    std::vector<std::pair<std::string, std::vector<double>>> starts;
    starts.push_back({"equal", {N0 / 3.0, N0 / 3.0, N0 / 3.0}});
    starts.push_back({"lowsig", {0.1, 2.0, N0 - 2.1}});
    starts.push_back({"highsig", {3.0, 3.0, N0 - 6.0}});
    starts.push_back({"accdom", {0.5, 0.5, N0 - 1.0}});
    starts.push_back({"rmddom", {0.5, 8.0, N0 - 8.5}});

    std::vector<double> fixed_nsig_tests{0.0, 1.0, fixed_nsig_user};

    std::cout << std::setprecision(10);
    std::cout << "================ input ================\n";
    std::cout << "datafile            = " << datafile << "\n";
    std::cout << "events_used         = " << events.size() << "\n";
    std::cout << "skipped_non5        = " << n_skipped_non5 << "\n";
    std::cout << "skipped_outside_win = " << n_skipped_outside << "\n";
    std::cout << "=======================================\n";

    SetYieldLowerBoundsEnabled(true);
    RunMode(datafile, events, components, starts, fixed_nsig_tests);

    SetYieldLowerBoundsEnabled(false);
    RunMode(datafile, events, components, starts, fixed_nsig_tests);

    SetYieldLowerBoundsEnabled(true);
    return 0;
}
