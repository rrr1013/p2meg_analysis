#include "p2meg/UpperLimit.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <numeric>
#include <random>
#include <vector>

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"

#include "p2meg/AnalysisWindow.h"
#include "p2meg/DetectorResolution.h"
#include "p2meg/MathUtils.h"

// ============================================================
// UpperLimit 実装
//
// 実装方針:
//  - 既存の NLL / PdfComponent / FitNLL を壊さずに流用する
//  - 固定 N_sig fit だけをこのファイル内で追加する
//  - toy 生成は PDF から直接サンプリングする
//
// 注意:
//  - sideband / normalisation の追加 nuisance は今回新設しない
//  - 現在の ConstraintNLL(yields) が入っていれば、それをそのまま使う
// ============================================================

struct AllowedPhiCell {
    double phi_e;  // [rad]
    double phi_g;  // [rad]
    double area;   // [rad^2] proposal の面積重み
};

static unsigned long long MixSeed(unsigned long long seed, unsigned long long salt)
{
    unsigned long long z = seed + 0x9e3779b97f4a7c15ULL + salt;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    z = z ^ (z >> 31);
    return z;
}

static FitResult MakeFailedFitResult(std::size_t npar)
{
    FitResult out;
    out.status = 1;
    out.nll_min = std::numeric_limits<double>::quiet_NaN();
    out.yields_hat.assign(npar, 0.0);
    out.yields_err.clear();
    return out;
}

static bool BuildAllowedPhiCells(std::vector<AllowedPhiCell>& cells,
                                 double& total_area)
{
    cells.clear();
    total_area = 0.0;

    const int N_phi_e = Math_GetNPhiE(detres);
    const int N_phi_g = Math_GetNPhiG(detres);
    if (!(N_phi_e >= 1) || !(N_phi_g >= 1)) return false;

    for (int ie = 0; ie <= N_phi_e; ++ie) {
        const double phi_e =
            Detector_PhiGridPoint(ie, detres.phi_e_min, detres.phi_e_max, N_phi_e);
        const double w_e =
            Detector_PhiBinWidth(ie, detres.phi_e_min, detres.phi_e_max, N_phi_e);
        if (!(w_e > 0.0) || !std::isfinite(w_e)) continue;

        for (int ig = 0; ig <= N_phi_g; ++ig) {
            if (!Detector_IsAllowedPhiPairIndex(ie, ig, detres)) continue;

            const double phi_g =
                Detector_PhiGridPoint(ig, detres.phi_g_min, detres.phi_g_max, N_phi_g);
            const double w_g =
                Detector_PhiBinWidth(ig, detres.phi_g_min, detres.phi_g_max, N_phi_g);
            if (!(w_g > 0.0) || !std::isfinite(w_g)) continue;

            const double theta = std::fabs(phi_e - phi_g);
            if (!(theta >= analysis_window.theta_min &&
                  theta <= analysis_window.theta_max)) {
                continue;
            }

            AllowedPhiCell cell;
            cell.phi_e = phi_e;
            cell.phi_g = phi_g;
            cell.area = w_e * w_g;
            if (!(cell.area > 0.0) || !std::isfinite(cell.area)) continue;

            cells.push_back(cell);
            total_area += cell.area;
        }
    }

    return (!cells.empty() && total_area > 0.0 && std::isfinite(total_area));
}

static bool ProposeUniformEvent(std::mt19937_64& rng,
                                const std::vector<AllowedPhiCell>& cells,
                                Event& ev_out)
{
    if (cells.empty()) return false;

    const double dEe = analysis_window.Ee_max - analysis_window.Ee_min;
    const double dEg = analysis_window.Eg_max - analysis_window.Eg_min;
    const double dt  = analysis_window.t_max  - analysis_window.t_min;
    if (!(dEe > 0.0) || !(dEg > 0.0) || !(dt > 0.0)) return false;

    std::vector<double> weights;
    weights.reserve(cells.size());
    for (const auto& cell : cells) weights.push_back(cell.area);

    std::uniform_real_distribution<double> uEe(analysis_window.Ee_min, analysis_window.Ee_max);
    std::uniform_real_distribution<double> uEg(analysis_window.Eg_min, analysis_window.Eg_max);
    std::uniform_real_distribution<double> uT(analysis_window.t_min, analysis_window.t_max);
    std::discrete_distribution<std::size_t> uCell(weights.begin(), weights.end());

    const std::size_t idx = uCell(rng);
    if (idx >= cells.size()) return false;

    ev_out.Ee = uEe(rng);
    ev_out.Eg = uEg(rng);
    ev_out.t  = uT(rng);
    ev_out.phi_detector_e = cells[idx].phi_e;
    ev_out.phi_detector_g = cells[idx].phi_g;
    return true;
}

static double EstimatePMax(const PdfComponent& component,
                           std::mt19937_64& rng,
                           const std::vector<AllowedPhiCell>& cells,
                           int scan_trials,
                           double safety)
{
    if (scan_trials <= 0) return 0.0;

    double pmax = 0.0;
    for (int i = 0; i < scan_trials; ++i) {
        Event ev{};
        if (!ProposeUniformEvent(rng, cells, ev)) continue;
        double p = component.eval ? component.eval(ev, component.ctx) : 0.0;
        if (!std::isfinite(p) || p < 0.0) p = 0.0;
        if (p > pmax) pmax = p;
    }

    if (!(pmax > 0.0) || !std::isfinite(pmax)) return 0.0;
    if (!(safety > 1.0) || !std::isfinite(safety)) safety = 5.0;
    return pmax * safety;
}

static bool GenerateEventsForComponent(const PdfComponent& component,
                                       long long n_events,
                                       std::mt19937_64& rng,
                                       const std::vector<AllowedPhiCell>& cells,
                                       const ToyGeneratorConfig& cfg,
                                       std::vector<Event>& out_events)
{
    out_events.clear();
    if (n_events <= 0) return true;
    if (cells.empty()) return false;

    std::mt19937_64 rng_scan(MixSeed(cfg.seed, static_cast<unsigned long long>(n_events + 11)));
    double safety = cfg.pmax_safety;
    if (!(safety > 1.0) || !std::isfinite(safety)) safety = 5.0;
    double update = cfg.pmax_update;
    if (!(update > 1.0) || !std::isfinite(update)) update = 1.2;

    double pmax = EstimatePMax(component, rng_scan, cells, cfg.pmax_scan_trials, safety);

    out_events.reserve(static_cast<std::size_t>(n_events));
    std::uniform_real_distribution<double> u01(0.0, 1.0);

    while (static_cast<long long>(out_events.size()) < n_events) {
        Event ev{};
        if (!ProposeUniformEvent(rng, cells, ev)) return false;

        double p = component.eval ? component.eval(ev, component.ctx) : 0.0;
        if (!std::isfinite(p) || p <= 0.0) continue;

        if (!(pmax > 0.0) || !std::isfinite(pmax)) {
            pmax = p * safety;
            continue;
        }

        if (p > pmax) {
            pmax = p * update;
        }

        const double accept_prob = p / pmax;
        if (!(accept_prob > 0.0) || !std::isfinite(accept_prob)) continue;

        if (u01(rng) < accept_prob) {
            out_events.push_back(ev);
        }
    }

    return true;
}

static FitResult FitNLLWithFixedMask(const std::vector<Event>& events,
                                     const std::vector<PdfComponent>& components,
                                     const FitConfig& cfg,
                                     const std::vector<bool>& fixed_mask,
                                     const std::vector<double>& fixed_values)
{
    const std::size_t npar = components.size();
    FitResult out = MakeFailedFitResult(npar);

    if (npar == 0) return out;
    if (cfg.start_yields.size() != npar) return out;
    if (fixed_mask.size() != npar) return out;
    if (fixed_values.size() != npar) return out;

    std::vector<std::size_t> free_indices;
    std::vector<double> start_full = cfg.start_yields;
    for (std::size_t i = 0; i < npar; ++i) {
        if (fixed_mask[i]) {
            if (!std::isfinite(fixed_values[i])) return out;
            start_full[i] = fixed_values[i];
        } else {
            free_indices.push_back(i);
            if (!std::isfinite(start_full[i])) start_full[i] = 0.0;
        }
    }

    if (free_indices.empty()) {
        out.yields_hat = start_full;
        out.nll_min = NLL(events, components, out.yields_hat);
        out.status = std::isfinite(out.nll_min) ? 0 : 2;
        out.yields_err.assign(npar, 0.0);
        return out;
    }

    std::unique_ptr<ROOT::Math::Minimizer> min(
        ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
    if (!min) return out;

    if (cfg.max_calls > 0) min->SetMaxFunctionCalls(cfg.max_calls);
    if (cfg.tol > 0.0) min->SetTolerance(cfg.tol);

    auto fcn = [&](const double* x) -> double {
        std::vector<double> yields = start_full;
        for (std::size_t j = 0; j < free_indices.size(); ++j) {
            yields[free_indices[j]] = x[j];
        }
        return NLL(events, components, yields);
    };

    ROOT::Math::Functor functor(fcn, static_cast<unsigned int>(free_indices.size()));
    min->SetFunction(functor);

    for (std::size_t j = 0; j < free_indices.size(); ++j) {
        const std::size_t i = free_indices[j];
        const double start = start_full[i];
        const double step = (std::abs(start) > 0.0) ? 0.1 * std::abs(start) : 1.0;

        // 既存 FitNLL と同じく、先頭2成分（N_sig, N_rmd）は下限 0 を付ける。
        // ただし fixed parameter はここには入らない。
        if (i == 0 || i == 1) {
            min->SetLowerLimitedVariable(static_cast<unsigned int>(j),
                                         components[i].name, start, step, 0.0);
        } else {
            min->SetVariable(static_cast<unsigned int>(j),
                             components[i].name, start, step);
        }
    }

    const bool ok = min->Minimize();
    out.status = ok ? 0 : 2;
    out.yields_hat = start_full;

    const double* xs = min->X();
    if (xs) {
        for (std::size_t j = 0; j < free_indices.size(); ++j) {
            out.yields_hat[free_indices[j]] = xs[j];
        }
    }
    out.nll_min = min->MinValue();

    out.yields_err.assign(npar, 0.0);
    const double* es = min->Errors();
    if (es) {
        for (std::size_t j = 0; j < free_indices.size(); ++j) {
            out.yields_err[free_indices[j]] = es[j];
        }
    }

    return out;
}

FitResult FitNLLFixedSignal(const std::vector<Event>& events,
                            const std::vector<PdfComponent>& components,
                            const FitConfig& cfg,
                            double N_sig_fixed)
{
    const std::size_t npar = components.size();
    FitResult out = MakeFailedFitResult(npar);
    if (npar == 0) return out;
    if (!std::isfinite(N_sig_fixed) || N_sig_fixed < 0.0) return out;

    std::vector<bool> fixed_mask(npar, false);
    std::vector<double> fixed_values(npar, 0.0);
    fixed_mask[0] = true;
    fixed_values[0] = N_sig_fixed;

    return FitNLLWithFixedMask(events, components, cfg, fixed_mask, fixed_values);
}

double EvaluateProfileLikelihoodQ(const std::vector<Event>& events,
                                  const std::vector<PdfComponent>& components,
                                  const FitConfig& free_fit_cfg,
                                  const FitConfig& prof_fit_cfg,
                                  double N_sig_fixed,
                                  FitResult& fit_free_out,
                                  FitResult& fit_prof_out)
{
    fit_free_out = FitNLL(events, components, free_fit_cfg);
    if (fit_free_out.status != 0 || !std::isfinite(fit_free_out.nll_min)) {
        fit_prof_out = MakeFailedFitResult(components.size());
        return std::numeric_limits<double>::quiet_NaN();
    }

    FitConfig prof_cfg = prof_fit_cfg;
    if (prof_cfg.start_yields.size() == components.size()) {
        prof_cfg.start_yields = fit_free_out.yields_hat;
    }
    if (prof_cfg.start_yields.size() == components.size()) {
        prof_cfg.start_yields[0] = N_sig_fixed;
    }

    fit_prof_out = FitNLLFixedSignal(events, components, prof_cfg, N_sig_fixed);
    if (fit_prof_out.status != 0 || !std::isfinite(fit_prof_out.nll_min)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    double q = 2.0 * (fit_prof_out.nll_min - fit_free_out.nll_min);
    if (!std::isfinite(q)) return std::numeric_limits<double>::quiet_NaN();
    if (q < 0.0) q = 0.0;
    return q;
}

bool GenerateToyDatasetFromModel(const std::vector<PdfComponent>& components,
                                 const std::vector<double>& mean_yields,
                                 const ToyGeneratorConfig& cfg,
                                 unsigned long long toy_index,
                                 std::vector<Event>& out_events,
                                 std::vector<double>* out_generated_yields)
{
    out_events.clear();
    if (components.empty()) return false;
    if (mean_yields.size() != components.size()) return false;

    std::vector<AllowedPhiCell> cells;
    double total_phi_area = 0.0;
    if (!BuildAllowedPhiCells(cells, total_phi_area)) return false;
    if (!(total_phi_area > 0.0)) return false;

    std::mt19937_64 rng(MixSeed(cfg.seed, toy_index + 1ULL));

    std::vector<double> generated(mean_yields.size(), 0.0);
    std::vector<Event> all_events;

    for (std::size_t k = 0; k < components.size(); ++k) {
        const double mean = mean_yields[k];
        if (!std::isfinite(mean) || mean < 0.0) return false;

        long long n_gen = 0;
        if (mean > 0.0) {
            std::poisson_distribution<long long> pois(mean);
            n_gen = pois(rng);
        }
        generated[k] = static_cast<double>(n_gen);

        std::vector<Event> comp_events;
        if (!GenerateEventsForComponent(components[k], n_gen, rng, cells, cfg, comp_events)) {
            return false;
        }
        all_events.insert(all_events.end(), comp_events.begin(), comp_events.end());
    }

    std::shuffle(all_events.begin(), all_events.end(), rng);
    out_events.swap(all_events);

    if (out_generated_yields) *out_generated_yields = generated;
    return true;
}

UpperLimitPointResult EvaluateUpperLimitPoint(const std::vector<Event>& events,
                                              const std::vector<PdfComponent>& components,
                                              const UpperLimitPointConfig& cfg)
{
    UpperLimitPointResult out;
    out.N_sig_test = cfg.N_sig_test;
    out.q_obs = std::numeric_limits<double>::quiet_NaN();
    out.p_value = 0.0;
    out.acceptance_threshold = 0.1;
    out.accepted = false;
    out.n_toys_requested = cfg.n_toys;
    out.n_toys_valid = 0;
    out.fit_free_obs = MakeFailedFitResult(components.size());
    out.fit_prof_obs = MakeFailedFitResult(components.size());

    if (cfg.cl > 0.0 && cfg.cl < 1.0 && std::isfinite(cfg.cl)) {
        out.acceptance_threshold = 1.0 - cfg.cl;
    }

    out.q_obs = EvaluateProfileLikelihoodQ(events, components,
                                           cfg.free_fit_cfg, cfg.prof_fit_cfg,
                                           cfg.N_sig_test,
                                           out.fit_free_obs, out.fit_prof_obs);
    if (!std::isfinite(out.q_obs)) {
        return out;
    }

    const std::vector<double>& mean_yields = out.fit_prof_obs.yields_hat;
    if (mean_yields.size() != components.size()) return out;

    int n_ge = 0;
    for (int itoy = 0; itoy < cfg.n_toys; ++itoy) {
        std::vector<Event> toy_events;
        std::vector<double> toy_generated_yields;
        if (!GenerateToyDatasetFromModel(components, mean_yields, cfg.toy_cfg,
                                         static_cast<unsigned long long>(itoy),
                                         toy_events, &toy_generated_yields)) {
            continue;
        }

        if (toy_events.empty()) continue;

        FitConfig free_cfg = cfg.free_fit_cfg;
        if (free_cfg.start_yields.size() == components.size()) {
            free_cfg.start_yields = toy_generated_yields;
            const double nsum =
                std::accumulate(free_cfg.start_yields.begin(), free_cfg.start_yields.end(), 0.0);
            if (!(nsum > 0.0)) {
                for (double& v : free_cfg.start_yields) v = 1.0;
            }
        }

        FitConfig prof_cfg = cfg.prof_fit_cfg;
        if (prof_cfg.start_yields.size() == components.size()) {
            prof_cfg.start_yields = toy_generated_yields;
            if (!prof_cfg.start_yields.empty()) prof_cfg.start_yields[0] = cfg.N_sig_test;
        }

        FitResult fit_free_toy;
        FitResult fit_prof_toy;
        const double q_toy = EvaluateProfileLikelihoodQ(toy_events, components,
                                                        free_cfg, prof_cfg,
                                                        cfg.N_sig_test,
                                                        fit_free_toy, fit_prof_toy);
        if (!std::isfinite(q_toy)) continue;

        ++out.n_toys_valid;
        if (q_toy >= out.q_obs) ++n_ge;
    }

    if (out.n_toys_valid > 0) {
        out.p_value = static_cast<double>(n_ge) /
                      static_cast<double>(out.n_toys_valid);
        out.accepted = (out.p_value >= out.acceptance_threshold);
    }

    return out;
}

UpperLimitScanResult EvaluateUpperLimitScan(const std::vector<Event>& events,
                                            const std::vector<PdfComponent>& components,
                                            const UpperLimitScanConfig& cfg)
{
    UpperLimitScanResult out;
    out.fit_free_obs = MakeFailedFitResult(components.size());
    out.points.clear();
    out.N_sig_90 = 0.0;
    out.BR_90 = 0.0;
    out.N_mu_eff = cfg.N_mu_eff;

    out.fit_free_obs = FitNLL(events, components, cfg.free_fit_cfg);

    for (std::size_t i = 0; i < cfg.N_sig_scan.size(); ++i) {
        UpperLimitPointConfig pcfg;
        pcfg.N_sig_test = cfg.N_sig_scan[i];
        pcfg.n_toys = cfg.n_toys_per_point;
        pcfg.cl = cfg.cl;
        pcfg.free_fit_cfg = cfg.free_fit_cfg;
        pcfg.prof_fit_cfg = cfg.prof_fit_cfg;
        pcfg.toy_cfg = cfg.toy_cfg;
        pcfg.toy_cfg.seed = MixSeed(cfg.toy_cfg.seed, static_cast<unsigned long long>(i + 101));

        UpperLimitPointResult pres = EvaluateUpperLimitPoint(events, components, pcfg);
        out.points.push_back(pres);

        if (pres.accepted && std::isfinite(pres.N_sig_test) &&
            pres.N_sig_test > out.N_sig_90) {
            out.N_sig_90 = pres.N_sig_test;
        }
    }

    if (std::isfinite(out.N_mu_eff) && out.N_mu_eff > 0.0) {
        out.BR_90 = out.N_sig_90 / out.N_mu_eff;
    }

    return out;
}
