// scripts/evaluate_rmd_aw_sm.cc
//
// 目的:
//   現在の detector model / AW 定義に整合する形で、
//   SM の RMD が AW に何事象入るかを MC 積分で見積もる。
//
// 出力:
//   - kappa_RMD_AW : 1 stopped muon あたりの AW 内再構成 RMD 期待事象数
//   - BR_eff_RMD   : signal 正規化 N_mu^eff で割った「実効分岐比」
//   - N_RMD_SM_AW  : 入力した N_mu^eff に対する AW 内期待 RMD 事象数
//   - 現在の尤度 fit の N_RMD_hat との比較
//
// 注意:
//   - 角度は現在の解析実装と同じく phi グリッドと許可マスクで扱う。
//   - 連続角の Jacobian は掛けず、phi のビン幅を測度として使う。
//     これは current detector model に整合した「実効的」な RMD 率であり、
//     連続 4pi 立体角での厳密な物理 branching ratio そのものではない。
//   - 真値窓は MakeRMDGridPdf と同様に energy response の 0.1 peak 点まで広げる。

#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "TRandom3.h"

#include "p2meg/AnalysisWindow.h"
#include "p2meg/AnalysisWindowUtils.h"
#include "p2meg/Constants.h"
#include "p2meg/DetectorResolution.h"
#include "p2meg/RMDSpectrum.h"

struct AllowedPair {
    double phi_e;
    double phi_g;
    double phi_weight; // dphi_e * dphi_g に対応する current model の測度
};

struct SummaryValues {
    double N_mu_eff = 0.0;
    double sigma_N_mu_eff = 0.0;
    double kappa_sig = 0.0;
    double sigma_kappa_sig = 0.0;
};

struct FitValues {
    double N_rmd_hat = 0.0;
    double err_rmd = 0.0;
};

static bool ParseKeyValueFile(const char* filepath,
                              const char* key,
                              double& value_out)
{
    value_out = 0.0;
    std::ifstream fin(filepath);
    if (!fin) return false;

    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') continue;

        std::istringstream iss(line);
        std::string k;
        double v = 0.0;
        if (!(iss >> k >> v)) continue;
        if (k == key) {
            value_out = v;
            return std::isfinite(value_out);
        }
    }
    return false;
}

static bool LoadSummaryValues(const char* filepath, SummaryValues& out)
{
    return ParseKeyValueFile(filepath, "N_mu_eff", out.N_mu_eff) &&
           ParseKeyValueFile(filepath, "sigma_N_mu_eff", out.sigma_N_mu_eff) &&
           ParseKeyValueFile(filepath, "kappa_sig", out.kappa_sig) &&
           ParseKeyValueFile(filepath, "sigma_kappa_sig", out.sigma_kappa_sig);
}

static bool LoadFitValues(const char* filepath, FitValues& out)
{
    std::ifstream fin(filepath);
    if (!fin) return false;

    std::string line;
    bool ok_rmd = false;
    bool ok_err = false;
    while (std::getline(fin, line)) {
        std::istringstream iss(line);
        std::string key;
        if (!(iss >> key)) continue;

        if (key == "N_rmd_hat") {
            char eq = '\0';
            double v = 0.0;
            if (iss >> eq >> v) {
                out.N_rmd_hat = v;
                ok_rmd = std::isfinite(v);
            }
        } else if (key == "err_rmd") {
            char eq = '\0';
            double v = 0.0;
            if (iss >> eq >> v) {
                out.err_rmd = v;
                ok_err = std::isfinite(v);
            }
        }
    }
    return ok_rmd && ok_err;
}

static double TimeAcceptanceInAW()
{
    if (!(detres.sigma_t > 0.0)) return 0.0;
    const double a = (analysis_window.t_min - detres.t_mean) /
                     (std::sqrt(2.0) * detres.sigma_t);
    const double b = (analysis_window.t_max - detres.t_mean) /
                     (std::sqrt(2.0) * detres.sigma_t);
    return 0.5 * (std::erf(b) - std::erf(a));
}

static AnalysisWindow4D ExpandTruthWindowByResponseLocal()
{
    AnalysisWindow4D out = analysis_window;

    const double Ee_low  = energy_response_offset_low_e(analysis_window.Ee_min);
    const double Ee_high = energy_response_offset_high_e(analysis_window.Ee_max);
    const double Eg_low  = energy_response_offset_low_g(analysis_window.Eg_min);
    const double Eg_high = energy_response_offset_high_g(analysis_window.Eg_max);

    out.Ee_min = analysis_window.Ee_min - Ee_high;
    out.Ee_max = analysis_window.Ee_max + Ee_low;
    out.Eg_min = analysis_window.Eg_min - Eg_high;
    out.Eg_max = analysis_window.Eg_max + Eg_low;

    const double mmu = kMassesPDG.m_mu;
    const double me  = kMassesPDG.m_e;
    const double r = (me * me) / (mmu * mmu);
    const double Ee_phys_max = 0.5 * mmu * (1.0 + r);
    const double Eg_phys_max = 0.5 * mmu * (1.0 - r);

    if (out.Ee_min < 0.0) out.Ee_min = 0.0;
    if (out.Eg_min < 1e-6) out.Eg_min = 1e-6;
    if (out.Ee_max > Ee_phys_max) out.Ee_max = Ee_phys_max;
    if (out.Eg_max > Eg_phys_max) out.Eg_max = Eg_phys_max;

    return out;
}

static void BuildAllowedPairs(std::vector<AllowedPair>& out_pairs,
                              double& total_phi_weight)
{
    out_pairs.clear();
    total_phi_weight = 0.0;

    const int N_phi_e = detres.N_phi_e;
    const int N_phi_g = detres.N_phi_g;

    for (int ie = 0; ie <= N_phi_e; ++ie) {
        for (int ig = 0; ig <= N_phi_g; ++ig) {
            if (!Detector_IsAllowedPhiPairIndex(ie, ig, detres)) continue;

            const double phi_e =
                Detector_PhiGridPoint(ie, detres.phi_e_min, detres.phi_e_max, N_phi_e);
            const double phi_g =
                Detector_PhiGridPoint(ig, detres.phi_g_min, detres.phi_g_max, N_phi_g);
            const double theta = std::fabs(phi_e - phi_g);
            if (!AnalysisWindow_In3D(analysis_window, analysis_window.Ee_min,
                                     analysis_window.Eg_min, theta)) {
                continue;
            }

            AllowedPair p{};
            p.phi_e = phi_e;
            p.phi_g = phi_g;
            p.phi_weight =
                Detector_PhiBinWidth(ie, detres.phi_e_min, detres.phi_e_max, N_phi_e) *
                Detector_PhiBinWidth(ig, detres.phi_g_min, detres.phi_g_max, N_phi_g);
            if (!(p.phi_weight > 0.0) || !std::isfinite(p.phi_weight)) continue;

            out_pairs.push_back(p);
            total_phi_weight += p.phi_weight;
        }
    }
}

int main(int argc, char** argv)
{
    const long long n_mc = (argc >= 2) ? std::atoll(argv[1]) : 800000;
    const unsigned long seed = (argc >= 3) ? std::strtoul(argv[2], nullptr, 10) : 20260312UL;

    if (!(n_mc > 0)) {
        std::cerr << "[evaluate_rmd_aw_sm] n_mc must be positive\n";
        return 1;
    }

    SummaryValues summary;
    if (!LoadSummaryValues("doc/finalanalysis/mu_stop_efficiency_cancellation_summary.txt", summary)) {
        std::cerr << "[evaluate_rmd_aw_sm] failed to load mu-stop summary\n";
        return 2;
    }

    FitValues fit;
    if (!LoadFitValues("doc/finalanalysis/run_nll_fit_step4f_run1to20_eg_doubleonly_allpatterns_500ns_accconstrained_current.txt", fit)) {
        std::cerr << "[evaluate_rmd_aw_sm] failed to load fit result\n";
        return 2;
    }

    const AnalysisWindow4D truth_win = ExpandTruthWindowByResponseLocal();
    const double dEe = truth_win.Ee_max - truth_win.Ee_min;
    const double dEg = truth_win.Eg_max - truth_win.Eg_min;
    if (!(dEe > 0.0) || !(dEg > 0.0)) {
        std::cerr << "[evaluate_rmd_aw_sm] invalid truth window\n";
        return 3;
    }

    std::vector<AllowedPair> pairs;
    double total_phi_weight = 0.0;
    BuildAllowedPairs(pairs, total_phi_weight);
    if (pairs.empty() || !(total_phi_weight > 0.0)) {
        std::cerr << "[evaluate_rmd_aw_sm] no allowed theta-window phi pairs\n";
        return 3;
    }

    TRandom3 rng(seed);

    const double time_acc = TimeAcceptanceInAW();
    const double volume_E = dEe * dEg;

    double sum_x = 0.0;
    double sum_x2 = 0.0;

    for (long long i = 0; i < n_mc; ++i) {
        double u = rng.Uniform(0.0, total_phi_weight);
        std::size_t ip = 0;
        double accw = 0.0;
        for (; ip < pairs.size(); ++ip) {
            accw += pairs[ip].phi_weight;
            if (u <= accw) break;
        }
        if (ip >= pairs.size()) ip = pairs.size() - 1;

        const AllowedPair& p = pairs[ip];
        const double Ee_true = rng.Uniform(truth_win.Ee_min, truth_win.Ee_max);
        const double Eg_true = rng.Uniform(truth_win.Eg_min, truth_win.Eg_max);

        const double cosThetaE  = std::cos(p.phi_e);
        const double cosThetaG  = std::cos(p.phi_g);
        const double cosThetaEG = std::cos(p.phi_e - p.phi_g);

        const double w = RMD_d6B_dEe_dEg_dOmegae_dOmegag(
            Ee_true, Eg_true,
            cosThetaEG, cosThetaE, cosThetaG,
            detres.P_mu, 1e-6);
        if (!(w > 0.0) || !std::isfinite(w)) continue;

        const double Ee_obs = smear_energy_trandom3_e(rng, Ee_true);
        const double Eg_obs = smear_energy_trandom3_g(rng, Eg_true);
        const double theta = std::fabs(p.phi_e - p.phi_g);

        double x = 0.0;
        if (AnalysisWindow_In4D(analysis_window, Ee_obs, Eg_obs, 0.0, theta)) {
            x = volume_E * total_phi_weight * time_acc * w;
        }

        sum_x += x;
        sum_x2 += x * x;
    }

    const double n = static_cast<double>(n_mc);
    const double kappa_rmd_aw = sum_x / n;
    double var_mean = 0.0;
    if (n_mc > 1) {
        const double ex2 = sum_x2 / n;
        const double ex = kappa_rmd_aw;
        double var = ex2 - ex * ex;
        if (var < 0.0) var = 0.0;
        var_mean = var / n;
    }
    const double sigma_kappa_rmd_aw = std::sqrt(var_mean);

    const double br_eff_sm =
        (summary.kappa_sig > 0.0) ? (kappa_rmd_aw / summary.kappa_sig) : 0.0;
    const double sigma_br_eff_sm =
        (br_eff_sm > 0.0)
            ? br_eff_sm * std::sqrt(
                  std::pow(sigma_kappa_rmd_aw / kappa_rmd_aw, 2) +
                  std::pow(summary.sigma_kappa_sig / summary.kappa_sig, 2))
            : 0.0;

    const double n_rmd_sm_aw = summary.N_mu_eff * br_eff_sm;
    const double sigma_n_rmd_sm_aw =
        (n_rmd_sm_aw > 0.0)
            ? n_rmd_sm_aw * std::sqrt(
                  std::pow(sigma_br_eff_sm / br_eff_sm, 2) +
                  std::pow(summary.sigma_N_mu_eff / summary.N_mu_eff, 2))
            : 0.0;

    const double br_eff_fit =
        (summary.N_mu_eff > 0.0) ? (fit.N_rmd_hat / summary.N_mu_eff) : 0.0;
    const double sigma_br_eff_fit =
        (summary.N_mu_eff > 0.0 && fit.err_rmd >= 0.0)
            ? br_eff_fit * std::sqrt(
                  std::pow(fit.err_rmd / fit.N_rmd_hat, 2) +
                  std::pow(summary.sigma_N_mu_eff / summary.N_mu_eff, 2))
            : 0.0;

    std::ofstream ofs("doc/finalanalysis/rmd_aw_sm_summary.txt");
    ofs << std::setprecision(10);
    ofs << "n_mc " << n_mc << "\n";
    ofs << "seed " << seed << "\n";
    ofs << "truth_Ee_min " << truth_win.Ee_min << "\n";
    ofs << "truth_Ee_max " << truth_win.Ee_max << "\n";
    ofs << "truth_Eg_min " << truth_win.Eg_min << "\n";
    ofs << "truth_Eg_max " << truth_win.Eg_max << "\n";
    ofs << "time_acceptance " << time_acc << "\n";
    ofs << "total_phi_weight " << total_phi_weight << "\n";
    ofs << "n_allowed_pairs " << pairs.size() << "\n";
    ofs << "kappa_sig " << summary.kappa_sig << "\n";
    ofs << "sigma_kappa_sig " << summary.sigma_kappa_sig << "\n";
    ofs << "N_mu_eff " << summary.N_mu_eff << "\n";
    ofs << "sigma_N_mu_eff " << summary.sigma_N_mu_eff << "\n";
    ofs << "kappa_RMD_AW " << kappa_rmd_aw << "\n";
    ofs << "sigma_kappa_RMD_AW " << sigma_kappa_rmd_aw << "\n";
    ofs << "BR_eff_RMD_SM " << br_eff_sm << "\n";
    ofs << "sigma_BR_eff_RMD_SM " << sigma_br_eff_sm << "\n";
    ofs << "N_RMD_SM_AW " << n_rmd_sm_aw << "\n";
    ofs << "sigma_N_RMD_SM_AW " << sigma_n_rmd_sm_aw << "\n";
    ofs << "N_RMD_hat_fit " << fit.N_rmd_hat << "\n";
    ofs << "sigma_N_RMD_hat_fit " << fit.err_rmd << "\n";
    ofs << "BR_eff_RMD_fit " << br_eff_fit << "\n";
    ofs << "sigma_BR_eff_RMD_fit " << sigma_br_eff_fit << "\n";

    std::cout << std::setprecision(10);
    std::cout << "[evaluate_rmd_aw_sm] n_mc = " << n_mc << "\n";
    std::cout << "[evaluate_rmd_aw_sm] truth window: Ee=[" << truth_win.Ee_min
              << "," << truth_win.Ee_max << "], Eg=[" << truth_win.Eg_min
              << "," << truth_win.Eg_max << "]\n";
    std::cout << "[evaluate_rmd_aw_sm] time acceptance in AW = " << time_acc << "\n";
    std::cout << "[evaluate_rmd_aw_sm] allowed phi pairs = " << pairs.size()
              << ", total_phi_weight = " << total_phi_weight << "\n";
    std::cout << "[evaluate_rmd_aw_sm] kappa_RMD_AW = " << kappa_rmd_aw
              << " +/- " << sigma_kappa_rmd_aw << "\n";
    std::cout << "[evaluate_rmd_aw_sm] BR_eff_RMD_SM = " << br_eff_sm
              << " +/- " << sigma_br_eff_sm << "\n";
    std::cout << "[evaluate_rmd_aw_sm] N_RMD_SM_AW = " << n_rmd_sm_aw
              << " +/- " << sigma_n_rmd_sm_aw << "\n";
    std::cout << "[evaluate_rmd_aw_sm] N_RMD_hat_fit = " << fit.N_rmd_hat
              << " +/- " << fit.err_rmd << "\n";
    std::cout << "[evaluate_rmd_aw_sm] BR_eff_RMD_fit = " << br_eff_fit
              << " +/- " << sigma_br_eff_fit << "\n";
    std::cout << "[evaluate_rmd_aw_sm] wrote doc/finalanalysis/rmd_aw_sm_summary.txt\n";

    return 0;
}
