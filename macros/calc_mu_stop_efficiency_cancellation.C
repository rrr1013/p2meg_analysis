// macros/calc_mu_stop_efficiency_cancellation.C
//
// 目的:
//   - Michel 高エネルギー窓を使った正規化で、有効停止ミューオン数を
//     少数の入力だけからまとめて計算する。
//
// このマクロで採用する式:
//   N_mu_eff
//     = N_HE^data
//       * (N_sig^MC / N_mu^stop,MC(sig))
//       / (N_HE^Michel,MC / N_mu^stop,MC(Michel))
//       * eta_gamma
//     = N_HE^data
//       * N_sig^MC
//       * N_mu^stop,MC(Michel)
//       / N_mu^stop,MC(sig)
//       / N_HE^Michel,MC
//       * eta_gamma
//
// 注意:
//   - 依頼に合わせて、入力は
//       N_HE^data,
//       N_mu^stop,MC(Michel),
//       N_HE^Michel,MC,
//       N_mu^stop,MC(sig),
//       N_sig^MC,
//       eta_gamma
//     をコード上部にハードコードする。
//   - 誤差は各カウントの統計誤差を sqrt(N) とする Gaussian 近似で評価する。
//   - eta_gamma の不確かさも入れたい場合は sigma_eta_gamma を設定する。
//
// 実行例:
//   root -l -q 'macros/calc_mu_stop_efficiency_cancellation.C'

R__ADD_INCLUDE_PATH(./include)

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "TString.h"
#include "TSystem.h"

// =========================
// ここを直接書き換えて使う
// =========================
static constexpr double kN_HE_data = 6542;                // N_HE^data
static constexpr double kN_mu_stop_MC_michel = 5e5 * 3;  // N_mu^stop,MC(Michel)
static constexpr double kN_HE_Michel_MC = 13678;         // N_HE^Michel,MC
static constexpr double kN_mu_stop_MC_sig = 5e7 * 3;     // N_mu^stop,MC(sig)
static constexpr double kN_sig_MC = 357173;               // N_sig^MC
static constexpr double kEta_gamma = 1;               // eta_gamma
static constexpr double kSigma_eta_gamma = 0.0;          // sigma(eta_gamma)

static bool IsValidCount(double x)
{
    return std::isfinite(x) && x > 0.0;
}

static double PoissonSigma(double n)
{
    if (!std::isfinite(n) || n < 0.0) return 0.0;
    return std::sqrt(n);
}

static double SafeRelativeSigma(double n)
{
    if (!IsValidCount(n)) return 0.0;
    return PoissonSigma(n) / n;
}

void calc_mu_stop_efficiency_cancellation()
{
    const double N_HE_data = kN_HE_data;
    const double N_mu_stop_MC_michel = kN_mu_stop_MC_michel;
    const double N_HE_Michel_MC = kN_HE_Michel_MC;
    const double N_mu_stop_MC_sig = kN_mu_stop_MC_sig;
    const double N_sig_MC = kN_sig_MC;
    const double eta_gamma = kEta_gamma;
    const double sigma_eta_gamma = kSigma_eta_gamma;

    if (!IsValidCount(N_HE_data) ||
        !IsValidCount(N_mu_stop_MC_michel) ||
        !IsValidCount(N_HE_Michel_MC) ||
        !IsValidCount(N_mu_stop_MC_sig) ||
        !IsValidCount(N_sig_MC) ||
        !std::isfinite(eta_gamma) ||
        !(eta_gamma > 0.0) ||
        !(eta_gamma <= 1.0) ||
        !std::isfinite(sigma_eta_gamma) ||
        sigma_eta_gamma < 0.0) {
        std::cerr << "[calc_mu_stop_efficiency_cancellation] "
                  << "入力値が不正です. count は正の有限値, eta_gamma は 0 < eta <= 1,"
                  << " sigma_eta_gamma は 0 以上にしてください." << std::endl;
        return;
    }

    // Michel 側高エネルギー窓の 1 stop muon あたり期待カウント
    const double kappa_michel = N_HE_Michel_MC / N_mu_stop_MC_michel;

    // signal 側の 1 stop muon あたり期待 event 数
    const double kappa_sig = N_sig_MC / N_mu_stop_MC_sig;

    const double N_mu_eff = N_HE_data * (kappa_sig / kappa_michel) * eta_gamma;

    // 各入力を独立とみなして相対誤差を伝播する。
    const double rel_sigma_sq =
        std::pow(SafeRelativeSigma(N_HE_data), 2) +
        std::pow(SafeRelativeSigma(N_sig_MC), 2) +
        std::pow(SafeRelativeSigma(N_mu_stop_MC_michel), 2) +
        std::pow(SafeRelativeSigma(N_mu_stop_MC_sig), 2) +
        std::pow(SafeRelativeSigma(N_HE_Michel_MC), 2) +
        ((eta_gamma > 0.0) ? std::pow(sigma_eta_gamma / eta_gamma, 2) : 0.0);
    const double sigma_N_mu_eff =
        (N_mu_eff > 0.0 && std::isfinite(rel_sigma_sq))
            ? N_mu_eff * std::sqrt(rel_sigma_sq)
            : 0.0;

    const double sigma_kappa_michel =
        (kappa_michel > 0.0)
            ? kappa_michel * std::sqrt(
                  std::pow(SafeRelativeSigma(N_HE_Michel_MC), 2) +
                  std::pow(SafeRelativeSigma(N_mu_stop_MC_michel), 2))
            : 0.0;

    const double sigma_kappa_sig =
        (kappa_sig > 0.0)
            ? kappa_sig * std::sqrt(
                  std::pow(SafeRelativeSigma(N_sig_MC), 2) +
                  std::pow(SafeRelativeSigma(N_mu_stop_MC_sig), 2))
            : 0.0;

    gSystem->mkdir("doc/finalanalysis", /*recursive=*/true);
    const TString outfile =
        "doc/finalanalysis/mu_stop_efficiency_cancellation_summary.txt";
    std::ofstream ofs(outfile.Data());
    if (!ofs) {
        std::cerr << "[calc_mu_stop_efficiency_cancellation] "
                  << "出力ファイルを開けません: " << outfile << std::endl;
        return;
    }

    ofs << std::setprecision(12);
    ofs << "# Michel HE normalization for effective stopped muons\n";
    ofs << "# Formula implemented in this macro:\n";
    ofs << "#   N_mu_eff = N_HE_data * (N_sig_MC / N_mu_stop_MC_sig)\n";
    ofs << "#              / (N_HE_Michel_MC / N_mu_stop_MC_michel)\n";
    ofs << "#              * eta_gamma\n";
    ofs << "#            = N_HE_data * N_sig_MC * N_mu_stop_MC_michel\n";
    ofs << "#              / N_mu_stop_MC_sig / N_HE_Michel_MC * eta_gamma\n";
    ofs << "# Statistical uncertainty: sigma(N)=sqrt(N), Gaussian propagation.\n";
    ofs << "\n";
    ofs << "N_HE_data            " << N_HE_data << "\n";
    ofs << "N_mu_stop_MC_michel  " << N_mu_stop_MC_michel << "\n";
    ofs << "N_HE_Michel_MC       " << N_HE_Michel_MC << "\n";
    ofs << "N_mu_stop_MC_sig     " << N_mu_stop_MC_sig << "\n";
    ofs << "N_sig_MC             " << N_sig_MC << "\n";
    ofs << "eta_gamma            " << eta_gamma << "\n";
    ofs << "sigma_eta_gamma      " << sigma_eta_gamma << "\n";
    ofs << "\n";
    ofs << "kappa_michel         " << kappa_michel << "\n";
    ofs << "sigma_kappa_michel   " << sigma_kappa_michel << "\n";
    ofs << "kappa_sig            " << kappa_sig << "\n";
    ofs << "sigma_kappa_sig      " << sigma_kappa_sig << "\n";
    ofs << "\n";
    ofs << "N_mu_eff             " << N_mu_eff << "\n";
    ofs << "sigma_N_mu_eff       " << sigma_N_mu_eff << "\n";
    ofs << "relative_sigma       "
        << ((N_mu_eff > 0.0) ? sigma_N_mu_eff / N_mu_eff : 0.0) << "\n";

    ofs.close();

    std::cout << std::setprecision(12);
    std::cout << "=== Michel HE normalization summary ===" << std::endl;
    std::cout << "N_HE^data              = " << N_HE_data << std::endl;
    std::cout << "N_mu^stop,MC(Michel)   = " << N_mu_stop_MC_michel << std::endl;
    std::cout << "N_HE^Michel,MC         = " << N_HE_Michel_MC << std::endl;
    std::cout << "N_mu^stop,MC(sig)      = " << N_mu_stop_MC_sig << std::endl;
    std::cout << "N_sig^MC               = " << N_sig_MC << std::endl;
    std::cout << "eta_gamma              = " << eta_gamma
              << " +/- " << sigma_eta_gamma << std::endl;
    std::cout << std::endl;
    std::cout << "kappa^Michel           = "
              << kappa_michel << " +/- " << sigma_kappa_michel << std::endl;
    std::cout << "kappa^sig              = "
              << kappa_sig << " +/- " << sigma_kappa_sig << std::endl;
    std::cout << std::endl;
    std::cout << "N_mu^eff               = "
              << N_mu_eff << " +/- " << sigma_N_mu_eff << std::endl;
    std::cout << "relative uncertainty   = "
              << ((N_mu_eff > 0.0) ? sigma_N_mu_eff / N_mu_eff : 0.0)
              << std::endl;
    std::cout << "summary file           = " << outfile << std::endl;
}
