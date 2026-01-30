// scripts/run_michel_pol_fit.cc
//
// Michel 偏極測定（別実験）: 4データセット（入替え運用）から P_mu を推定する
//
// 使い方:
//   ./build/run_michel_pol_fit  A_plus.dat  A_minus.dat  B_plus.dat  B_minus.dat
//
// 入力データ形式（Ee-only）:
//   - 1行1事象、1列: Ee [MeV]
//   - 空行・#コメント行は無視
//
// 出力:
//   - 設定値（MichelPolConfig）
//   - 推定 P_mu と 1σ
//   - chi2/ndf, 使用ビン数
//   - 読み込み統計（read/skipped）
#include <iostream>
#include <iomanip>
#include <string>

#include "p2meg/MichelPolFit.h"
#include "p2meg/MichelPolConfig.h"

static void PrintUsage(const char* prog) {
    std::cout
        << "Usage:\n"
        << "  " << prog << " A_plus.dat A_minus.dat B_plus.dat B_minus.dat\n";
}

static void PrintConfig() {
    const auto& cfg = michel_pol_config;

    std::cout << "=== MichelPolConfig ===\n";
    std::cout << "cos_theta_abs  = " << cfg.cos_theta_abs << "\n";
    std::cout << "Ee_min..max    = " << cfg.Ee_min << " .. " << cfg.Ee_max << " [MeV]\n";
    std::cout << "nbins_Ee       = " << cfg.nbins_Ee << "\n";
    std::cout << "fit_Ee_min..max= " << cfg.fit_Ee_min << " .. " << cfg.fit_Ee_max << " [MeV]\n";
    std::cout << "min_counts_each= " << cfg.min_counts_each << "\n";
    std::cout << "\n";
}

int main(int argc, char** argv) {
    if (argc != 5) {
        PrintUsage(argv[0]);
        return 2;
    }

    const std::string path_A_plus  = argv[1];
    const std::string path_A_minus = argv[2];
    const std::string path_B_plus  = argv[3];
    const std::string path_B_minus = argv[4];

    PrintConfig();

    const MichelPolFitResult res =
        FitMichelPolarizationFrom4Files(path_A_plus, path_A_minus, path_B_plus, path_B_minus);

    std::cout << "=== Fit Result ===\n";
    std::cout << "status         = " << res.status << "\n";
    std::cout << "n_read_total   = " << res.n_read_total << "\n";
    std::cout << "n_skipped_total= " << res.n_skipped_total << "\n";
    std::cout << "n_bins_used    = " << res.n_bins_used << "\n";

    if (res.status != 0) {
        std::cout << "fit failed.\n";
        return 1;
    }

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "P_mu_hat       = " << res.P_mu_hat << "\n";
    std::cout << "err_P_mu       = " << res.err_P_mu << "\n";

    if (res.ndf > 0) {
        std::cout << "chi2/ndf       = " << res.chi2 << " / " << res.ndf
                  << " = " << (res.chi2 / static_cast<double>(res.ndf)) << "\n";
    } else {
        std::cout << "chi2           = " << res.chi2 << "\n";
        std::cout << "ndf            = " << res.ndf << "\n";
    }

    return 0;
}
