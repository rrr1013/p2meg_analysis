#include "TString.h"
#include "TSystem.h"

// ============================================================
// 感度計算の実行ラッパ
//
// 使い方:
//   root -l -q 'macros/run_sensitivity_br_toymc.C(
//     "data/finaldata/step4f_run1to20_eg_doubleonly_allpatterns_500ns.txt",
//     1708.77883387, 25.8915493204,
//     0.0, 0.0012, 0.00005,
//     400, 1000, 246813579ULL, 8, 200,
//     "doc/finalanalysis/sensitivity_br_scan_blind50ns.txt")'
//
// 注意:
//   - 実計算本体は build/run_sensitivity_br_toymc
//   - outfile を与えると標準出力をそのまま保存する
// ============================================================

void run_sensitivity_br_toymc(
    const char* datafile =
        "data/finaldata/step4f_run1to20_eg_doubleonly_allpatterns_500ns.txt",
    double N_mu_eff_nom = 1708.77883387,
    double N_mu_eff_sigma = 25.8915493204,
    double br_min = 0.0,
    double br_max = 0.0012,
    double br_step = 0.00005,
    int n_calib_toys = 400,
    int n_sens_toys = 1000,
    unsigned long long seed = 246813579ULL,
    int n_validation_toys = 8,
    int n_validation_inner = 200,
    const char* outfile = "")
{
    TString cmd;
    cmd.Form("./build/run_sensitivity_br_toymc "
             "\"%s\" %.17g %.17g %.17g %.17g %.17g %d %d %llu %d %d",
             datafile,
             N_mu_eff_nom, N_mu_eff_sigma,
             br_min, br_max, br_step,
             n_calib_toys, n_sens_toys,
             seed, n_validation_toys, n_validation_inner);

    if (outfile && outfile[0] != '\0') {
        cmd += " > \"";
        cmd += outfile;
        cmd += "\"";
    }

    gSystem->Exec(cmd.Data());
}
