// scripts/make_acc_mockdata.cc
//
// ============================================================
// p2MEG: ACC（accidental）擬似データ（.dat）を生成する
//
// 出力仕様（make_pdfmix_mockdata.cc に準拠）:
//  - 引数: n_acc のみ
//  - 出力: data/mockdata/acc_<n_acc>.dat
//  - 出力列(5列): Ee  Eg  t  phi_detector_e  phi_detector_g
//
// 生成方針（簡略ACCモデル）:
//  - e 側（Ee, phi_e）: Michel の 2D shape（Ee, cosθe）から棄却法で生成
//      * 観測角 phi_detector_e は DetectorResolution の格子点で提案
//      * cosθe は cos(phi_detector_e) と同一視（同一平面近似）
//  - γ 側（Eg, phi_g）: RMD の完全式（偏極込み）を用い、(Ee_dummy,Eg,phi_e_dummy,phi_g) を棄却法で生成し、
//                      受理されたら (Eg, phi_g) のみ採用（Ee_dummy,phi_e_dummy は捨てる）
//      * cosθe = cos(phi_e_dummy), cosθg = cos(phi_g), cosθeg = cos(phi_e_dummy-phi_g)
//  - エネルギースメア:
//      * DetectorResolution.h の energy_response_shape_e/g(E_res,E_true) に従う
//      * 観測値が解析窓外なら、その側だけ再生成（窓内条件付け）
//  - t: 解析窓内で一様分布（analysis_window.t_min..t_max）
//
// 進捗表示:
//  - acc について進捗バーを stderr に表示
//
// 自動設定:
//  - 解析窓: analysis_window
//  - 分解能: detres（N_phi_e/g, P_mu など）
//  - 質量:   kMassesPDG（内部で使う場合）
//
// 実行例:
//   ./build/make_acc_mockdata 5000
// ============================================================

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cerrno>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "TError.h"
#include "TRandom3.h"

#include "../include/p2meg/Event.h"
#include "../include/p2meg/AnalysisWindow.h"
#include "../include/p2meg/DetectorResolution.h"
#include "../include/p2meg/Constants.h"
#include "../include/p2meg/MathUtils.h"
#include "../include/p2meg/MichelSpectrum.h"
#include "../include/p2meg/RMDSpectrum.h"

// ---- 自動設定（必要ならここだけ変更）----
static const char* kOutDir = "data/mockdata";

// pmax 初期推定のプリスキャン回数（0でも動くが、あると少し速い）
static constexpr int kScanTrialsMichel = 200000;
static constexpr int kScanTrialsRmd6   = 200000;

// pmax 安全係数（初期化/更新の“余裕”）
static constexpr double kSafetyFactor = 5.0;

// pmax が不足していたときの更新係数（局所的に引き上げる）
static constexpr double kUpdateFactor = 1.2;

// 進捗バー幅
static constexpr int kBarWidth = 40;

// 進捗更新間隔 [ms]
static constexpr int kProgressIntervalMs = 200;

// ------------------------------------------------------------
// 文字列→long long
// ------------------------------------------------------------
static bool ParseLL(const char* s, long long& out)
{
    if (!s) return false;
    char* end = nullptr;
    errno = 0;
    long long v = std::strtoll(s, &end, 10);
    if (errno != 0) return false;
    if (end == s || *end != '\0') return false;
    if (v < 0) return false;
    out = v;
    return true;
}

// ------------------------------------------------------------
// 出力用: 検出器角を DetectorResolution の格子点に丸める
// 物理カットではなく出力整形のための丸め
// ------------------------------------------------------------
static double SnapPhiToGrid(double phi, double phi_min, double phi_max, int N_phi)
{
    return Detector_PhiSnapToGrid(phi, phi_min, phi_max, N_phi);
}

// ------------------------------------------------------------
// 進捗バー表示（stderr）
// ------------------------------------------------------------
static void PrintProgressBar(const char* label,
                             long long done, long long total,
                             long long tries_e, long long tries_g,
                             double pmax_e, double pmax_g,
                             double acc_rate)
{
    if (total <= 0) return;

    const double frac = (done <= 0) ? 0.0 : (done >= total ? 1.0 : (double(done) / double(total)));
    const int filled = static_cast<int>(std::round(frac * kBarWidth));

    std::ostringstream oss;
    oss << "\r" << label << " [";
    for (int i = 0; i < kBarWidth; ++i) oss << (i < filled ? '#' : '-');
    oss << "] " << std::fixed << std::setprecision(1) << (frac * 100.0)
        << "% (" << done << "/" << total << ")"
        << " tries_e=" << tries_e
        << " tries_g=" << tries_g
        << " acc=" << std::setprecision(3) << (acc_rate * 100.0) << "%"
        << " pmax_e=" << std::setprecision(6) << pmax_e
        << " pmax_g=" << std::setprecision(6) << pmax_g;

    std::cerr << oss.str() << std::flush;

    if (done >= total) {
        std::cerr << "\n";
    }
}

// ------------------------------------------------------------
// 一様提案（Michel: Ee, phi_e(grid)）
// ------------------------------------------------------------
struct MichelCand {
    double Ee_true;
    double phi_e;   // grid point [0,pi]
    double w;       // Michel shape (non-normalized)
};

static MichelCand ProposeMichel2D(TRandom3& rng, double Pmu)
{
    const int N_phi = Math_GetNPhiE(detres);
    const double step = Detector_PhiStep(detres.phi_e_min, detres.phi_e_max, N_phi);

    const double Ee = rng.Uniform(analysis_window.Ee_min, analysis_window.Ee_max);
    const int ui = static_cast<int>(rng.Integer(static_cast<ULong_t>(N_phi + 1)));
    const double phi = detres.phi_e_min + step * static_cast<double>(ui);
    const double costh = std::cos(phi);

    MichelCand c;
    c.Ee_true = Ee;
    c.phi_e   = phi;
    c.w       = Michel_d2Shape_dE_dCosTheta(Ee, costh, Pmu, kMichelSM, kMassesPDG);
    return c;
}

static double EstimatePMaxMichel2D(TRandom3& rng, int trials, double Pmu)
{
    double pmax = 0.0;
    for (int i = 0; i < trials; ++i) {
        const MichelCand c = ProposeMichel2D(rng, Pmu);
        if (std::isfinite(c.w) && c.w > pmax) pmax = c.w;
    }
    if (!(pmax > 0.0)) return 0.0;
    return pmax * kSafetyFactor;
}

static bool SampleMichel2D_AR(TRandom3& rng, double Pmu, double& pmax_inout,
                             double& Ee_true_out, double& phi_e_out, long long& tries_inout)
{
    while (true) {
        ++tries_inout;
        const MichelCand c = ProposeMichel2D(rng, Pmu);

        const double w = c.w;
        if (!std::isfinite(w) || w <= 0.0) continue;

        if (!(pmax_inout > 0.0) || !std::isfinite(pmax_inout)) {
            pmax_inout = w * kSafetyFactor;
            continue; // 初期化した回は判定しない（安全側）
        }
        if (w > pmax_inout) {
            pmax_inout = w * kUpdateFactor;
        }

        const double u = rng.Rndm();
        if (u < (w / pmax_inout)) {
            Ee_true_out = c.Ee_true;
            phi_e_out   = c.phi_e;
            return true;
        }
    }
}

// ------------------------------------------------------------
// 一様提案（RMD6: Ee_dummy, Eg_true, phi_e_dummy(grid), phi_g(grid)）
// 受理後に (Eg_true, phi_g) だけ使う
// ------------------------------------------------------------
struct Rmd6Cand {
    double Ee_dummy;
    double Eg_true;
    double phi_e_dummy;
    double phi_g;
    double w; // RMD_d6B...
};

static Rmd6Cand ProposeRmd6_4D(TRandom3& rng, double Pmu)
{
    const int N_phi_e = Math_GetNPhiE(detres);
    const int N_phi_g = Math_GetNPhiG(detres);
    const double step_e = Detector_PhiStep(detres.phi_e_min, detres.phi_e_max, N_phi_e);
    const double step_g = Detector_PhiStep(detres.phi_g_min, detres.phi_g_max, N_phi_g);

    const double Ee = rng.Uniform(analysis_window.Ee_min, analysis_window.Ee_max);
    const double Eg = rng.Uniform(analysis_window.Eg_min, analysis_window.Eg_max);
    const int ui_e = static_cast<int>(rng.Integer(static_cast<ULong_t>(N_phi_e + 1)));
    const int ui_g = static_cast<int>(rng.Integer(static_cast<ULong_t>(N_phi_g + 1)));
    const double phie = detres.phi_e_min + step_e * static_cast<double>(ui_e);
    const double phig = detres.phi_g_min + step_g * static_cast<double>(ui_g);

    const double cosE  = std::cos(phie);
    const double cosG  = std::cos(phig);
    const double cosEG = std::cos(phie - phig);

    Rmd6Cand c;
    c.Ee_dummy    = Ee;
    c.Eg_true     = Eg;
    c.phi_e_dummy = phie;
    c.phi_g       = phig;
    c.w = RMD_d6B_dEe_dEg_dOmegae_dOmegag(Ee, Eg, cosEG, cosE, cosG, Pmu);
    return c;
}

static double EstimatePMaxRmd6_4D(TRandom3& rng, int trials, double Pmu)
{
    double pmax = 0.0;
    for (int i = 0; i < trials; ++i) {
        const Rmd6Cand c = ProposeRmd6_4D(rng, Pmu);
        if (std::isfinite(c.w) && c.w > pmax) pmax = c.w;
    }
    if (!(pmax > 0.0)) return 0.0;
    return pmax * kSafetyFactor;
}

static bool SampleGammaFromRmd6_AR(TRandom3& rng, double Pmu, double& pmax_inout,
                                  double& Eg_true_out, double& phi_g_out, long long& tries_inout)
{
    while (true) {
        ++tries_inout;
        const Rmd6Cand c = ProposeRmd6_4D(rng, Pmu);

        const double w = c.w;
        if (!std::isfinite(w) || w <= 0.0) continue;

        if (!(pmax_inout > 0.0) || !std::isfinite(pmax_inout)) {
            pmax_inout = w * kSafetyFactor;
            continue; // 初期化した回は判定しない（安全側）
        }
        if (w > pmax_inout) {
            pmax_inout = w * kUpdateFactor;
        }

        const double u = rng.Rndm();
        if (u < (w / pmax_inout)) {
            Eg_true_out = c.Eg_true;
            phi_g_out   = c.phi_g;
            return true;
        }
    }
}

// ------------------------------------------------------------
// 出力ファイル名
// ------------------------------------------------------------
static std::string MakeOutPath(long long nacc)
{
    std::ostringstream oss;
    oss << kOutDir << "/acc_" << nacc << ".dat";
    return oss.str();
}

int main(int argc, char** argv)
{
    if (argc != 2) {
        std::cerr
            << "Usage: " << argv[0] << " <n_acc>\n"
            << "Output: " << kOutDir << "/acc_<n_acc>.dat\n";
        return 1;
    }

    long long nacc = 0;
    if (!ParseLL(argv[1], nacc)) {
        std::cerr << "ERROR: failed to parse n_acc (non-negative integer required)\n";
        return 1;
    }

    try {
        std::filesystem::create_directories(kOutDir);
    } catch (...) {
        std::cerr << "ERROR: failed to create directory: " << kOutDir << "\n";
        return 1;
    }

    const std::string outpath = MakeOutPath(nacc);

    // RNG（引数なしなので自動seed）
    const uint64_t seed =
        static_cast<uint64_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    TRandom3 rng(static_cast<ULong_t>(seed));
    TRandom3 rng_scan(static_cast<ULong_t>(seed ^ 0x9e3779b97f4a7c15ULL));

    const double Pmu = detres.P_mu;

    // ---- pmax 初期推定（0でもOK）----
    double pmax_michel = 0.0;
    double pmax_rmd6   = 0.0;
    if (nacc > 0) {
        pmax_michel = EstimatePMaxMichel2D(rng_scan, kScanTrialsMichel, Pmu);
        pmax_rmd6   = EstimatePMaxRmd6_4D(rng_scan, kScanTrialsRmd6, Pmu);
    }

    // ---- t 一様 ----
    const double tmin = analysis_window.t_min;
    const double tmax = analysis_window.t_max;

    // ---- 生成 ----
    std::vector<Event> out;
    out.reserve(static_cast<size_t>(nacc));

    long long tries_e = 0;
    long long tries_g = 0;

    auto last_print = std::chrono::steady_clock::now();
    PrintProgressBar("acc", 0, nacc, 0, 0, pmax_michel, pmax_rmd6, 0.0);

    while (static_cast<long long>(out.size()) < nacc) {
        // e 側を窓内（観測）になるまで作る
        double Ee_true = 0.0, phi_e = 0.0;
        double Ee_obs  = std::numeric_limits<double>::quiet_NaN();
        while (true) {
            SampleMichel2D_AR(rng, Pmu, pmax_michel, Ee_true, phi_e, tries_e);

            Ee_obs = smear_energy_trandom3_e(rng, Ee_true);
            if (std::isfinite(Ee_obs) &&
                Ee_obs >= analysis_window.Ee_min && Ee_obs <= analysis_window.Ee_max) {
                break;
            }
        }

        // γ 側を窓内（観測）になるまで作る
        double Eg_true = 0.0, phi_g = 0.0;
        double Eg_obs  = std::numeric_limits<double>::quiet_NaN();
        while (true) {
            SampleGammaFromRmd6_AR(rng, Pmu, pmax_rmd6, Eg_true, phi_g, tries_g);

            Eg_obs = smear_energy_trandom3_g(rng, Eg_true);
            if (std::isfinite(Eg_obs) &&
                Eg_obs >= analysis_window.Eg_min && Eg_obs <= analysis_window.Eg_max) {
                break;
            }
        }

        int idx_e = Detector_PhiIndexFromValue(phi_e, detres.phi_e_min, detres.phi_e_max,
                                               Math_GetNPhiE(detres));
        int idx_g = Detector_PhiIndexFromValue(phi_g, detres.phi_g_min, detres.phi_g_max,
                                               Math_GetNPhiG(detres));
        if (idx_e < 0 || idx_g < 0 || !Detector_IsAllowedPhiPairIndex(idx_e, idx_g, detres)) {
            continue;
        }

        // t 一様
        const double t = rng.Uniform(tmin, tmax);

        Event ev;
        ev.Ee = Ee_obs;
        ev.Eg = Eg_obs;
        ev.t  = t;
        ev.phi_detector_e = phi_e;
        ev.phi_detector_g = phi_g;
        out.push_back(ev);

        // 進捗更新
        const auto now = std::chrono::steady_clock::now();
        const auto dt_ms = std::chrono::duration_cast<std::chrono::milliseconds>(now - last_print).count();
        if (dt_ms >= kProgressIntervalMs || static_cast<long long>(out.size()) == nacc) {
            const long long done = static_cast<long long>(out.size());
            const double acc_rate = (done > 0) ? (double(done) / double(std::max<long long>(1, done))) : 0.0;
            PrintProgressBar("acc", done, nacc, tries_e, tries_g, pmax_michel, pmax_rmd6, acc_rate);
            last_print = now;
        }
    }

    Info("make_acc_mockdata", "generated acc: N=%lld tries_e=%lld tries_g=%lld pmax_michel=%.6g pmax_rmd6=%.6g",
         nacc, tries_e, tries_g, pmax_michel, pmax_rmd6);

    // ---- 書き出し（5列）----
    std::ofstream fout(outpath);
    if (!fout) {
        std::cerr << "ERROR: failed to open output: " << outpath << "\n";
        return 1;
    }

    fout << "# acc mockdata generated from Michel(2D AR) + RMD_d6(4D AR -> gamma singles)\n";
    fout << "# t is uniform in analysis window\n";
    fout << "# energies are smeared by energy_response_shape_e/g(E_res,E_true)\n";
    fout << "# phi_detector_e/g are snapped to detector phi grid at output\n";
    fout << "# n_acc=" << nacc << " seed=" << seed << " P_mu=" << detres.P_mu << "\n";
    fout << "Ee\tEg\tt\tphi_detector_e\tphi_detector_g\n";

    fout << std::setprecision(15);
    for (const auto& ev : out) {
        const double phi_e_out = SnapPhiToGrid(ev.phi_detector_e,
                                               detres.phi_e_min, detres.phi_e_max,
                                               Math_GetNPhiE(detres));
        const double phi_g_out = SnapPhiToGrid(ev.phi_detector_g,
                                               detres.phi_g_min, detres.phi_g_max,
                                               Math_GetNPhiG(detres));
        fout << ev.Ee << "\t"
             << ev.Eg << "\t"
             << ev.t  << "\t"
             << phi_e_out << "\t"
             << phi_g_out << "\n";
    }
    fout.close();

    std::cout << "wrote: " << outpath << "\n";
    std::cout << "generated: n_acc=" << nacc << " total=" << nacc << "\n";
    return 0;
}
