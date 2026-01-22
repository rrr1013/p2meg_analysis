// scripts/make_pdfmix_mockdata.cc
//
// ============================================================
// p2MEG: PDF から擬似データ（.dat）を生成する（t は後付け）
//
// 仕様:
//  - 引数: nsig nrmd のみ（個数固定）
//  - 出力: data/mockdata/pdfmix_<nsig>_<nrmd>.dat
//  - 出力列(5列): Ee  Eg  t  phi_detector_e  phi_detector_g
//
// 生成方針:
//  - (Ee, Eg, phi_e, phi_g) は 4D 空間で一様提案 → PdfComponent.eval を用いた棄却法
//      phi_detector_e/g は DetectorResolution の格子点を一様に提案
//  - t は棄却後に、窓内トランケート正規（detres.t_mean, detres.sigma_t）で後付け
//  - 出力時に phi_detector_e/g は DetectorResolution の格子点へ丸める
//
// 進捗表示:
//  - sig / rmd それぞれについて進捗バーを stderr に表示
//
// 自動設定:
//  - 解析窓: analysis_window
//  - 分解能: detres
//  - 質量:   kMassesPDG
//  - RMD 格子: data/pdf_cache/rmd_grid_gen.root, key="rmd_grid"
//
// 実行例:
//   ./build/make_pdfmix_mockdata 3000 0
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
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "TError.h"

#include "../include/p2meg/Event.h"
#include "../include/p2meg/AnalysisWindow.h"
#include "../include/p2meg/DetectorResolution.h"
#include "../include/p2meg/Constants.h"
#include "../include/p2meg/MathUtils.h"
#include "../include/p2meg/PdfWrappers.h"
#include "../include/p2meg/RMDGridPdf.h"

// ---- 自動設定（必要ならここだけ変更）----
static const char* kOutDir = "data/mockdata";
static const char* kRmdGridFile = "data/pdf_cache/rmd_grid_gen.root";
static const char* kRmdGridKey  = "rmd_grid";

// pmax 初期推定のプリスキャン回数（0でも動くが、あると少し速い）
static constexpr int kScanTrials = 200000;

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
// 標準正規 CDF Φ(z)
// ------------------------------------------------------------
static double StdNormCDF(double z)
{
    static const double inv_sqrt2 = 1.0 / std::sqrt(2.0);
    return 0.5 * (1.0 + std::erf(z * inv_sqrt2));
}

// ------------------------------------------------------------
// 窓内トランケート正規 PDF（窓内で正規化）
// ------------------------------------------------------------
static double TruncNormPdf(double x, double mu, double sigma, double a, double b)
{
    if (!(sigma > 0.0)) return 0.0;
    if (x < a || x > b) return 0.0;

    const double za = (a - mu) / sigma;
    const double zb = (b - mu) / sigma;
    const double Z  = StdNormCDF(zb) - StdNormCDF(za);
    if (!(Z > 0.0)) return 0.0;

    const double z = (x - mu) / sigma;
    const double phi = std::exp(-0.5 * z * z) / (std::sqrt(2.0 * std::acos(-1.0)) * sigma);
    return phi / Z;
}

// ------------------------------------------------------------
// 窓内トランケート正規サンプル（単純リジェクト）
// ------------------------------------------------------------
static double TruncNormSample(std::mt19937_64& rng, double mu, double sigma, double a, double b)
{
    std::normal_distribution<double> nd(mu, sigma);
    while (true) {
        const double x = nd(rng);
        if (x >= a && x <= b) return x;
    }
}

// ------------------------------------------------------------
// (Ee,Eg,phi_e,phi_g) の一様提案
// ------------------------------------------------------------
static Event ProposeUniform4(std::mt19937_64& rng)
{
    std::uniform_real_distribution<double> uEe(analysis_window.Ee_min, analysis_window.Ee_max);
    std::uniform_real_distribution<double> uEg(analysis_window.Eg_min, analysis_window.Eg_max);
    const int N_phi_e = Math_GetNPhiE(detres);
    const int N_phi_g = Math_GetNPhiG(detres);
    const double step_e = Detector_PhiStep(detres.phi_e_min, detres.phi_e_max, N_phi_e);
    const double step_g = Detector_PhiStep(detres.phi_g_min, detres.phi_g_max, N_phi_g);
    std::uniform_int_distribution<int> ui_e(0, N_phi_e);
    std::uniform_int_distribution<int> ui_g(0, N_phi_g);

    while (true) {
        const int ie = ui_e(rng);
        const int ig = ui_g(rng);
        if (!Detector_IsAllowedPhiPairIndex(ie, ig, detres)) continue;

        Event ev;
        ev.Ee = uEe(rng);
        ev.Eg = uEg(rng);
        ev.t  = detres.t_mean; // ダミー（評価時に t0 へ上書き）
        ev.phi_detector_e = detres.phi_e_min + step_e * static_cast<double>(ie);
        ev.phi_detector_g = detres.phi_g_min + step_g * static_cast<double>(ig);
        return ev;
    }
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
                             long long tries,
                             double pmax,
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
        << " tries=" << tries
        << " acc=" << std::setprecision(3) << (acc_rate * 100.0) << "%"
        << " pmax4=" << std::setprecision(6) << pmax;

    std::cerr << oss.str() << std::flush;

    if (done >= total) {
        std::cerr << "\n";
    }
}

// ------------------------------------------------------------
// comp.eval には t が入っているので、t=t0 に固定して評価し Pt(t0) で割って
// 4D 部分の“密度”に戻す（p(other)*Pt(t) で因子分解している前提）
// ------------------------------------------------------------
static double EvalNoT4D(const PdfComponent& comp, const Event& ev4, double t0, double pt0)
{
    if (!(pt0 > 0.0) || !std::isfinite(pt0)) return 0.0;

    Event ev = ev4;
    ev.t = t0;

    const double p = comp.eval(ev, comp.ctx);
    if (!std::isfinite(p) || p <= 0.0) return 0.0;

    return p / pt0;
}

// ------------------------------------------------------------
// 4D の pmax をプリスキャンで推定（0でも可）
// ------------------------------------------------------------
static double EstimatePMax4D(const PdfComponent& comp, std::mt19937_64& rng, int trials,
                            double t0, double pt0)
{
    double pmax = 0.0;
    for (int i = 0; i < trials; ++i) {
        const Event ev4 = ProposeUniform4(rng);
        const double p4 = EvalNoT4D(comp, ev4, t0, pt0);
        if (std::isfinite(p4) && p4 > pmax) pmax = p4;
    }
    if (!(pmax > 0.0)) return 0.0;
    return pmax * kSafetyFactor;
}

// ------------------------------------------------------------
// 4D 棄却法で N 個生成（上限なしで最後まで回す）
//  - pmax が 0 の場合: 最初に p4>0 を見つけた時点で pmax を初期化して続行
//  - p4>pmax の場合: pmax をその場で引き上げて続行（再スタートしない）
//  - 進捗バー付き
// ------------------------------------------------------------
static void Generate4D_AcceptReject_NoStop(const PdfComponent& comp,
                                          long long N,
                                          std::mt19937_64& rng,
                                          std::vector<Event>& out4,
                                          double& pmax_inout,
                                          double t0,
                                          double pt0,
                                          const char* label)
{
    out4.clear();
    out4.reserve(static_cast<size_t>(N));
    if (N == 0) return;

    std::uniform_real_distribution<double> u01(0.0, 1.0);

    long long tries = 0;
    auto last_print = std::chrono::steady_clock::now();

    PrintProgressBar(label, 0, N, 0, pmax_inout, 0.0);

    while (static_cast<long long>(out4.size()) < N) {
        ++tries;

        const Event ev4 = ProposeUniform4(rng);
        const double p4 = EvalNoT4D(comp, ev4, t0, pt0);
        if (!std::isfinite(p4) || p4 <= 0.0) {
            // 進捗更新だけ（一定間隔）
        } else {
            // pmax がまだ無いなら初期化して続行（この試行では受理判定しない＝余分に捨てるだけ）
            if (!(pmax_inout > 0.0) || !std::isfinite(pmax_inout)) {
                pmax_inout = p4 * kSafetyFactor;
            } else {
                // pmax が不足していたら更新して続行（再スタートしない）
                if (p4 > pmax_inout) {
                    pmax_inout = p4 * kUpdateFactor;
                }

                const double accept_prob = p4 / pmax_inout;
                const double u = u01(rng);
                if (u < accept_prob) {
                    out4.push_back(ev4);
                }
            }
        }

        // 進捗更新（一定間隔 or 完了時）
        const auto now = std::chrono::steady_clock::now();
        const auto dt_ms = std::chrono::duration_cast<std::chrono::milliseconds>(now - last_print).count();
        if (dt_ms >= kProgressIntervalMs || static_cast<long long>(out4.size()) == N) {
            const double acc_rate = (tries > 0) ? (double(out4.size()) / double(tries)) : 0.0;
            PrintProgressBar(label, static_cast<long long>(out4.size()), N, tries, pmax_inout, acc_rate);
            last_print = now;
        }
    }

    const double acc = (tries > 0) ? (double(N) / double(tries)) : 0.0;
    Info("make_pdfmix_mockdata", "generated 4D: comp=%s N=%lld tries=%lld acc=%.6g pmax4=%.6g",
         comp.name, N, tries, acc, pmax_inout);
}

// ------------------------------------------------------------
// 出力ファイル名
// ------------------------------------------------------------
static std::string MakeOutPath(long long nsig, long long nrmd)
{
    std::ostringstream oss;
    oss << kOutDir << "/pdfmix_" << nsig << "_" << nrmd << ".dat";
    return oss.str();
}

int main(int argc, char** argv)
{
    if (argc != 3) {
        std::cerr
            << "Usage: " << argv[0] << " <nsig> <nrmd>\n"
            << "Output: " << kOutDir << "/pdfmix_<nsig>_<nrmd>.dat\n";
        return 1;
    }

    long long nsig = 0, nrmd = 0;
    if (!ParseLL(argv[1], nsig) || !ParseLL(argv[2], nrmd)) {
        std::cerr << "ERROR: failed to parse nsig/nrmd (non-negative integers required)\n";
        return 1;
    }

    try {
        std::filesystem::create_directories(kOutDir);
    } catch (...) {
        std::cerr << "ERROR: failed to create directory: " << kOutDir << "\n";
        return 1;
    }

    const std::string outpath = MakeOutPath(nsig, nrmd);

    // RNG（引数なしなので自動seed）
    const uint64_t seed =
        static_cast<uint64_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    std::mt19937_64 rng(seed);
    std::mt19937_64 rng_scan(seed ^ 0x9e3779b97f4a7c15ULL);

    // ---- RMD 格子ロード（nrmd>0 のとき必須）----
    if (nrmd > 0) {
        const bool ok = RMDGridPdf_Load(kRmdGridFile, kRmdGridKey);
        if (!ok || !RMDGridPdf_IsLoaded()) {
            std::cerr << "ERROR: failed to load RMD grid pdf: file=" << kRmdGridFile
                      << " key=" << kRmdGridKey << "\n";
            return 1;
        }
    }

    // ---- PdfComponent 準備 ----
    SignalPdfContext sigctx;
    sigctx.win = analysis_window;
    sigctx.res = detres;
    sigctx.ms  = kMassesPDG;

    PdfComponent sig = MakeSignalComponent(&sigctx);
    PdfComponent rmd = MakeRMDComponent();

    // ---- t 後付けのための Pt(t) ----
    // t0 は窓内にクランプ（pt0>0 を確実にする）
    double t0 = detres.t_mean;
    if (t0 < analysis_window.t_min) t0 = analysis_window.t_min;
    if (t0 > analysis_window.t_max) t0 = analysis_window.t_max;

    const double pt0 = TruncNormPdf(t0, detres.t_mean, detres.sigma_t,
                                   analysis_window.t_min, analysis_window.t_max);
    if (!(pt0 > 0.0) || !std::isfinite(pt0)) {
        std::cerr << "ERROR: pt0 invalid. check detres.sigma_t / t window.\n";
        return 1;
    }

    // ---- pmax(4D) 初期推定（0でもOK）----
    double pmax_sig4 = 0.0;
    double pmax_rmd4 = 0.0;

    if (nsig > 0) pmax_sig4 = EstimatePMax4D(sig, rng_scan, kScanTrials, t0, pt0);
    if (nrmd > 0) pmax_rmd4 = EstimatePMax4D(rmd, rng_scan, kScanTrials, t0, pt0);

    // ---- 4D 生成（止めずに最後まで回す）----
    std::vector<Event> ev_sig4;
    std::vector<Event> ev_rmd4;

    if (nsig > 0) {
        Generate4D_AcceptReject_NoStop(sig, nsig, rng, ev_sig4, pmax_sig4, t0, pt0, "sig");
    }
    if (nrmd > 0) {
        Generate4D_AcceptReject_NoStop(rmd, nrmd, rng, ev_rmd4, pmax_rmd4, t0, pt0, "rmd");
    }

    // ---- t を後付け（signal/rmd とも同じ Pt）----
    for (auto& ev : ev_sig4) {
        ev.t = TruncNormSample(rng, detres.t_mean, detres.sigma_t, analysis_window.t_min, analysis_window.t_max);
    }
    for (auto& ev : ev_rmd4) {
        ev.t = TruncNormSample(rng, detres.t_mean, detres.sigma_t, analysis_window.t_min, analysis_window.t_max);
    }

    // ---- 混合してシャッフル ----
    std::vector<Event> all;
    all.reserve(static_cast<size_t>(nsig + nrmd));
    all.insert(all.end(), ev_sig4.begin(), ev_sig4.end());
    all.insert(all.end(), ev_rmd4.begin(), ev_rmd4.end());
    std::shuffle(all.begin(), all.end(), rng);

    // ---- 書き出し（5列）----
    std::ofstream fout(outpath);
    if (!fout) {
        std::cerr << "ERROR: failed to open output: " << outpath << "\n";
        return 1;
    }

    fout << "# pdfmix mockdata generated from PdfComponent (AR in 4D + t post attach)\n";
    fout << "# phi_detector_e/g are snapped to detector phi grid at output\n";
    fout << "# nsig=" << nsig << " nrmd=" << nrmd << " seed=" << seed << "\n";
    if (nrmd > 0) {
        fout << "# rmd_grid: file=" << kRmdGridFile << " key=" << kRmdGridKey << "\n";
    }
    fout << "Ee\tEg\tt\tphi_detector_e\tphi_detector_g\n";

    fout << std::setprecision(15);
    for (const auto& ev : all) {
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
    std::cout << "generated: nsig=" << nsig << " nrmd=" << nrmd << " total=" << (nsig + nrmd) << "\n";
    return 0;
}
