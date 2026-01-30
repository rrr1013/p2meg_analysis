// scripts/make_michel_E_mockdata.cc
//
// Michel 偏極測定（別実験）用 Ee-only 模擬データ生成（露光=同一、率に整合）
//
// 入力（コマンドライン引数）:
//   1) detector_angle_deg : 検出器を置く角度 [deg]
//   2) output_path        : 出力ファイルパス
//
// 出力フォーマット:
//   - 1行1事象、1列: Ee [MeV]
//
// 生成モデル:
//   - 真の形: Michel_d2Shape_dE_dCosTheta(E_true, cosθ, P_mu_true)
//   - P_mu_true は detres.P_mu を使用
//   - エネルギー応答: energy_response_shape_e(E_obs, E_true)
//   - E_obs は michel_pol_config の Ee_min..Ee_max の範囲内に生成
//
// 改良点（重要）:
//   - 角度ごとに固定イベント数を出すのをやめる
//   - 畳み込み済みの「率」 w_i(θ) を作り、総率 sumw(θ)=Σ_i w_i(θ) に比例して
//     期待事象数 N(θ) を決める（同一露光を模擬）
//   - N(θ) はポアソン揺らぎを入れる
//   - 参照角（cosθ=0）で平均 N_ref になるように露光スケールを内部で決める
//
// 注意:
//   - Run1/Run2 を別データとして作りたい場合は output_path を変える（seed が変わる）
//

#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "TRandom3.h"

#include "p2meg/DetectorResolution.h"
#include "p2meg/MichelPolConfig.h"
#include "p2meg/MichelSpectrum.h"
#include "p2meg/Constants.h"

static inline bool IsFinite(double x) { return std::isfinite(x); }
static inline double DegToRad(double deg) { return deg * (pi / 180.0); }

// 文字列から安定な 32-bit seed を作る（FNV-1a）
static uint32_t Hash32FNV1a(const std::string& s)
{
    uint32_t h = 2166136261u;
    for (unsigned char c : s) {
        h ^= static_cast<uint32_t>(c);
        h *= 16777619u;
    }
    // ROOT TRandom3 は seed=0 を特別扱い（時間依存）するので避ける
    if (h == 0u) h = 1u;
    return h;
}

static double IntegrateResponseShapeOverInterval(double E_true,
                                                 double E_low,
                                                 double E_high,
                                                 int n_steps)
{
    if (!(E_high > E_low)) return 0.0;
    if (n_steps < 1) n_steps = 1;

    const double h = (E_high - E_low) / static_cast<double>(n_steps);
    double sum = 0.0;
    for (int k = 0; k < n_steps; ++k) {
        const double Eobs = E_low + (k + 0.5) * h;
        const double sh = energy_response_shape_e(Eobs, E_true);
        if (IsFinite(sh) && sh > 0.0) sum += sh;
    }
    return sum * h;
}

// 観測ビン重み w_i(θ) を作る（率に比例、正規化しない）
static std::vector<double> BuildObservedRateBinned(const MichelPolConfig& cfg,
                                                   double cos_theta,
                                                   double P_mu_true,
                                                   double& sumw_out)
{
    sumw_out = 0.0;

    std::vector<double> w;
    if (!(cfg.Ee_max > cfg.Ee_min)) return w;
    if (!(cfg.nbins_Ee > 0)) return w;

    const int nb = cfg.nbins_Ee;
    w.assign(nb, 0.0);

    const double binw = (cfg.Ee_max - cfg.Ee_min) / static_cast<double>(nb);

    // 数値積分設定
    const double Etrue_min = cfg.Ee_min;
    const double Etrue_max = cfg.Ee_max;
    const int N_true = 2500;
    const double dE = (Etrue_max - Etrue_min) / static_cast<double>(N_true);

    const int N_obs_in_bin = 6;

    std::vector<double> Ibin(nb, 0.0);

    for (int it = 0; it < N_true; ++it) {
        const double E_true = Etrue_min + (it + 0.5) * dE;
        if (!(E_true >= 0.0) || !IsFinite(E_true)) continue;

        const double s_true = Michel_d2Shape_dE_dCosTheta(E_true, cos_theta, P_mu_true);
        if (!IsFinite(s_true) || s_true <= 0.0) continue;

        // 解析窓内での応答 shape をビン積分し、窓内で正規化した p_bin を使う
        // （MichelPolTemplate.cc 側と同じ扱いにして closure を取る）
        double Z = 0.0;
        for (int j = 0; j < nb; ++j) {
            const double Elow  = cfg.Ee_min + j * binw;
            const double Ehigh = Elow + binw;
            const double I = IntegrateResponseShapeOverInterval(E_true, Elow, Ehigh, N_obs_in_bin);
            Ibin[j] = I;
            Z += I;
        }
        if (!(Z > 0.0) || !IsFinite(Z)) continue;

        for (int j = 0; j < nb; ++j) {
            const double p = Ibin[j] / Z;
            if (!IsFinite(p) || p <= 0.0) continue;
            w[j] += s_true * p * dE;
        }
    }

    // 非負ガード + sumw
    for (int j = 0; j < nb; ++j) {
        if (!IsFinite(w[j]) || w[j] < 0.0) w[j] = 0.0;
        sumw_out += w[j];
    }
    if (!IsFinite(sumw_out) || sumw_out < 0.0) sumw_out = 0.0;

    return w;
}

static bool BuildCdf(const std::vector<double>& w, std::vector<double>& cdf)
{
    cdf.clear();
    double sumw = 0.0;
    for (double x : w) if (IsFinite(x) && x > 0.0) sumw += x;
    if (!(sumw > 0.0) || !IsFinite(sumw)) return false;

    cdf.resize(w.size(), 0.0);
    double acc = 0.0;
    for (size_t i = 0; i < w.size(); ++i) {
        const double wi = (IsFinite(w[i]) && w[i] > 0.0) ? w[i] : 0.0;
        acc += wi / sumw;
        cdf[i] = acc;
    }
    cdf.back() = 1.0;
    return true;
}

static int SampleIndexFromCdf(const std::vector<double>& cdf, double u01)
{
    for (int i = 0; i < (int)cdf.size(); ++i) {
        if (u01 <= cdf[i]) return i;
    }
    return (int)cdf.size() - 1;
}

static void PrintUsage(const char* prog)
{
    std::cout
        << "Usage:\n"
        << "  " << prog << " detector_angle_deg output_path\n"
        << "\n"
        << "Example:\n"
        << "  " << prog << " 120 data/mockdata_michel/run1_A_minus.dat\n";
}

int main(int argc, char** argv)
{
    if (argc != 3) {
        PrintUsage(argv[0]);
        return 2;
    }

    const double theta_deg = std::atof(argv[1]);
    const std::string out_path = argv[2];

    if (!IsFinite(theta_deg)) {
        std::cout << "[make_michel_E_mockdata] invalid angle.\n";
        return 2;
    }

    const auto& cfg = michel_pol_config;

    const double theta_rad = DegToRad(theta_deg);
    const double cos_theta = std::cos(theta_rad);

    const double P_mu_true = detres.P_mu;

    // seed は output_path と angle から決める（run名を変えると独立になる）
    const uint32_t seed = Hash32FNV1a(out_path + "|" + std::to_string(theta_deg));
    TRandom3 rng(seed);

    // 参照角（cosθ=0）での総率から露光スケールを決める
    // 参照角の平均イベント数（目安）
    static constexpr long long N_ref_mean = 200000;

    double sumw_ref = 0.0;
    const std::vector<double> w_ref = BuildObservedRateBinned(cfg, /*cosθ=*/0.0, P_mu_true, sumw_ref);
    if (!(sumw_ref > 0.0)) {
        std::cout << "[make_michel_E_mockdata] failed to compute reference rate.\n";
        return 1;
    }
    const double exposure = (double)N_ref_mean / sumw_ref;  // 露光スケール（任意単位）

    // 指定角での率
    double sumw = 0.0;
    const std::vector<double> w = BuildObservedRateBinned(cfg, cos_theta, P_mu_true, sumw);
    if (!(sumw > 0.0)) {
        std::cout << "[make_michel_E_mockdata] failed to build rate for this angle.\n";
        return 1;
    }

    // 同一露光での期待事象数 -> ポアソン揺らぎ
    const double lambda = exposure * sumw;
    if (!(lambda > 0.0) || !IsFinite(lambda)) {
        std::cout << "[make_michel_E_mockdata] invalid expected count.\n";
        return 1;
    }
    const long long N_events = (long long)rng.PoissonD(lambda);

    std::vector<double> cdf;
    if (!BuildCdf(w, cdf)) {
        std::cout << "[make_michel_E_mockdata] failed to build CDF.\n";
        return 1;
    }

    std::ofstream ofs(out_path);
    if (!ofs) {
        std::cout << "[make_michel_E_mockdata] cannot open output file.\n";
        return 1;
    }

    const double binw = (cfg.Ee_max - cfg.Ee_min) / static_cast<double>(cfg.nbins_Ee);

    std::cout << "=== Michel Ee-only mockdata (rate-consistent) ===\n";
    std::cout << "theta_deg    = " << theta_deg << " [deg]\n";
    std::cout << "cos(theta)   = " << std::setprecision(10) << cos_theta << "\n";
    std::cout << "P_mu_true    = " << P_mu_true << "\n";
    std::cout << "Ee range     = " << cfg.Ee_min << " .. " << cfg.Ee_max << " [MeV]\n";
    std::cout << "nbins_Ee     = " << cfg.nbins_Ee << "\n";
    std::cout << "seed         = " << seed << "\n";
    std::cout << "sumw(angle)  = " << std::setprecision(12) << sumw << "\n";
    std::cout << "sumw(ref)    = " << std::setprecision(12) << sumw_ref << "  (cosθ=0)\n";
    std::cout << "N_expected   = " << std::fixed << std::setprecision(2) << lambda << "\n";
    std::cout << "N_generated  = " << N_events << "\n";
    std::cout << "output_path  = " << out_path << "\n\n";

    ofs << std::fixed << std::setprecision(8);

    const long long report_every = (N_events >= 200000) ? 20000 : 10000;

    for (long long n = 0; n < N_events; ++n) {
        const double u = rng.Uniform();
        const int ib = SampleIndexFromCdf(cdf, u);

        const double Elow  = cfg.Ee_min + ib * binw;
        const double Ehigh = Elow + binw;

        // ビン内一様
        const double Ee_obs = Elow + rng.Uniform() * (Ehigh - Elow);

        ofs << Ee_obs << "\n";

        if (report_every > 0 && (n + 1) % report_every == 0) {
            std::cout << "generated " << (n + 1) << " / " << N_events << "\n";
        }
    }

    std::cout << "\nDone.\n";
    return 0;
}
