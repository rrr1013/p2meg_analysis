#include "p2meg/MakeRMDGridPdf.h"

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TNamed.h"
#include "TParameter.h"
#include "TRandom3.h"
#include "THn.h"
#include "TAxis.h"

#include "p2meg/AnalysisWindow.h"
#include "p2meg/DetectorResolution.h"
#include "p2meg/RMDSpectrum.h"

//============================================================
// 内部設定
//
// 単位：Ee, Eg は MeV、t は ns、theta は rad
//============================================================

// ---- 3D格子ビニング（Ee, Eg, theta）----
static constexpr int kNBins_Ee = 40;
static constexpr int kNBins_Eg = 40;
static constexpr int kNBins_th = 30;

// ---- 生成統計 ----
static constexpr long kNTruthSamples  = 1000000L; // 真値サンプル数（Ee,Eg を一様）
static constexpr int  kNSmearPerTruth = 50;       // 1真値あたりのスメア回数
static constexpr unsigned long kSeed  = 20251216UL;

// ---- RMD 理論関数の d_min ----
static constexpr double kDMin = 1e-6;

// ---- cos の許容範囲推定（粗い走査）----
static constexpr int kCosScanN = 121; // [-1,1] をこの点数で走査（奇数推奨）
static constexpr double kCosMargin = 0.0;

//============================================================
// 内部補助関数
//============================================================

static bool IsFinite(double x) {
    return std::isfinite(x);
}

static double Clamp(double x, double lo, double hi) {
    if (x < lo) return lo;
    if (x > hi) return hi;
    return x;
}

static bool IsInsideWindow3D(double Ee, double Eg, double theta) {
    if (Ee < analysis_window_rmd.Ee_min || Ee > analysis_window_rmd.Ee_max) return false;
    if (Eg < analysis_window_rmd.Eg_min || Eg > analysis_window_rmd.Eg_max) return false;
    if (theta < analysis_window_rmd.theta_min || theta > analysis_window_rmd.theta_max) return false;
    return true;
}

// theta の上限 pi 境界歪みを避けるためのスメア
// δ = pi - theta（δ>=0）でガウシアン→折り返し→theta に戻す
static double SmearThetaWithReflection(double theta_true, double sigma_theta, TRandom3& rng) {
    const double pi = 3.14159265358979323846;

    const double delta_true = pi - theta_true; // >=0 を想定
    double delta_obs = delta_true + rng.Gaus(0.0, sigma_theta);

    // 物理範囲 δ>=0 に折り返し
    if (delta_obs < 0.0) delta_obs = -delta_obs;

    // 極端に delta_obs が pi を超えると theta_obs<0 になるが、その場合は窓チェックで落ちる
    return pi - delta_obs;
}

// 与えた (Ee,Eg) に対して w>0 になる cos 範囲を粗く推定する。
// 返り値：範囲が見つかれば true、無ければ false。
static bool EstimateAllowedCosRange(double Ee, double Eg, double& cos_min, double& cos_max) {
    bool found = false;
    double cmin = 0.0, cmax = 0.0;

    for (int i = 0; i < kCosScanN; ++i) {
        const double u = (kCosScanN == 1) ? 0.0 : static_cast<double>(i) / (kCosScanN - 1);
        const double c = -1.0 + 2.0 * u; // [-1,1]
        const double w = RMD_d3B_dEe_dEg_dcos(Ee, Eg, c, kDMin);
        if (w > 0.0 && IsFinite(w)) {
            if (!found) {
                found = true;
                cmin = c;
                cmax = c;
            } else {
                if (c < cmin) cmin = c;
                if (c > cmax) cmax = c;
            }
        }
    }

    if (!found) return false;

    cmin -= kCosMargin;
    cmax += kCosMargin;

    cmin = Clamp(cmin, -1.0, 1.0);
    cmax = Clamp(cmax, -1.0, 1.0);
    if (!(cmax > cmin)) return false;

    cos_min = cmin;
    cos_max = cmax;
    return true;
}

// THnD（3D）の全ビン総和
static double SumAllBins3(const THnD& h) {
    const int n0 = h.GetAxis(0)->GetNbins();
    const int n1 = h.GetAxis(1)->GetNbins();
    const int n2 = h.GetAxis(2)->GetNbins();

    std::vector<int> idx(3, 1);
    double sum = 0.0;

    for (int i0 = 1; i0 <= n0; ++i0) {
        idx[0] = i0;
        for (int i1 = 1; i1 <= n1; ++i1) {
            idx[1] = i1;
            for (int i2 = 1; i2 <= n2; ++i2) {
                idx[2] = i2;
                const Long64_t bin = h.GetBin(idx.data());
                sum += h.GetBinContent(bin);
            }
        }
    }
    return sum;
}

// THnD（3D）を密度に変換して、解析窓内で積分=1に正規化する
// density3 = (C / total_mass) / volume3(bin)
static int ConvertToDensityAndNormalize3(THnD& h) {
    const TAxis* ax0 = h.GetAxis(0);
    const TAxis* ax1 = h.GetAxis(1);
    const TAxis* ax2 = h.GetAxis(2);

    const int n0 = ax0->GetNbins();
    const int n1 = ax1->GetNbins();
    const int n2 = ax2->GetNbins();

    const double total_mass = SumAllBins3(h);
    if (!(total_mass > 0.0) || !IsFinite(total_mass)) return 1;

    std::vector<int> idx(3, 1);

    for (int i0 = 1; i0 <= n0; ++i0) {
        idx[0] = i0;
        const double w0 = ax0->GetBinWidth(i0);
        for (int i1 = 1; i1 <= n1; ++i1) {
            idx[1] = i1;
            const double w1 = ax1->GetBinWidth(i1);
            for (int i2 = 1; i2 <= n2; ++i2) {
                idx[2] = i2;
                const double w2 = ax2->GetBinWidth(i2);

                const double vol = w0 * w1 * w2; // [MeV^2 * rad]
                if (!(vol > 0.0) || !IsFinite(vol)) continue;

                const Long64_t bin = h.GetBin(idx.data());
                const double C = h.GetBinContent(bin);

                double density = 0.0;
                if (C > 0.0 && IsFinite(C)) {
                    density = (C / total_mass) / vol;
                }
                h.SetBinContent(bin, density);
            }
        }
    }
    return 0;
}

static std::string BuildMetaString(long n_range_ok, long n_cos_ok, long n_filled3) {
    std::ostringstream oss;
    oss << "MakeRMDGridPdf meta (3D grid only)\n";
    oss << "bins: Ee=" << kNBins_Ee << ", Eg=" << kNBins_Eg << ", theta=" << kNBins_th << "\n";
    oss << "truth_samples=" << kNTruthSamples << "\n";
    oss << "smear_per_truth=" << kNSmearPerTruth << "\n";
    oss << "seed=" << kSeed << "\n";
    oss << "d_min=" << kDMin << "\n";
    oss << "cos_scan_N=" << kCosScanN << "\n";
    oss << "cos_range_ok=" << n_range_ok << "\n";
    oss << "cos_wpos_ok=" << n_cos_ok << "\n";
    oss << "filled_entries_3d=" << n_filled3 << "\n";

    oss << "window: Ee=[" << analysis_window_rmd.Ee_min << "," << analysis_window_rmd.Ee_max << "] MeV\n";
    oss << "window: Eg=[" << analysis_window_rmd.Eg_min << "," << analysis_window_rmd.Eg_max << "] MeV\n";
    oss << "window: t=[" << analysis_window_rmd.t_min << "," << analysis_window_rmd.t_max << "] ns (applied at evaluation)\n";
    oss << "window: theta=[" << analysis_window_rmd.theta_min << "," << analysis_window_rmd.theta_max << "] rad\n";

    oss << "res: sigma_Ee=" << detres_rmd.sigma_Ee << " MeV\n";
    oss << "res: sigma_Eg=" << detres_rmd.sigma_Eg << " MeV\n";
    oss << "res: sigma_t=" << detres_rmd.sigma_t << " ns (used analytically at evaluation)\n";
    oss << "res: sigma_theta=" << detres_rmd.sigma_theta << " rad\n";
    oss << "res: t_mean=" << detres_rmd.t_mean << " ns (used analytically at evaluation)\n";

    return oss.str();
}

//============================================================
// 本体
//============================================================

int MakeRMDGridPdf(const char* out_filepath, const char* key) {
    if (!out_filepath || !key) {
        std::cerr << "[MakeRMDGridPdf] invalid arguments\n";
        return 1;
    }

    // soft photon 発散対策：Eg > Eg_min
    if (!(analysis_window_rmd.Eg_min > 0.0)) {
        std::cerr << "[MakeRMDGridPdf] Eg_min must be > 0 to avoid soft photon divergence.\n";
        return 2;
    }

    // ---- 3Dヒスト（Ee, Eg, theta）----
    const int ndim3 = 3;
    int nbins3[ndim3] = {kNBins_Ee, kNBins_Eg, kNBins_th};
    double xmin3[ndim3] = {
        analysis_window_rmd.Ee_min,
        analysis_window_rmd.Eg_min,
        analysis_window_rmd.theta_min
    };
    double xmax3[ndim3] = {
        analysis_window_rmd.Ee_max,
        analysis_window_rmd.Eg_max,
        analysis_window_rmd.theta_max
    };

    THnD h3("rmd_grid_tmp", "RMD smeared grid (3D);Ee;Eg;theta", ndim3, nbins3, xmin3, xmax3);
    h3.Sumw2();
    h3.GetAxis(0)->SetTitle("Ee [MeV]");
    h3.GetAxis(1)->SetTitle("Eg [MeV]");
    h3.GetAxis(2)->SetTitle("theta [rad]");

    TRandom3 rng(kSeed);

    const double Ee_min = analysis_window_rmd.Ee_min;
    const double Ee_max = analysis_window_rmd.Ee_max;
    const double Eg_min = analysis_window_rmd.Eg_min;
    const double Eg_max = analysis_window_rmd.Eg_max;

    const double sigma_Ee = detres_rmd.sigma_Ee;
    const double sigma_Eg = detres_rmd.sigma_Eg;
    const double sigma_th = detres_rmd.sigma_theta;

    long n_range_ok = 0;   // (Ee,Eg) で cos 範囲が見つかった回数
    long n_cos_ok   = 0;   // その範囲で cos を投げて w>0 だった回数
    long n_filled3  = 0;   // 3D に fill した回数
    long n_range_ng = 0;   // cos 範囲推定に失敗した回数

    for (long it = 0; it < kNTruthSamples; ++it) {
        const double Ee_true = rng.Uniform(Ee_min, Ee_max);
        const double Eg_true = rng.Uniform(Eg_min, Eg_max);

        // 許される cos 範囲を推定し、その範囲でだけ cos をサンプル
        double cmin = 0.0, cmax = 0.0;
        if (!EstimateAllowedCosRange(Ee_true, Eg_true, cmin, cmax)) {
            ++n_range_ng;
            continue;
        }
        ++n_range_ok;

        const double L = cmax - cmin;
        if (!(L > 0.0) || !IsFinite(L)) continue;

        // 範囲内一様
        const double cos_true = rng.Uniform(cmin, cmax);

        // 理論重み（運動学的に許されない領域では 0）
        const double w0 = RMD_d3B_dEe_dEg_dcos(Ee_true, Eg_true, cos_true, kDMin);
        if (!(w0 > 0.0) || !IsFinite(w0)) continue;
        ++n_cos_ok;

        // 重要度補正：proposal が [cmin,cmax] 一様なので、w = w0 / (1/L) = w0 * L
        const double w = w0 * L;

        const double c = Clamp(cos_true, -1.0, 1.0);
        const double theta_true = std::acos(c);

        // 同じ真値点から複数回スメアして観測分布を埋める
        for (int is = 0; is < kNSmearPerTruth; ++is) {
            const double Ee_obs = Ee_true + rng.Gaus(0.0, sigma_Ee);
            const double Eg_obs = Eg_true + rng.Gaus(0.0, sigma_Eg);
            const double theta_obs = SmearThetaWithReflection(theta_true, sigma_th, rng);

            if (!IsInsideWindow3D(Ee_obs, Eg_obs, theta_obs)) continue;

            double x3[3] = {Ee_obs, Eg_obs, theta_obs};
            h3.Fill(x3, w);
            ++n_filled3;
        }
    }

    std::cout << "[MakeRMDGridPdf] (Ee,Eg) with allowed cos-range: " << n_range_ok
              << " / " << kNTruthSamples << "\n";
    std::cout << "[MakeRMDGridPdf] cos sampled with w>0: " << n_cos_ok
              << " / " << n_range_ok << "\n";
    std::cout << "[MakeRMDGridPdf] filled entries (3D): " << n_filled3 << "\n";
    std::cout << "[MakeRMDGridPdf] cos-range estimation failed: " << n_range_ng << "\n";

    // 正規化（∫ pdf3 dEe dEg dtheta = 1）
    const int norm3 = ConvertToDensityAndNormalize3(h3);
    if (norm3 != 0) {
        std::cerr << "[MakeRMDGridPdf] normalization (3D) failed\n";
        return 3;
    }

    // ---- 保存 ----
    TFile fout(out_filepath, "RECREATE");
    if (fout.IsZombie()) {
        std::cerr << "[MakeRMDGridPdf] cannot open output file: " << out_filepath << "\n";
        return 4;
    }

    h3.SetName(key);
    h3.Write(key);

    const std::string meta = BuildMetaString(n_range_ok, n_cos_ok, n_filled3);
    TNamed meta_obj((std::string(key) + "_meta").c_str(), meta.c_str());
    meta_obj.Write();

    TParameter<double> par_dmin((std::string(key) + "_d_min").c_str(), kDMin);
    par_dmin.Write();

    TParameter<Long64_t> par_seed((std::string(key) + "_seed").c_str(),
                                  static_cast<Long64_t>(kSeed));
    par_seed.Write();

    fout.Close();

    std::cout << "[MakeRMDGridPdf] saved (3D): " << out_filepath << " (key=" << key << ")\n";
    return 0;
}
