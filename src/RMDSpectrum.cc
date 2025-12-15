// src/RMDSpectrum.cc
#include "p2meg/RMDSpectrum.h"
#include "p2meg/Constants.h"

#include <cmath>

//============================================================
// 実装メモ
//  - 停止ミューオン静止系。
//  - 角度の「全体回転」を平均した後の三重微分
//      d^3B/(dx dy dcosθ) = (α/(8π)) * β * (1/y) * F(x,y,d)
//    を実装。
//  - F は Kuno & Okada Appendix の F^(0),F^(1),F^(2) を使用。
//  - 物理定数 α, π は Constants.h の alpha, pi を使用。
//  - r=(me/mmu)^2 のような依存量はこのファイル内で計算する。
//  - 重要：ニュートリノ系の不変質量 q^2=(p_mu - p_e - k)^2 >= 0 を必ず満たす必要がある。
//    これを満たさない (x,y,cosθ) は物理的に不可能なので 0 を返す。
//============================================================

//------------------------------
// r = (me/mmu)^2
//------------------------------
static inline double r_mass_ratio2()
{
    const double mmu = kMassesPDG.m_mu;
    const double me  = kMassesPDG.m_e;
    return (me * me) / (mmu * mmu);
}

//------------------------------
// 運動学領域チェック（x,y）
//   ※これは「ある角度が存在する」ための範囲で、cosθ まで固定したときの領域はさらに狭くなる。
//------------------------------
static bool IsAllowedXY(double x, double y)
{
    if (!(y > 0.0)) return false;

    const double r = r_mass_ratio2();
    const double sqrt_r = std::sqrt(r);

    // y の最大は 1 - r
    if (y > 1.0 - r) return false;

    // x の最大は 1 + r
    const double x_max = 1.0 + r;

    // y によって x_min が変わる
    if (y <= 1.0 - sqrt_r) {
        const double x_min = 2.0 * sqrt_r;
        return (x > x_min && x < x_max);
    } else {
        const double denom = 1.0 - y;
        if (!(denom > 0.0)) return false;
        const double x_min = denom + r / denom;
        return (x >= x_min && x <= x_max);
    }
}

//------------------------------
// β(x) と d
//------------------------------
static double BetaFromX(double x)
{
    // β = sqrt(1 - 4r/x^2)
    const double r = r_mass_ratio2();
    const double t = 1.0 - 4.0 * r / (x * x);
    return (t > 0.0) ? std::sqrt(t) : 0.0;
}

static double d_from(double beta, double cosTheta)
{
    return 1.0 - beta * cosTheta;
}

//------------------------------
// ニュートリノ系の不変質量チェック
//   q = p_mu - p_e - k (k: photon)
//   q^2 >= 0 が必要
//
//   muon rest frame で、x=2Ee/mmu, y=2Eg/mmu, d=1-βcosθ とすると
//     q^2 / mmu^2 = 1 + r - x - y + (x y / 2) d
//------------------------------
static bool IsAllowedQ2(double x, double y, double d)
{
    const double r = r_mass_ratio2();
    const double q2_over_m2 = 1.0 + r - x - y + 0.5 * x * y * d;

    // 数値誤差で -1e-15 程度になることを避けるためのマージン
    return (q2_over_m2 >= -1e-12);
}

//------------------------------
// F^(0), F^(1), F^(2)
//------------------------------
static double F0(double x, double y, double d)
{
    const double x2 = x * x;
    const double x3 = x2 * x;

    const double y2 = y * y;
    const double d2 = d * d;

    const double term1 = (8.0 / d) *
        ( y2 * (3.0 - 2.0 * y)
        + 6.0 * x * y * (1.0 - y)
        + 2.0 * x2 * (3.0 - 4.0 * y)
        - 4.0 * x3 );

    const double term2 = 8.0 *
        ( -x * y * (3.0 - y - y2)
        - x2 * (3.0 - y - 4.0 * y2)
        + 2.0 * x3 * (1.0 + 2.0 * y) );

    const double term3 = 2.0 * d *
        ( x2 * y * (6.0 - 5.0 * y - 2.0 * y2)
        - 2.0 * x3 * y * (4.0 + 3.0 * y) );

    const double term4 = 2.0 * d2 * ( x3 * y2 * (2.0 + y) );

    return term1 + term2 + term3 + term4;
}

static double F1(double x, double y, double d)
{
    const double x2 = x * x;
    const double d2 = d * d;

    const double term1 = (32.0 / d2) *
        ( -(y * (3.0 - 2.0 * y)) / x
          - (3.0 - 4.0 * y)
          + 2.0 * x );

    const double term2 = (8.0 / d) *
        ( y * (6.0 - 5.0 * y)
          - 2.0 * x * (4.0 + y)
          + 6.0 * x2 );

    const double term3 = 8.0 *
        ( x * (4.0 - 3.0 * y + y * y)
          - 3.0 * x2 * (1.0 + y) );

    const double term4 = 6.0 * d * ( x2 * y * (2.0 + y) );

    return term1 + term2 + term3 + term4;
}

static double F2(double x, double y, double d)
{
    const double d2 = d * d;

    const double term1 = (32.0 / d2) * ( (4.0 - 3.0 * y) / x - 3.0 );
    const double term2 = (48.0 * y / d);

    return term1 + term2;
}

static double F(double x, double y, double d)
{
    const double r  = r_mass_ratio2();
    const double r2 = r * r;
    return F0(x, y, d) + r * F1(x, y, d) + r2 * F2(x, y, d);
}

//============================================================
// 公開関数（既存）
//============================================================
double RMD_d3B_dxdy_dcos(double x, double y, double cosTheta, double d_min)
{
    // 入力の基本チェック
    if (!IsAllowedXY(x, y)) return 0.0;
    if (cosTheta < -1.0 || cosTheta > 1.0) return 0.0;

    const double beta = BetaFromX(x);
    if (!(beta > 0.0)) return 0.0;

    double d = d_from(beta, cosTheta);

    // 物理領域チェック：q^2 >= 0
    if (!IsAllowedQ2(x, y, d)) return 0.0;

    // 共線領域の数値ガード（物理カットではない）
    if (d < d_min) d = d_min;

    // soft photon で 1/y 発散するので、呼び出し側で必ず y>y_min を課すこと
    const double f = F(x, y, d);

    // d^3B/(dx dy dcosθ) = (α/(8π)) * β * (1/y) * F
    return (alpha / (8.0 * pi)) * beta * (1.0 / y) * f;
}

double RMD_d3B_dEe_dEg_dcos(double Ee_MeV, double Eg_MeV, double cosTheta, double d_min)
{
    if (!(Ee_MeV > 0.0) || !(Eg_MeV > 0.0)) return 0.0;

    const double mmu = kMassesPDG.m_mu;

    const double x = 2.0 * Ee_MeV / mmu;
    const double y = 2.0 * Eg_MeV / mmu;

    // 変数変換: x=2Ee/mmu, y=2Eg/mmu より
    //   dx dy = (2/mmu)^2 dEe dEg  => d^3B/(dEe dEg dcos) = (4/mmu^2) * d^3B/(dx dy dcos)
    const double val = RMD_d3B_dxdy_dcos(x, y, cosTheta, d_min);
    return (4.0 / (mmu * mmu)) * val;
}

//============================================================
// ここから追加：完全式（F,G,H + 偏極）
//
// Kuno & Okada (Rev.Mod.Phys.73 (2001) 151) Appendix A の式に従い、
//   d^6B/(dx dy dΩ_e dΩ_γ)
//   = (α/(64π^3)) * β * (1/y) * [ F(x,y,d) - β (P·p̂_e) G(x,y,d) - (P·p̂_γ) H(x,y,d) ]
//
// ここでは入力を「系統Ⅱ」:
//   cosThetaEG = p̂_e · p̂_γ
//   cosThetaE  = P̂ · p̂_e
//   cosThetaG  = P̂ · p̂_γ
// とし、Pmu を偏極度（スカラー；符号込み）として掛ける。
//============================================================

//------------------------------
// 角度三つ組の整合性チェック
//  cosΘ_eγ = cosΘ_e cosΘ_γ + sinΘ_e sinΘ_γ cosΔφ
//  したがって
//    |cosΘ_eγ - cosΘ_e cosΘ_γ| <= sinΘ_e sinΘ_γ
// を満たさない入力は幾何学的に不可能なので 0 を返す。
//------------------------------
static bool IsConsistentAngleTriplet(double cosEG, double cosE, double cosG)
{
    if (cosEG < -1.0 || cosEG > 1.0) return false;
    if (cosE  < -1.0 || cosE  > 1.0) return false;
    if (cosG  < -1.0 || cosG  > 1.0) return false;

    const double se2 = 1.0 - cosE * cosE;
    const double sg2 = 1.0 - cosG * cosG;

    // 数値誤差ガード
    double prod = se2 * sg2;
    if (prod < 0.0) prod = 0.0;

    const double rhs = std::sqrt(prod);
    const double lhs = cosEG - cosE * cosG;

    const double tol = 1e-12;
    return (lhs * lhs <= (rhs + tol) * (rhs + tol));
}

//------------------------------
// G^(0), G^(1), G^(2)
//------------------------------
static double G0(double x, double y, double d)
{
    const double x2 = x * x;
    const double x3 = x2 * x;

    const double y2 = y * y;

    const double term1 = (8.0 / d) *
        ( x * y * (1.0 - 2.0 * y)
        + 2.0 * x2 * (1.0 - 3.0 * y)
        - 4.0 * x3 );

    const double term2 = 4.0 *
        ( -x2 * (2.0 - 3.0 * y - 4.0 * y2)
          + 2.0 * x3 * (2.0 + 3.0 * y) );

    const double term3 = -4.0 * d * ( x3 * y * (2.0 + y) );

    return term1 + term2 + term3;
}

static double G1(double x, double y, double d)
{
    const double x2 = x * x;
    const double d2 = d * d;

    const double term1 = (32.0 / d2) * ( -1.0 + 2.0 * y + 2.0 * x );
    const double term2 = (8.0 / d)  * ( -x * y + 6.0 * x2 );
    const double term3 = -12.0 * x2 * (2.0 + y);

    return term1 + term2 + term3;
}

static double G2(double /*x*/, double /*y*/, double d)
{
    const double d2 = d * d;
    return -96.0 / d2;
}

static double G(double x, double y, double d)
{
    const double r  = r_mass_ratio2();
    const double r2 = r * r;
    return G0(x, y, d) + r * G1(x, y, d) + r2 * G2(x, y, d);
}

//------------------------------
// H^(0), H^(1), H^(2)
//------------------------------
static double H0(double x, double y, double d)
{
    const double x2 = x * x;
    const double x3 = x2 * x;

    const double y2 = y * y;
    const double y3 = y2 * y;

    const double d2 = d * d;

    const double term1 = (8.0 / d) *
        ( y2 * (1.0 - 2.0 * y)
        + x * y * (1.0 - 4.0 * y)
        - 2.0 * x2 * y );

    const double term2 = 4.0 *
        ( 2.0 * x * y2 * (1.0 + y)
        - x2 * y * (1.0 - 4.0 * y)
        + 2.0 * x3 * y );

    const double term3 = 2.0 * d *
        ( x2 * y2 * (1.0 - 2.0 * y)
        - 4.0 * x3 * y2 );

    const double term4 = 2.0 * d2 * ( x3 * y3 );

    return term1 + term2 + term3 + term4;
}

static double H1(double x, double y, double d)
{
    const double x2 = x * x;
    const double d2 = d * d;

    const double term1 = (32.0 / d2) *
        ( -(y * (1.0 - 2.0 * y)) / x
          + 2.0 * y );

    const double term2 = (8.0 / d) *
        ( y * (2.0 - 5.0 * y)
          - x * y );

    const double term3 = 4.0 * x * y * (2.0 * y - 3.0 * x);

    const double term4 = 6.0 * d * ( x2 * y * y );

    return term1 + term2 + term3 + term4;
}

static double H2(double x, double y, double d)
{
    const double d2 = d * d;

    const double term1 = -96.0 * y / (d2 * x);
    const double term2 =  48.0 * y / d;

    return term1 + term2;
}

static double H(double x, double y, double d)
{
    const double r  = r_mass_ratio2();
    const double r2 = r * r;
    return H0(x, y, d) + r * H1(x, y, d) + r2 * H2(x, y, d);
}

//============================================================
// 公開関数（追加：完全式）
//============================================================
double RMD_d6B_dxdy_dOmegae_dOmegag(double x, double y,
                                   double cosThetaEG,
                                   double cosThetaE,
                                   double cosThetaG,
                                   double Pmu,
                                   double d_min)
{
    // 入力の基本チェック
    if (!IsAllowedXY(x, y)) return 0.0;
    if (!IsConsistentAngleTriplet(cosThetaEG, cosThetaE, cosThetaG)) return 0.0;

    const double beta = BetaFromX(x);
    if (!(beta > 0.0)) return 0.0;

    const double d0 = d_from(beta, cosThetaEG);

    // 物理領域チェック：q^2 >= 0（チェックは d0 を使う）
    if (!IsAllowedQ2(x, y, d0)) return 0.0;

    // 共線領域の数値ガード（物理カットではない）
    double d = d0;
    if (d < d_min) d = d_min;

    // soft photon で 1/y 発散（実際に分岐比を定義するには Eg > Eg_min が必須）
    const double f = F(x, y, d);
    const double g = G(x, y, d);
    const double h = H(x, y, d);

    // d^6B/(dx dy dΩe dΩγ) = (α/(64π^3)) * β * (1/y) * [ F - β(P·p̂e)G - (P·p̂γ)H ]
    // ここでは (P·p̂e)=Pmu*cosThetaE, (P·p̂γ)=Pmu*cosThetaG とする。
    const double pref = (alpha / (64.0 * pi * pi * pi)) * beta * (1.0 / y);
    const double bracket = f - beta * (Pmu * cosThetaE) * g - (Pmu * cosThetaG) * h;

    return pref * bracket;
}

double RMD_d6B_dEe_dEg_dOmegae_dOmegag(double Ee_MeV, double Eg_MeV,
                                      double cosThetaEG,
                                      double cosThetaE,
                                      double cosThetaG,
                                      double Pmu,
                                      double d_min)
{
    if (!(Ee_MeV > 0.0) || !(Eg_MeV > 0.0)) return 0.0;

    const double mmu = kMassesPDG.m_mu;

    const double x = 2.0 * Ee_MeV / mmu;
    const double y = 2.0 * Eg_MeV / mmu;

    // 変数変換: x=2Ee/mmu, y=2Eg/mmu より
    //   dx dy = (2/mmu)^2 dEe dEg  => d^6B/(dEe dEg dΩe dΩγ) = (4/mmu^2) * d^6B/(dx dy dΩe dΩγ)
    const double val = RMD_d6B_dxdy_dOmegae_dOmegag(x, y, cosThetaEG, cosThetaE, cosThetaG, Pmu, d_min);
    return (4.0 / (mmu * mmu)) * val;
}
