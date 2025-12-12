#include <iostream>
#include <iomanip>
#include <cmath>

#include "p2meg/RMDSpectrum.h"
#include "p2meg/Constants.h"

//============================================================
// RMDSpectrum の簡易テスト（更新版）
//
// 目的：
//  1) コンパイル・リンクが通るか
//  2) q^2>=0 の物理領域チェックが効いて、領域外が 0 になるか
//  3) x,y 版と Ee,Eg 版が変数変換で一致するか
//============================================================

static double r_mass_ratio2()
{
    const double mmu = kMassesPDG.m_mu;
    const double me  = kMassesPDG.m_e;
    return (me * me) / (mmu * mmu);
}

static double beta_from_x(double x)
{
    const double r = r_mass_ratio2();
    const double t = 1.0 - 4.0 * r / (x * x);
    return (t > 0.0) ? std::sqrt(t) : 0.0;
}

static bool is_allowed_q2(double x, double y, double cosTheta)
{
    const double r = r_mass_ratio2();
    const double beta = beta_from_x(x);
    if (!(beta > 0.0)) return false;

    const double d = 1.0 - beta * cosTheta;
    const double q2_over_m2 = 1.0 + r - x - y + 0.5 * x * y * d;
    return (q2_over_m2 >= -1e-12);
}

static void PrintOnePoint(double Ee, double Eg, double c)
{
    const double mmu = kMassesPDG.m_mu;

    const double x = 2.0 * Ee / mmu;
    const double y = 2.0 * Eg / mmu;

    const bool ok = is_allowed_q2(x, y, c);

    const double v_xy = RMD_d3B_dxdy_dcos(x, y, c);
    const double v_E  = RMD_d3B_dEe_dEg_dcos(Ee, Eg, c);

    const double v_E_from_xy = (4.0 / (mmu * mmu)) * v_xy;

    std::cout << std::setprecision(10);
    std::cout << "Ee=" << Ee << " MeV, Eg=" << Eg << " MeV, cos=" << c
              << "  allowed(q^2>=0)=" << (ok ? "yes" : "no") << "\n";
    std::cout << "  x=" << x << ", y=" << y << "\n";
    std::cout << "  d3B/dx dy dcos = " << v_xy << "\n";
    std::cout << "  d3B/dEe dEg dcos (direct)     = " << v_E << "\n";
    std::cout << "  d3B/dEe dEg dcos (from x,y)   = " << v_E_from_xy << "\n";
    std::cout << "  diff (direct - from)          = " << (v_E - v_E_from_xy) << "\n\n";
}

int main()
{
    std::cout << "==== RMDSpectrum simple test (with q^2 check) ====\n";
    std::cout << "alpha=" << alpha << ", pi=" << pi << "\n";
    std::cout << "m_mu=" << kMassesPDG.m_mu << " MeV, m_e=" << kMassesPDG.m_e << " MeV\n\n";

    // 代表点（物理的に可能な点を選ぶ）
    // Ee=50, Eg=40 は角度がほぼ反平行でないと q^2<0 になりやすいので cos=-0.99 にする
    PrintOnePoint(50.0, 40.0, -0.99);

    // 角度スキャン：Ee=50, Eg=40（物理領域が狭い）
    std::cout << "---- scan cos(theta_eγ)  (Ee=50, Eg=40) ----\n";
    for (int i = 0; i <= 20; ++i) {
        const double c = -1.0 + 0.1 * i;
        const double x = 2.0 * 50.0 / kMassesPDG.m_mu;
        const double y = 2.0 * 40.0 / kMassesPDG.m_mu;
        const bool ok = is_allowed_q2(x, y, c);

        const double v = RMD_d3B_dEe_dEg_dcos(50.0, 40.0, c);
        std::cout << "cos=" << std::setw(5) << c
                  << "  allowed=" << (ok ? "yes" : " no")
                  << "  d3B=" << std::setprecision(10) << v << "\n";
    }
    std::cout << "\n";

    // もう少し広い領域が出る例：Eg を下げる
    std::cout << "---- scan cos(theta_eγ)  (Ee=50, Eg=10) ----\n";
    for (int i = 0; i <= 20; ++i) {
        const double c = -1.0 + 0.1 * i;
        const double x = 2.0 * 50.0 / kMassesPDG.m_mu;
        const double y = 2.0 * 10.0 / kMassesPDG.m_mu;
        const bool ok = is_allowed_q2(x, y, c);

        const double v = RMD_d3B_dEe_dEg_dcos(50.0, 10.0, c);
        std::cout << "cos=" << std::setw(5) << c
                  << "  allowed=" << (ok ? "yes" : " no")
                  << "  d3B=" << std::setprecision(10) << v << "\n";
    }
    std::cout << "\n";

    std::cout << "==== done ====\n";
    return 0;
}
