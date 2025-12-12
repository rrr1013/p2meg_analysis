#include "p2meg/MichelSpectrum.h"

#include <cmath>  // 必要なら std::fabs などに使う

// ------------------------------------------------------------
// 補助関数
// ------------------------------------------------------------
double Michel_Emax(const ParticleMasses& ms)
{
    // μ静止系での端点エネルギー（簡略化）
    // 本来 me を含めると端点はわずかに変わるが、ここでは Emax = m_mu/2 を採用する
    return ms.m_mu / 2.0;
}

double Michel_x(double Ee, const ParticleMasses& ms)
{
    // 無次元変数 x = Ee / Emax
    const double Emax = Michel_Emax(ms);
    return Ee / Emax;
}

double Michel_x0(const ParticleMasses& ms)
{
    // η 項に現れる x0 = me / Emax
    const double Emax = Michel_Emax(ms);
    return ms.m_e / Emax;
}

// ------------------------------------------------------------
// ミシェル崩壊スペクトル（非正規化 shape）
// ------------------------------------------------------------
double Michel_d2Shape_dE_dCosTheta(double Ee,
                                   double costh,
                                   double P_mu,
                                   const MichelParams& mp,
                                   const ParticleMasses& ms)
{
    // cosθ の範囲チェック
    if (costh < -1.0 || costh > 1.0) return 0.0;

    // x の範囲チェック（0 < x < 1 の領域のみを有効とする）
    const double x = Michel_x(Ee, ms);
    if (x <= 0.0 || x >= 1.0) return 0.0;

    const double x0 = Michel_x0(ms);

    // よく使われる Michel パラメータ表式（形のみ）：
    //
    // d^2Γ / (x^2 dx dcosθ) ∝
    //   (3 - 3x) + (2/3)ρ(4x-3) + 3η x0 (1-x)/x
    // + P_mu ξ cosθ [ (1-x) + (2/3)δ(4x-3) ]
    //
    // この関数は、全体定数を落として
    //   shape ∝ x^2 * [ ... ]
    // を返す（尤度側で規格化やイベント数を扱う想定）。
    const double term_iso =
        (3.0 - 3.0 * x)
      + (2.0 / 3.0) * mp.rho * (4.0 * x - 3.0)
      + 3.0 * mp.eta * x0 * (1.0 - x) / x;

    const double term_aniso =
        P_mu * mp.xi * costh *
        ((1.0 - x) + (2.0 / 3.0) * mp.delta * (4.0 * x - 3.0));

    const double shape = x * x * (term_iso + term_aniso);

    // パラメータ走査やフィットで負になるケースへの安全策として 0 に丸める
    return (shape > 0.0) ? shape : 0.0;
}

double Michel_dShape_dE(double Ee,
                        const MichelParams& mp,
                        const ParticleMasses& ms)
{
    // 角度積分したエネルギースペクトルの形（非正規化）
    //
    // 角度依存項は cosθ に比例するので、
    // ∫_{-1}^{1} cosθ d(cosθ) = 0 により消える。
    // したがって、角度積分後の形は isotropic 項のみで決まる（定数因子は省略）。
    const double x = Michel_x(Ee, ms);
    if (x <= 0.0 || x >= 1.0) return 0.0;

    const double x0 = Michel_x0(ms);

    const double term_iso =
        (3.0 - 3.0 * x)
      + (2.0 / 3.0) * mp.rho * (4.0 * x - 3.0)
      + 3.0 * mp.eta * x0 * (1.0 - x) / x;

    const double shape = x * x * term_iso;

    // 安全策として負を 0 に丸める
    return (shape > 0.0) ? shape : 0.0;
}
