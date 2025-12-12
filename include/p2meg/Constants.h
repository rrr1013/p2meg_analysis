#ifndef P2MEG_CONSTANTS_H
#define P2MEG_CONSTANTS_H

// ============================================================
// p2MEG 解析用の定数・パラメータ定義
//
// ・単位は MeV（エネルギー・質量）
// ・ヘッダオンリーで使うため、inline constexpr で多重定義を防ぐ
// ============================================================

// ミシェルパラメータ（相互作用の形だけを表す）
struct MichelParams {
    double rho;
    double eta;
    double xi;
    double delta;
};

// 粒子質量（運動学の入力）
struct ParticleMasses {
    double m_mu; // [MeV]
    double m_e;  // [MeV]
};

// ---- 既定値：標準模型 (V-A) のミシェルパラメータ ----
inline constexpr MichelParams kMichelSM{
    0.75, // rho
    0.0,  // eta
    1.0,  // xi
    0.75  // delta
};

// ---- 既定値：質量（“PDG値相当”） ----
inline constexpr ParticleMasses kMassesPDG{
    105.658, // m_mu [MeV]
    0.511    // m_e  [MeV]
};

// 円周率
inline constexpr double pi =
    3.141592653589793238462643383279502884;

// 微細構造定数
inline constexpr double alpha =
    1.0 / 137.035999084;

#endif // P2MEG_CONSTANTS_H
