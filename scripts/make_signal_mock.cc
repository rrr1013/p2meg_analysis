// scripts/make_signal_mock.cc
// 出力: Ee Egamma t theta（Ee,Egamma: MeV / t: ns / theta: rad）
//
// signal（μ+→e+γ）の擬似データを生成する
// 角度は δ = pi - theta (δ>=0) を「ガウシアン→折り返し（half-normal）→thetaに戻す」
// さらに解析窓 [theta_min, theta_max] に対応する δ 範囲でトランケート（窓内）して生成する。

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>
#include <string>

#include "p2meg/Constants.h"
#include "p2meg/AnalysisWindow.h"
#include "p2meg/DetectorResolution.h"

int main(int argc, char** argv) {
  // 使い方: ./make_signal_mock [N] [out] [seed]
  const int N = (argc >= 2) ? std::atoi(argv[1]) : 10000;
  const std::string out = (argc >= 3) ? argv[2] : "-";
  const std::uint64_t seed =
      (argc >= 4) ? static_cast<std::uint64_t>(std::strtoull(argv[3], nullptr, 10))
                  : 12345ULL;
  if (N <= 0) return 1;

  // ---- 解析窓・分解能（プロジェクトで定義している定数を使用）----
  const AnalysisWindow4D& win = analysis_window;
  const DetectorResolutionConst& res = detres;

  // ---- 真値（停止μの2体崩壊：運動学）----
  const double m_mu = kMassesPDG.m_mu;  // [MeV]
  const double m_e  = kMassesPDG.m_e;   // [MeV]
  const double Ee0 = (m_mu * m_mu + m_e * m_e) / (2.0 * m_mu);  // [MeV]
  const double Eg0 = (m_mu * m_mu - m_e * m_e) / (2.0 * m_mu);  // [MeV]
  const double t0  = res.t_mean;                                // [ns]
  const double th0 = pi;                                        // [rad] back-to-back

  // ---- 分解能（標準偏差）----
  const double sigmaEe = res.sigma_Ee;      // [MeV]
  const double sigmaEg = res.sigma_Eg;      // [MeV]
  const double sigmat  = res.sigma_t;       // [ns]
  const double sigmath = res.sigma_theta;   // [rad]（delta の幅）

  if (!(sigmaEe > 0.0) || !(sigmaEg > 0.0) || !(sigmat > 0.0) || !(sigmath > 0.0)) return 1;

  // ---- 角度のトランケート条件（theta窓 <-> delta窓）----
  // theta in [win.theta_min, win.theta_max]
  // delta = pi - theta なので、delta in [pi-theta_max, pi-theta_min]
  const double dmin = pi - win.theta_max;
  const double dmax = pi - win.theta_min;

  // half-normal は delta>=0 のみなので、下限は 0 にクリップ
  const double dlo = (dmin < 0.0) ? 0.0 : dmin;
  const double dhi = dmax;

  if (!(dlo < dhi)) return 1;

  std::mt19937_64 rng(seed);
  std::normal_distribution<double> gEe(Ee0, sigmaEe);
  std::normal_distribution<double> gEg(Eg0, sigmaEg);
  std::normal_distribution<double> gt (t0,  sigmat);
  std::normal_distribution<double> g0 (0.0, sigmath);  // delta の元になるガウシアン

  std::ofstream fout;
  std::ostream* os = &std::cout;
  if (out != "-" && !out.empty()) {
    fout.open(out);
    if (!fout) return 1;
    os = &fout;
  }

  (*os) << "Ee Egamma t theta\n";

  for (int i = 0; i < N; ++i) {
    // Ee, Eg, t は単純にガウシアン（解析窓で切るなら下でチェック）
    const double Ee = gEe(rng);
    const double Eg = gEg(rng);
    const double t  = gt(rng);

    // 角度: δ = |N(0, sigmath)| を生成し、解析窓に対応する δ 範囲でトランケート
    double delta = 0.0;
    while (true) {
      delta = std::fabs(g0(rng));  // 折り返し（half-normal）
      if (delta >= dlo && delta <= dhi) break;
    }
    const double theta = th0 - delta; // theta = pi - delta（delta>=0 なので theta<=pi）

    // 解析窓外は捨てて「窓内サンプル」を作る（窓内正規化PDFと整合）
    if (Ee < win.Ee_min || Ee > win.Ee_max) { --i; continue; }
    if (Eg < win.Eg_min || Eg > win.Eg_max) { --i; continue; }
    if (t  < win.t_min  || t  > win.t_max ) { --i; continue; }
    if (theta < win.theta_min || theta > win.theta_max) { --i; continue; }

    (*os) << Ee << " " << Eg << " " << t << " " << theta << "\n";
  }

  return 0;
}
