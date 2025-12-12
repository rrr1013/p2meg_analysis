// scripts/make_signal_mock.cc
// 出力: Ee Egamma t theta（Ee,Egamma: MeV / t: ns / theta: rad）

// signalの擬似データを生成する

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>
#include <string>

#include "p2meg/Constants.h"

static double Clamp(double x, double lo, double hi) {
  if (x < lo) return lo;
  if (x > hi) return hi;
  return x;
}

int main(int argc, char** argv) {
  // 使い方: ./make_signal_mock [N] [out] [seed]
  const int N = (argc >= 2) ? std::atoi(argv[1]) : 10000;
  const std::string out = (argc >= 3) ? argv[2] : "-";
  const std::uint64_t seed = (argc >= 4) ? static_cast<std::uint64_t>(std::strtoull(argv[3], nullptr, 10)) : 12345ULL;
  if (N <= 0) return 1;

  // ---- 真値（停止μの2体崩壊：運動学）----
  const double m_mu = kMassesPDG.m_mu;  // [MeV]
  const double m_e  = kMassesPDG.m_e;   // [MeV]
  const double Ee0 = (m_mu*m_mu + m_e*m_e) / (2.0*m_mu);  // [MeV]
  const double Eg0 = (m_mu*m_mu - m_e*m_e) / (2.0*m_mu);  // [MeV]
  const double t0  = 0.0;                                   // [ns]
  const double th0 = pi;                                    // [rad] back-to-back

  // ---- 分解能（標準偏差）----
  const double sigmaEe = 0.10 * Ee0;              // [MeV]
  const double sigmaEg = 0.10 * Eg0;              // [MeV]
  const double sigmat  = 1.0;                     // [ns]
  const double sigmath = 10.0 * pi / 180.0;       // [rad] = 10°

  std::mt19937_64 rng(seed);
  std::normal_distribution<double> gEe(Ee0, sigmaEe);
  std::normal_distribution<double> gEg(Eg0, sigmaEg);
  std::normal_distribution<double> gt (t0,  sigmat);
  std::normal_distribution<double> gdt(0.0, sigmath);  // 角度の“ずれ” δθ

  std::ofstream fout;
  std::ostream* os = &std::cout;
  if (out != "-" && !out.empty()) {
    fout.open(out);
    if (!fout) return 1;
    os = &fout;
  }

  (*os) << "Ee Egamma t theta\n";
  for (int i = 0; i < N; ++i) {
    const double Ee = std::max(0.0, gEe(rng));  // 再構成ガード（物理カットではない）
    const double Eg = std::max(0.0, gEg(rng));  // 再構成ガード（物理カットではない）
    const double t  = gt(rng);

    // θ は acos(dot) 由来で本来 0..pi：ここでは π 近傍のずれを |δθ| で表して折り返す
    const double th = Clamp(th0 - std::fabs(gdt(rng)), 0.0, pi);

    (*os) << Ee << " " << Eg << " " << t << " " << th << "\n";
  }
  return 0;
}
