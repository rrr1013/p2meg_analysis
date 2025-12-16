// scripts/make_rmd_mock.cc
// 出力: Ee Egamma t theta（Ee,Egamma: MeV / t: ns / theta: rad）
//
// RMD の擬似データを生成する（解析窓内）
// 手法:
//  - 4D 一様提案 q(Ee,Eg,t,theta) を解析窓内で投げる
//  - ターゲット p = RMDGridPdf(Ee,Eg,t,theta) に対し棄却法でサンプル
//
// 注意:
//  - RMDGridPdf_Load が ROOT ファイルを読むので、ビルドには ROOT が必要
//  - pmax は事前に乱数点で概算して安全係数を掛ける（保守的）

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>
#include <string>

#include "p2meg/AnalysisWindow.h"
#include "p2meg/RMDGridPdf.h"

// 解析窓は include/p2meg/AnalysisWindow.h の analysis_window を使う想定
// 例: inline constexpr AnalysisWindow4D analysis_window{...};

static double Uniform(std::mt19937_64& rng, double a, double b) {
  std::uniform_real_distribution<double> u(a, b);
  return u(rng);
}

int main(int argc, char** argv) {
  // 使い方:
  // ./make_rmd_mock [N] [out] [seed] [grid_root] [key]
  //
  // 例:
  // ./make_rmd_mock 20000 data/mockdata/rmd_mock.dat 12345 data/pdf_cache/rmd_grid.root rmd_grid
  const int N = (argc >= 2) ? std::atoi(argv[1]) : 10000;
  const std::string out = (argc >= 3) ? argv[2] : "-";
  const std::uint64_t seed =
      (argc >= 4) ? static_cast<std::uint64_t>(std::strtoull(argv[3], nullptr, 10))
                  : 12345ULL;
  const std::string grid_root = (argc >= 5) ? argv[4] : "data/pdf_cache/rmd_grid.root";
  const std::string key       = (argc >= 6) ? argv[5] : "rmd_grid";

  if (N <= 0) return 1;

  // ---- PDFロード ----
  if (!RMDGridPdf_Load(grid_root.c_str(), key.c_str())) {
    std::cerr << "RMDGridPdf_Load failed: " << grid_root << " key=" << key << "\n";
    return 1;
  }
  if (!RMDGridPdf_IsLoaded()) {
    std::cerr << "RMDGridPdf_IsLoaded() is false\n";
    return 1;
  }

  // ---- 解析窓 ----
  const AnalysisWindow4D& win = analysis_window;

  // ---- RNG ----
  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> u01(0.0, 1.0);

  // ---- 出力先 ----
  std::ofstream fout;
  std::ostream* os = &std::cout;
  if (out != "-" && !out.empty()) {
    fout.open(out);
    if (!fout) return 1;
    os = &fout;
  }

  // ---- pmax の概算（乱数点で最大値を拾い、安全係数を掛ける）----
  // 解析窓内で均一に点を打って、pmax_est = max(p) を探す
  // 安全係数 safety を掛けて pmax とする（棄却法で p/pmax <= 1 を保証するため）
  const int Nscan = 200000;      // 必要なら増やす
  const double safety = 1.30;    // 保守的係数
  double pmax_est = 0.0;

  for (int i = 0; i < Nscan; ++i) {
    const double Ee = Uniform(rng, win.Ee_min, win.Ee_max);
    const double Eg = Uniform(rng, win.Eg_min, win.Eg_max);
    const double t  = Uniform(rng, win.t_min,  win.t_max);
    const double th = Uniform(rng, win.theta_min, win.theta_max);
    const double p  = RMDGridPdf(Ee, Eg, t, th);
    if (std::isfinite(p) && p > pmax_est) pmax_est = p;
  }

  if (!(pmax_est > 0.0) || !std::isfinite(pmax_est)) {
    std::cerr << "pmax_est is not positive. scan failed?\n";
    return 1;
  }
  const double pmax = safety * pmax_est;

  // ---- ヘッダ ----
  (*os) << "Ee Egamma t theta\n";

  // ---- 棄却法で生成 ----
  int accepted = 0;
  std::uint64_t trials = 0;

  while (accepted < N) {
    ++trials;
    const double Ee = Uniform(rng, win.Ee_min, win.Ee_max);
    const double Eg = Uniform(rng, win.Eg_min, win.Eg_max);
    const double t  = Uniform(rng, win.t_min,  win.t_max);
    const double th = Uniform(rng, win.theta_min, win.theta_max);

    const double p = RMDGridPdf(Ee, Eg, t, th);
    if (!(p > 0.0) || !std::isfinite(p)) continue;

    const double r = u01(rng);
    if (r * pmax < p) {
      (*os) << Ee << " " << Eg << " " << t << " " << th << "\n";
      ++accepted;
    }
  }

  const double acc = static_cast<double>(accepted) / static_cast<double>(trials);
  std::cerr << "[make_rmd_mock] accepted " << accepted << " / trials " << trials
            << " (acceptance=" << acc << ")\n";
  std::cerr << "[make_rmd_mock] pmax_est=" << pmax_est << " pmax=" << pmax
            << " (safety=" << safety << ")\n";

  return 0;
}
