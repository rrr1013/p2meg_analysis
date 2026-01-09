// ------------------------------------------------------------
// 使い方
//  - RMD の「格子 PDF」を生成するだけの最小 main。
//  - g++ やビルドスクリプトでこのファイルをコンパイルし、生成された実行ファイル
//    （例：`./build/make_rmd_grid_pdf`）を動かすだけでよい。
//  - 実行前に `mkdir -p data/pdf_cache` で出力先ディレクトリを作っておく。
//  - 成功すると `data/pdf_cache/rmd_grid.root` に key=`rmd_grid` の格子 PDF
//    とメタ情報が保存され、解析側は `RMDGridPdf_Load` でこのファイルを読む。
//
// 重要（旧3D版からの変更点）
//  - 生成・保存する格子は 4D（Ee, Eg, i_e, i_g）。
//    i_e, i_g は検出器の離散角（N_theta 分割）に対応するインデックスで、
//      i = 0, 1, 2, ... , N_theta
//    を表す（角度は θ = i * π/N_theta）。
//  - RMDGridPdf 側では Event の (cos_detector_e, cos_detector_g) を
//    N_theta 格子へ最近傍に丸めて (i_e, i_g) を決め、cosΔφ=+1 固定の平面仮定で
//    θ_{eγ} を再構成して解析窓の角度カットに用いる。
//  - 以前の 3D（Ee, Eg, theta）格子 root を残したままだと、ロード側が
//    次元不一致で失敗するので、必ず作り直すこと。
//
// 解析条件を変えて格子を作り直すとき
//  1) `include/p2meg/DetectorResolution.h` や `include/p2meg/AnalysisWindow.h` を編集
//     - 例: Ee/Eg/t の分解能、t_mean、偏極度 P_mu、角度分割数 N_theta、解析窓
//  2) 本ファイル（main）＋格子生成バイナリを再コンパイル
//  3) `./build/make_rmd_grid_pdf` で新しい `rmd_grid.root` を作成
//  4) 続けて解析コード実行（例: `./build/run_nll_fit ...`）
//
// 参考: 一括ビルド用ワンライナー
//  - 前提: カレントが repo 直下、`build/` と `data/pdf_cache/` が存在、root-config が使える
//      mkdir -p build data/pdf_cache
//  - コマンド（make_rmd_grid_pdf と run_nll_fit をまとめてビルド）:
//      g++ -O2 -std=c++17 -Wall -Wextra -pedantic -Iinclude $(root-config --cflags) -o build/make_rmd_grid_pdf macros/make_rmd_grid_pdf_main.cc src/MakeRMDGridPdf.cc src/RMDSpectrum.cc $(root-config --libs) && g++ -O2 -std=c++17 -Wall -Wextra -pedantic -Iinclude $(root-config --cflags) -o build/run_nll_fit scripts/run_nll_fit.cc src/RMDGridPdf.cc src/Likelihood.cc src/NLLFit.cc src/PdfWrappers.cc src/SignalPdf.cc src/ConstraintNLL.cc $(root-config --libs)
// ------------------------------------------------------------

#include <iostream>
#include "p2meg/MakeRMDGridPdf.h"

int main() {
  const char* out = "data/pdf_cache/rmd_grid.root";
  const char* key = "rmd_grid";

  const int rc = MakeRMDGridPdf(out, key);
  std::cout << "MakeRMDGridPdf returned " << rc << std::endl;
  return rc;
}
