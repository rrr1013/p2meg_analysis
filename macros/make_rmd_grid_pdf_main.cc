// ------------------------------------------------------------
// 使い方
//  - RMD の「格子 PDF」を生成するだけの最小 main。
//  - g++ やビルドスクリプトでこのファイルをコンパイルし、生成された実行ファイル
//    （例：`./build/make_rmd_grid_pdf`）を動かすだけでよい。
//  - 実行前に `mkdir -p data/pdf_cache` で出力先ディレクトリを作っておく。
//  - 成功すると `data/pdf_cache/rmd_grid.root` に key=`rmd_grid` の格子 PDF
//    とメタ情報が保存され、解析側は `RMDGridPdf_Load` でこのファイルを読む。
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
