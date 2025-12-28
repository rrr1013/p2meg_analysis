// ------------------------------------------------------------
// 使い方
//  - RMD の格子 PDF（Ee, Eg, theta）を生成するだけの最小 main。
//  - g++ やビルドスクリプトでこのファイルをコンパイルし、生成された実行ファイル
//    （例：`./build/make_rmd_grid_pdf`）を動かすだけでよい。
//  - 実行前に `mkdir -p data/pdf_cache` で出力先ディレクトリを作っておく。
//  - 成功すると `data/pdf_cache/rmd_grid.root` に key=`rmd_grid` の 3D 格子 PDF
//    とメタ情報が保存され、解析側は `RMDGridPdf_Load` でこのファイルを読む。
//  - 分解能・解析窓を変えて格子を作り直すとき:
//      1) `include/p2meg/DetectorResolution.h` や `include/p2meg/AnalysisWindow.h` を編集
//      2) 本ファイル＋解析側バイナリを再コンパイル（下のワンライナー）
//      3) `./build/make_rmd_grid_pdf` で新しい `rmd_grid.root` を作成
//      4) 続けて NLL などの解析コード実行（例: `./build/run_nll_fit ...`）
//  - 一括ビルド用ワンライナー（カレントが repo 直下、build/data/pdf_cache が既存前提、root-config 必須）:
//      g++ -std=c++17 macros/make_rmd_grid_pdf_main.cc src/MakeRMDGridPdf.cc src/RMDSpectrum.cc -Iinclude $(root-config --cflags --libs) -o build/make_rmd_grid_pdf && g++ -std=c++17 scripts/run_nll_fit.cc src/RMDGridPdf.cc src/MakeRMDGridPdf.cc src/RMDSpectrum.cc src/Likelihood.cc src/NLLFit.cc src/PdfWrappers.cc src/SignalPdf.cc src/ConstraintNLL.cc -Iinclude $(root-config --cflags --libs) -o build/run_nll_fit && g++ -std=c++17 scripts/make_rmd_mock.cc src/RMDGridPdf.cc src/RMDSpectrum.cc src/SignalPdf.cc src/PdfWrappers.cc -Iinclude $(root-config --cflags --libs) -o build/make_rmd_mock && g++ -std=c++17 scripts/make_signal_mock.cc -Iinclude -o build/make_signal_mock && g++ -std=c++17 scripts/test_signal_pdf.cc src/SignalPdf.cc -Iinclude -o build/test_signal_pdf
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
