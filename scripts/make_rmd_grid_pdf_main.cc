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
// ------------------------------------------------------------

#include <iostream>
#include <cstdlib>

#include "p2meg/MakeRMDGridPdf.h"
#include "p2meg/AnalysisWindow.h"

static bool ParseDouble(const char* s, double& out) {
  if (!s) return false;
  char* end = nullptr;
  const double v = std::strtod(s, &end);
  if (end == s || *end != '\0') return false;
  out = v;
  return true;
}

int main(int argc, char** argv) {
  const char* out = "data/pdf_cache/rmd_grid.root";
  const char* key = "rmd_grid";

  AnalysisWindow4D truth_win = analysis_window;
  if (argc == 5) {
    if (!ParseDouble(argv[1], truth_win.Ee_min) ||
        !ParseDouble(argv[2], truth_win.Ee_max) ||
        !ParseDouble(argv[3], truth_win.Eg_min) ||
        !ParseDouble(argv[4], truth_win.Eg_max)) {
      std::cerr << "Usage: " << argv[0] << " [Ee_min Ee_max Eg_min Eg_max]\n";
      return 1;
    }
  } else if (argc != 1) {
    std::cerr << "Usage: " << argv[0] << " [Ee_min Ee_max Eg_min Eg_max]\n";
    return 1;
  }

  const int rc = MakeRMDGridPdfWithTruthWindow(out, key, truth_win);
  std::cout << "MakeRMDGridPdf returned " << rc << std::endl;
  return rc;
}
