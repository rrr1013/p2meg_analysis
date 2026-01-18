// scripts/make_acc_grid_pdf.cc
//
// ============================================================
// p2MEG: ACC 4D 格子 PDF の生成スクリプト
//
// 入力:
//  - 5列のイベント列（Ee Eg t phi_detector_e phi_detector_g）
//  - 先頭/途中にメタ情報行が混ざっていてもよい（5列doubleの行のみ採用）
//
// 出力:
//  - 4D 格子 PDF (Ee, Eg, phi_e, phi_g) を ROOT に保存
//  - 時間因子は解析窓内一様として評価側で掛ける
//
// ビルド例（repo直下で）:
//  g++ -O2 -std=c++17 -Wall -Wextra -pedantic -Iinclude \
//    $(root-config --cflags) \
//    -o build/make_acc_grid_pdf \
//    scripts/make_acc_grid_pdf.cc \
//    src/*.cc \
//    $(root-config --libs)
//
// 実行例:
//  ./build/make_acc_grid_pdf data/mockdata/acc_5000.dat
//  ./build/make_acc_grid_pdf data/mockdata/acc_5000.dat data/pdf_cache/acc_grid.root acc_grid
// ============================================================

#include <cerrno>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "../include/p2meg/Event.h"
#include "../include/p2meg/MakeACCGridPdf.h"

// ------------------------------------------------------------
// 行から double 列を抽出（先頭が数値でない行はスキップ）
// ------------------------------------------------------------
static bool ParseDoublesFromLine(const std::string& line, std::vector<double>& out)
{
  out.clear();

  std::istringstream iss(line);

  std::string first;
  if (!(iss >> first)) return false;

  // コメント行
  if (!first.empty() && first[0] == '#') return false;

  // 1列目が数値でない（ヘッダ行など）ならスキップ
  char* endptr = nullptr;
  const double v0 = std::strtod(first.c_str(), &endptr);
  if (endptr == first.c_str() || *endptr != '\0') return false;

  out.push_back(v0);

  double v = 0.0;
  while (iss >> v) out.push_back(v);

  return true;
}

static bool LoadEventsFromDat(const char* filepath,
                              std::vector<Event>& events,
                              long& n_skipped_non5)
{
  events.clear();
  n_skipped_non5 = 0;

  std::ifstream fin(filepath);
  if (!fin) {
    std::cerr << "[make_acc_grid_pdf] cannot open: " << filepath << "\n";
    return false;
  }

  std::string line;
  std::vector<double> cols;
  while (std::getline(fin, line)) {
    if (!ParseDoublesFromLine(line, cols)) continue;

    // 5列必須：違う行は飛ばす（数も数える）
    if (cols.size() != 5) {
      ++n_skipped_non5;
      continue;
    }

    Event ev{};
    ev.Ee = cols[0];
    ev.Eg = cols[1];
    ev.t  = cols[2];
    ev.phi_detector_e = cols[3];
    ev.phi_detector_g = cols[4];

    events.push_back(ev);
  }

  if (events.empty()) {
    std::cerr << "[make_acc_grid_pdf] no valid 5-column events loaded from: "
              << filepath << "\n";
    return false;
  }

  return true;
}

static void PrintUsage(const char* prog)
{
  std::cerr << "Usage:\n";
  std::cerr << "  " << prog << " <input.dat> [out.root] [key]\n";
}

int main(int argc, char** argv)
{
  if (argc < 2) {
    PrintUsage(argv[0]);
    return 1;
  }

  const char* infile = argv[1];
  const char* outfile = (argc >= 3) ? argv[2] : "data/pdf_cache/acc_grid.root";
  const char* key = (argc >= 4) ? argv[3] : "acc_grid";

  std::vector<Event> events;
  long n_skipped_non5 = 0;

  if (!LoadEventsFromDat(infile, events, n_skipped_non5)) return 1;

  std::cout << "[make_acc_grid_pdf] loaded events=" << events.size()
            << " skipped_non5=" << n_skipped_non5 << "\n";

  std::filesystem::path outpath(outfile);
  if (!outpath.parent_path().empty()) {
    std::error_code ec;
    std::filesystem::create_directories(outpath.parent_path(), ec);
    if (ec) {
      std::cerr << "[make_acc_grid_pdf] cannot create output dir: "
                << outpath.parent_path().string() << "\n";
      return 1;
    }
  }

  const int ret = MakeACCGridPdf(events, outfile, key);
  if (ret != 0) {
    std::cerr << "[make_acc_grid_pdf] MakeACCGridPdf failed (code=" << ret << ")\n";
    return ret;
  }

  return 0;
}
