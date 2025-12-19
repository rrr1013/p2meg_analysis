// scripts/run_nll_fit.cc
//
// 使い方:
//   ./build/run_nll_fit data/mockdata/signal_mock.dat data/pdf_cache/rmd_grid.root rmd_grid
//
// 入力データ形式（空白区切り）:
//   Ee Egamma t theta
//   Ee Egamma t theta cos_detector_e cos_detector_g   （6列でも可）
//
// 先頭行のヘッダや # コメント行は自動でスキップします。

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cmath>

#include "p2meg/Event.h"
#include "p2meg/Likelihood.h"
#include "p2meg/NLLFit.h"
#include "p2meg/PdfWrappers.h"

#include "p2meg/AnalysisWindow.h"
#include "p2meg/DetectorResolution.h"
#include "p2meg/Constants.h"
#include "p2meg/RMDGridPdf.h"

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

static bool LoadEventsFromDat(const char* filepath, std::vector<Event>& events, int max_events = -1)
{
  events.clear();

  std::ifstream fin(filepath);
  if (!fin) {
    std::cerr << "[run_nll_fit] cannot open: " << filepath << "\n";
    return false;
  }

  std::string line;
  std::vector<double> cols;
  while (std::getline(fin, line)) {
    if (!ParseDoublesFromLine(line, cols)) continue;
    if (cols.size() < 4) continue;

    Event ev{};
    ev.Ee    = cols[0];
    ev.Eg    = cols[1];
    ev.t     = cols[2];
    ev.theta = cols[3];

    // 将来用（なければ 0）
    ev.cos_detector_e = (cols.size() >= 6) ? cols[4] : 0.0;
    ev.cos_detector_g = (cols.size() >= 6) ? cols[5] : 0.0;

    events.push_back(ev);

    if (max_events > 0 && static_cast<int>(events.size()) >= max_events) break;
  }

  if (events.empty()) {
    std::cerr << "[run_nll_fit] no events loaded from: " << filepath << "\n";
    return false;
  }
  return true;
}

int main(int argc, char** argv)
{
  const char* datafile = (argc >= 2) ? argv[1] : "data/mockdata/signal_mock.dat";
  const char* rmd_root = (argc >= 3) ? argv[2] : "data/pdf_cache/rmd_grid.root";
  const char* rmd_key  = (argc >= 4) ? argv[3] : "rmd_grid";

  // 1) 入力データ読み込み
  std::vector<Event> events;
  if (!LoadEventsFromDat(datafile, events)) return 1;

  std::cout << "[run_nll_fit] loaded events: " << events.size()
            << " from " << datafile << "\n";

  // 2) RMD 格子PDFロード
  if (!RMDGridPdf_Load(rmd_root, rmd_key)) {
    std::cerr << "[run_nll_fit] RMDGridPdf_Load failed: "
              << rmd_root << " key=" << rmd_key << "\n";
    return 2;
  }

  // 3) Signal コンテキスト（既定の解析窓・分解能・質量）
  static SignalPdfContext sigctx{
    analysis_window,
    detres,
    kMassesPDG
  };

  // 4) 成分（sig, rmd）
  std::vector<PdfComponent> components;
  components.push_back(MakeSignalComponent(&sigctx));
  components.push_back(MakeRMDComponent());

  // 5) 初期値（データを差し替えても動くようにNに比例させる）
  const double N = static_cast<double>(events.size());
  FitConfig cfg;
  cfg.start_yields = {0.9 * N, 0.1 * N}; // {N_sig, N_rmd}
  cfg.max_calls = 20000;
  cfg.tol = 1e-3;

  // 6) フィット
  const FitResult res = FitNLL(events, components, cfg);

  std::cout << "==================== Fit Result ====================\n";
  std::cout << "status   = " << res.status << "\n";
  std::cout << "nll_min  = " << res.nll_min << "\n";
  if (res.yields_hat.size() >= 2) {
    std::cout << "N_sig_hat = " << res.yields_hat[0] << "\n";
    std::cout << "N_rmd_hat = " << res.yields_hat[1] << "\n";
  }
  if (res.yields_err.size() >= 2) {
    std::cout << "err_sig   = " << res.yields_err[0] << "\n";
    std::cout << "err_rmd   = " << res.yields_err[1] << "\n";
  }
  std::cout << "====================================================\n";

  return 0;
}
