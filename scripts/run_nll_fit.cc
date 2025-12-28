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

// デバッグ用: 解析窓内イベントの値を表示したいときはコメントを外す
// #define P2MEG_DEBUG_PRINT_WINDOW_EVENTS

// データが解析窓に入っているかの簡易判定
static bool IsInsideAnalysisWindow(const Event& ev, const AnalysisWindow4D& win)
{
  if (!std::isfinite(ev.Ee) || !std::isfinite(ev.Eg) ||
      !std::isfinite(ev.t)  || !std::isfinite(ev.theta)) {
    return false;
  }

  if (ev.Ee < win.Ee_min || ev.Ee > win.Ee_max) return false;
  if (ev.Eg < win.Eg_min || ev.Eg > win.Eg_max) return false;
  if (ev.t  < win.t_min  || ev.t  > win.t_max ) return false;
  if (ev.theta < win.theta_min || ev.theta > win.theta_max) return false;

  return true;
}

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

#ifdef P2MEG_DEBUG_PRINT_WINDOW_EVENTS
    // 窓判定に入ったイベントを、ファイル値(raw)と Event 値の両方で表示
    const bool in_window = IsInsideAnalysisWindow(ev, analysis_window);
    if (in_window) {
      std::cout << "[run_nll_fit][window] raw: Ee=" << cols[0]
                << " Eg=" << cols[1]
                << " t=" << cols[2]
                << " theta=" << cols[3];
      if (cols.size() >= 6) {
        std::cout << " cos_e=" << cols[4] << " cos_g=" << cols[5];
      }
      std::cout << " | event: Ee=" << ev.Ee
                << " Eg=" << ev.Eg
                << " t=" << ev.t
                << " theta=" << ev.theta
                << " cos_e=" << ev.cos_detector_e
                << " cos_g=" << ev.cos_detector_g
                << "\n";
    }
#endif

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

  // 解析窓内のイベント数を数えておく（窓外があるとNLLはペナルティになる）
  int n_in_window = 0;
  for (const auto& ev : events) {
    if (IsInsideAnalysisWindow(ev, analysis_window)) ++n_in_window;
  }

  std::cout << "[run_nll_fit] events in analysis window: "
            << n_in_window << " / " << events.size() << "\n";

  // 解析窓内のみをフィットに使う
  std::vector<Event> events_in_window;
  events_in_window.reserve(events.size());
  for (const auto& ev : events) {
    if (IsInsideAnalysisWindow(ev, analysis_window)) {
      events_in_window.push_back(ev);
    }
  }
  if (events_in_window.empty()) {
    std::cerr << "[run_nll_fit] no events fall inside the analysis window. Aborting fit.\n";
    return 3;
  }

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
  const double N = static_cast<double>(events_in_window.size());
  FitConfig cfg;
  cfg.start_yields = {0.9 * N, 0.1 * N}; // {N_sig, N_rmd}
  cfg.max_calls = 20000;
  cfg.tol = 1e-3;

  // 6) フィット
  const FitResult res = FitNLL(events_in_window, components, cfg);

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
