// scripts/run_nll_fit.cc
//
// 使い方:
//   ./build/run_nll_fit data/mockdata/signal_mock.dat data/pdf_cache/rmd_grid.root rmd_grid
//
// 入力データ形式（空白区切り）:
//   Ee Egamma t theta cos_detector_e cos_detector_g   （6列必須）
//
// 先頭行のヘッダや # コメント行は自動でスキップします。
// 注意:
//  - 本スクリプトは RMD 成分（角度離散化）を使うため、cos_detector_e/g が必須です。
//  - 4列(thetaのみ) から cos_detector_e/g を一意に復元できないため、4列入力は受け付けません。
//  - 読み込んだ theta は、(cos_detector_e, cos_detector_g) を N_theta 格子へ最近傍丸めした後、
//    cosΔφ=+1 固定の平面仮定で再構成した θ_eγ（離散）に置き換えて扱います。

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
 #define P2MEG_DEBUG_PRINT_WINDOW_EVENTS

//============================================================
// 角度（cos_detector_*）→ 離散化 → θ_eγ 再構成
//  - RMDGridPdf.cc と同じロジック（cosΔφ=+1 固定）
//============================================================

static bool IsFinite(double x) { return std::isfinite(x); }

static double Clamp(double x, double lo, double hi)
{
  if (x < lo) return lo;
  if (x > hi) return hi;
  return x;
}

static int ThetaNearestIndex(double theta, int N_theta)
{
  if (!(N_theta >= 1)) return -1;
  if (!IsFinite(theta)) return -1;
  if (theta < 0.0 || theta > pi) return -1;

  const double step = pi / static_cast<double>(N_theta);
  if (!(step > 0.0) || !IsFinite(step)) return -1;

  long long i_ll = std::llround(theta / step);
  if (i_ll < 0LL) i_ll = 0LL;
  if (i_ll > static_cast<long long>(N_theta)) i_ll = static_cast<long long>(N_theta);
  return static_cast<int>(i_ll);
}

static int CosToIndex(double cosv, int N_theta)
{
  if (!(N_theta >= 1)) return -1;
  if (!IsFinite(cosv)) return -1;

  const double c = Clamp(cosv, -1.0, 1.0);
  const double theta = std::acos(c);
  if (!IsFinite(theta)) return -1;

  return ThetaNearestIndex(theta, N_theta);
}

static bool ThetaEG_FromIndices(int ie, int ig, int N_theta, double& thetaEG_out)
{
  if (!(N_theta >= 1)) return false;
  if (ie < 0 || ie > N_theta) return false;
  if (ig < 0 || ig > N_theta) return false;

  const double thetaE = pi * static_cast<double>(ie) / static_cast<double>(N_theta);
  const double thetaG = pi * static_cast<double>(ig) / static_cast<double>(N_theta);

  // cosΔφ=+1（Δφ=0）固定 → cos(thetaE - thetaG)
  const double cEG = std::cos(thetaE - thetaG);
  const double c = Clamp(cEG, -1.0, 1.0);

  const double th = std::acos(c);
  if (!IsFinite(th)) return false;

  thetaEG_out = th;
  return true;
}

static bool ReconstructThetaEG_Discrete(double cos_detector_e, double cos_detector_g,
                                       int N_theta, double& thetaEG_out)
{
  const int ie = CosToIndex(cos_detector_e, N_theta);
  const int ig = CosToIndex(cos_detector_g, N_theta);
  if (ie < 0 || ig < 0) return false;

  return ThetaEG_FromIndices(ie, ig, N_theta, thetaEG_out);
}

//============================================================
// データが解析窓に入っているかの簡易判定（Event.theta は再構成後を想定）
//============================================================

static bool IsInsideAnalysisWindow(const Event& ev, const AnalysisWindow4D& win)
{
  if (!IsFinite(ev.Ee) || !IsFinite(ev.Eg) ||
      !IsFinite(ev.t)  || !IsFinite(ev.theta)) {
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

  const int N_theta = detres.N_theta;
  if (!(N_theta >= 1)) {
    std::cerr << "[run_nll_fit] detres.N_theta must be >= 1\n";
    return false;
  }

  std::string line;
  std::vector<double> cols;

  long line_no = 0;
  long skipped_bad_cols = 0;
  while (std::getline(fin, line)) {
    ++line_no;
    if (!ParseDoublesFromLine(line, cols)) continue;

    // 必須: Ee Eg t theta cos_e cos_g（6列）
    if (cols.size() < 6) {
      ++skipped_bad_cols;
      std::cerr << "[run_nll_fit] WARN: need 6 columns (Ee Eg t theta cos_detector_e cos_detector_g). "
                << "skipping line=" << line_no << " cols=" << cols.size() << "\n";
      continue;
    }

    Event ev{};
    ev.Ee    = cols[0];
    ev.Eg    = cols[1];
    ev.t     = cols[2];

    // 入力theta（raw）は保持しない（後で再構成thetaに置き換える）
    const double theta_raw = cols[3];

    ev.cos_detector_e = cols[4];
    ev.cos_detector_g = cols[5];

    // cos_detector_e/g から、N_theta 最近傍格子で離散化した θ_eγ を再構成して ev.theta に入れる
    double theta_reco = 0.0;
    if (!ReconstructThetaEG_Discrete(ev.cos_detector_e, ev.cos_detector_g, N_theta, theta_reco)) {
      std::cerr << "[run_nll_fit] ERROR: failed to reconstruct theta from cos detectors. "
                << "line=" << line_no
                << " cos_e=" << ev.cos_detector_e
                << " cos_g=" << ev.cos_detector_g << "\n";
      return false;
    }
    ev.theta = theta_reco;

#ifdef P2MEG_DEBUG_PRINT_WINDOW_EVENTS
    const bool in_window = IsInsideAnalysisWindow(ev, analysis_window);
    if (in_window) {
      std::cout << "[run_nll_fit][window] raw: Ee=" << cols[0]
                << " Eg=" << cols[1]
                << " t=" << cols[2]
                << " theta_raw=" << theta_raw
                << " cos_e=" << cols[4]
                << " cos_g=" << cols[5]
                << " | event: Ee=" << ev.Ee
                << " Eg=" << ev.Eg
                << " t=" << ev.t
                << " theta(reco)=" << ev.theta
                << " cos_e=" << ev.cos_detector_e
                << " cos_g=" << ev.cos_detector_g
                << "\n";
    }
#else
    (void)theta_raw;
#endif

    events.push_back(ev);

    if (max_events > 0 && static_cast<int>(events.size()) >= max_events) break;
  }

  if (events.empty()) {
    std::cerr << "[run_nll_fit] no events loaded from: " << filepath << "\n";
    return false;
  }
  if (skipped_bad_cols > 0) {
    std::cerr << "[run_nll_fit] WARN: skipped " << skipped_bad_cols
              << " line(s) with <6 columns in: " << filepath << "\n";
  }
  return true;
}

int main(int argc, char** argv)
{
  const char* datafile = (argc >= 2) ? argv[1] : "data/mockdata/signal_mock.dat";
  const char* rmd_root = (argc >= 3) ? argv[2] : "data/pdf_cache/rmd_grid.root";
  const char* rmd_key  = (argc >= 4) ? argv[3] : "rmd_grid";

  // 1) 入力データ読み込み（theta は再構成した離散 θ_eγ に置き換えられる）
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
