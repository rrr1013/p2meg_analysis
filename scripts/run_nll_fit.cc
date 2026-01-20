// scripts/run_nll_fit.cc
//
// 使い方:
//  g++ -O2 -std=c++17 -Wall -Wextra -pedantic -Iinclude $(root-config --cflags) -o build/run_nll_fit scripts/run_nll_fit.cc src/*.cc $(root-config --libs)
//   ./build/run_nll_fit data/mockdata/testdata1.dat
//   Ee Egamma t phi_detector_e phi_detector_g
//
// 先頭行のヘッダや # コメント行は自動でスキップします。
//
//  - phi_detector_e/g は 0..pi の範囲を想定する（解析側でクランプされる）。
//  - 5列でない行はスキップし、スキップ数を表示します。

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <limits>

#include "p2meg/Event.h"
#include "p2meg/Likelihood.h"
#include "p2meg/NLLFit.h"
#include "p2meg/PdfWrappers.h"

#include "p2meg/AnalysisWindow.h"
#include "p2meg/DetectorResolution.h"
#include "p2meg/Constants.h"
#include "p2meg/RMDGridPdf.h"
#include "p2meg/ACCGridPdf.h"

// デバッグ用: 解析窓内イベントの値を表示したいときはコメントを外す
// #define P2MEG_DEBUG_PRINT_WINDOW_EVENTS

// デバッグ用: フィット失敗時に pi<=0 を起こすイベントを表示したいときはコメントを外す
#define P2MEG_DEBUG_PRINT_ZERO_PI_EVENTS

// デバッグ用: pi<=0 イベントの q^2 を表示したいときはコメントを外す
#define P2MEG_DEBUG_PRINT_ZERO_PI_Q2

// 固定設定（必要ならコードを書き換える方針）
static const char* kDefaultRmdRoot = "data/pdf_cache/rmd_grid.root";
static const char* kDefaultRmdKey  = "rmd_grid";
static const char* kDefaultAccRoot = "data/pdf_cache/acc_grid.root";
static const char* kDefaultAccKey  = "acc_grid";

static bool IsFinite(double x) { return std::isfinite(x); }

static double Clamp(double x, double lo, double hi)
{
  if (x < lo) return lo;
  if (x > hi) return hi;
  return x;
}

static double ThetaFromPhi(double phi_e, double phi_g)
{
  if (!IsFinite(phi_e) || !IsFinite(phi_g)) return 0.0;
  const double pe = Clamp(phi_e, 0.0, pi);
  const double pg = Clamp(phi_g, 0.0, pi);
  return std::fabs(pe - pg);
}

// RMD の q^2/m_mu^2 を評価（参考: RMDSpectrum の定義）
static double RmdQ2OverM2(double Ee, double Eg, double phi_e, double phi_g)
{
  if (!IsFinite(Ee) || !IsFinite(Eg)) return std::numeric_limits<double>::quiet_NaN();
  if (!(Ee > 0.0) || !(Eg > 0.0)) return std::numeric_limits<double>::quiet_NaN();

  const double mmu = kMassesPDG.m_mu; // [MeV]
  const double me  = kMassesPDG.m_e;  // [MeV]
  if (!(mmu > 0.0) || !IsFinite(mmu) || !IsFinite(me)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  // x=2Ee/m_mu, y=2Eg/m_mu（無次元）
  const double x = 2.0 * Ee / mmu;
  const double y = 2.0 * Eg / mmu;
  if (!(x > 0.0) || !(y > 0.0)) return std::numeric_limits<double>::quiet_NaN();

  // r=(me/m_mu)^2（無次元）
  const double r = (me * me) / (mmu * mmu);

  // beta = sqrt(1 - 4r/x^2)（無次元）
  const double tt = 1.0 - 4.0 * r / (x * x);
  if (!(tt > 0.0) || !IsFinite(tt)) return std::numeric_limits<double>::quiet_NaN();
  const double beta = std::sqrt(tt);

  // theta_eg = |phi_e - phi_g|（rad）
  const double theta_eg = ThetaFromPhi(phi_e, phi_g);
  if (!(theta_eg >= 0.0)) return std::numeric_limits<double>::quiet_NaN();

  // d = 1 - beta * cos(theta_eg)（無次元）
  const double d = 1.0 - beta * std::cos(theta_eg);

  // q^2/m_mu^2 = 1 + r - x - y + (x*y/2)*d（無次元）
  const double q2_over_m2 = 1.0 + r - x - y + 0.5 * x * y * d;
  return IsFinite(q2_over_m2) ? q2_over_m2 : std::numeric_limits<double>::quiet_NaN();
}

// データが解析窓に入っているかの簡易判定
static bool IsInsideAnalysisWindow(const Event& ev, const AnalysisWindow4D& win)
{
  if (!IsFinite(ev.Ee) || !IsFinite(ev.Eg) ||
      !IsFinite(ev.t)  || !IsFinite(ev.phi_detector_e) || !IsFinite(ev.phi_detector_g)) {
    return false;
  }

  const double theta_eg = ThetaFromPhi(ev.phi_detector_e, ev.phi_detector_g);

  if (ev.Ee < win.Ee_min || ev.Ee > win.Ee_max) return false;
  if (ev.Eg < win.Eg_min || ev.Eg > win.Eg_max) return false;
  if (ev.t  < win.t_min  || ev.t  > win.t_max ) return false;
  if (theta_eg < win.theta_min || theta_eg > win.theta_max) return false;

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

static double EvalPiWithStartYields(const Event& ev,
                                    const std::vector<PdfComponent>& components,
                                    const std::vector<double>& yields)
{
  const std::size_t npar = components.size();
  double pi_sum = 0.0;
  for (std::size_t k = 0; k < npar; ++k) {
    double pk = components[k].eval ? components[k].eval(ev, components[k].ctx) : 0.0;
    if (!std::isfinite(pk) || pk < 0.0) pk = 0.0;
    const double Nk = (k < yields.size()) ? yields[k] : 0.0;
    pi_sum += Nk * pk;
  }
  return pi_sum;
}

static bool LoadEventsFromDat(const char* filepath,
                             std::vector<Event>& events,
                             long& n_skipped_non5,
                             int max_events = -1)
{
  events.clear();
  n_skipped_non5 = 0;

  std::ifstream fin(filepath);
  if (!fin) {
    std::cerr << "[run_nll_fit] cannot open: " << filepath << "\n";
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

    // 生値
    const double Ee_raw    = cols[0];
    const double Eg_raw    = cols[1];
    const double t_raw     = cols[2];
    const double phi_e_raw = cols[3];
    const double phi_g_raw = cols[4];

    Event ev{};
    ev.Ee = Ee_raw;
    ev.Eg = Eg_raw;
    ev.t  = t_raw;
    ev.phi_detector_e = phi_e_raw;
    ev.phi_detector_g = phi_g_raw;

#ifdef P2MEG_DEBUG_PRINT_WINDOW_EVENTS
    const bool in_window = IsInsideAnalysisWindow(ev, analysis_window);
    if (in_window) {
      const double theta_eg = ThetaFromPhi(ev.phi_detector_e, ev.phi_detector_g);
      std::cout << "[run_nll_fit][window] raw: Ee=" << Ee_raw
                << " Eg=" << Eg_raw
                << " t=" << t_raw
                << " phi_e=" << phi_e_raw
                << " phi_g=" << phi_g_raw
                << " | event: Ee=" << ev.Ee
                << " Eg=" << ev.Eg
                << " t=" << ev.t
                << " theta_eg=" << theta_eg
                << " phi_e=" << ev.phi_detector_e
                << " phi_g=" << ev.phi_detector_g
                << "\n";
    }
#endif

    events.push_back(ev);

    if (max_events > 0 && static_cast<int>(events.size()) >= max_events) break;
  }

  if (events.empty()) {
    std::cerr << "[run_nll_fit] no valid 5-column events loaded from: " << filepath << "\n";
    return false;
  }
  return true;
}

int main(int argc, char** argv)
{
  const char* datafile = (argc >= 2) ? argv[1] : "data/mockdata/testdata1.dat";

  // 1) 入力データ読み込み（5列必須）
  std::vector<Event> events;
  long n_skipped_non5 = 0;
  if (!LoadEventsFromDat(datafile, events, n_skipped_non5)) return 1;

  std::cout << "[run_nll_fit] loaded events: " << events.size()
            << " from " << datafile << "\n";
  std::cout << "[run_nll_fit] skipped non-5-column lines: " << n_skipped_non5 << "\n";

  // 解析窓内のイベント数を数えておく
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

  // 2) RMD 格子PDFロード（固定設定）
  if (!RMDGridPdf_Load(kDefaultRmdRoot, kDefaultRmdKey)) {
    std::cerr << "[run_nll_fit] RMDGridPdf_Load failed: "
              << kDefaultRmdRoot << " key=" << kDefaultRmdKey << "\n";
    return 2;
  }

  // 3) ACC 格子PDFロード（固定設定）
  if (!ACCGridPdf_Load(kDefaultAccRoot, kDefaultAccKey)) {
    std::cerr << "[run_nll_fit] ACCGridPdf_Load failed: "
              << kDefaultAccRoot << " key=" << kDefaultAccKey << "\n";
    return 2;
  }

  // 4) Signal コンテキスト（既定の解析窓・分解能・質量）
  static SignalPdfContext sigctx{
    analysis_window,
    detres,
    kMassesPDG
  };

  // 5) 成分（sig, rmd, acc）
  std::vector<PdfComponent> components;
  components.push_back(MakeSignalComponent(&sigctx));
  components.push_back(MakeRMDComponent());
  components.push_back(MakeACCComponent());

  // 6) 初期値（均等割り）
  const double N0 = static_cast<double>(events_in_window.size());
  FitConfig cfg;
  cfg.start_yields = {N0 / 3.0, N0 / 3.0, N0 / 3.0}; // {N_sig, N_rmd, N_acc}
  cfg.max_calls = 20000;
  cfg.tol = 1e-3;

  // 7) pi<=0 のイベントはフィットから除外（デバッグ表示は後段のまま）
  std::vector<Event> events_fit;
  events_fit.reserve(events_in_window.size());
  for (const auto& ev : events_in_window) {
    const double pi_sum = EvalPiWithStartYields(ev, components, cfg.start_yields);
    if (pi_sum > 0.0 && std::isfinite(pi_sum)) {
      events_fit.push_back(ev);
    }
  }
  if (events_fit.empty()) {
    std::cerr << "[run_nll_fit] no events remain after pi<=0 cut. Aborting fit.\n";
    return 4;
  }

  const double N = static_cast<double>(events_fit.size());
  cfg.start_yields = {N / 3.0, N / 3.0, N / 3.0}; // {N_sig, N_rmd, N_acc}

  // 8) フィット
  const FitResult res = FitNLL(events_fit, components, cfg);

#ifdef P2MEG_DEBUG_PRINT_ZERO_PI_EVENTS
    for (double y : cfg.start_yields) std::cout << " " << y;
    std::cout << "\n";

    int n_zero_pi = 0;
    const std::size_t npar = components.size();
    std::vector<double> pks(npar, 0.0);

    for (const auto& ev : events_in_window) {
      double pi_sum = 0.0;
      for (std::size_t k = 0; k < npar; ++k) {
        double pk = components[k].eval ? components[k].eval(ev, components[k].ctx) : 0.0;
        if (!std::isfinite(pk) || pk < 0.0) pk = 0.0;
        pks[k] = pk;
        const double Nk = (k < cfg.start_yields.size()) ? cfg.start_yields[k] : 0.0;
        pi_sum += Nk * pk;
      }

      if (!(pi_sum > 0.0) || !std::isfinite(pi_sum)) {
        ++n_zero_pi;
        std::cout << "[run_nll_fit][zero-pi] Ee=" << ev.Ee
                  << " Eg=" << ev.Eg
                  << " t=" << ev.t
                  << " phi_e=" << ev.phi_detector_e
                  << " phi_g=" << ev.phi_detector_g
                  << " |";
        for (std::size_t k = 0; k < npar; ++k) {
          std::cout << " " << components[k].name << "=" << pks[k];
        }
#ifdef P2MEG_DEBUG_PRINT_ZERO_PI_Q2
        const double q2_over_m2 = RmdQ2OverM2(ev.Ee, ev.Eg, ev.phi_detector_e, ev.phi_detector_g);
        const double mmu = kMassesPDG.m_mu; // [MeV]
        const double q2 = (std::isfinite(q2_over_m2) && std::isfinite(mmu)) ? (q2_over_m2 * mmu * mmu)
                                                                            : std::numeric_limits<double>::quiet_NaN();
        std::cout << " q2_over_m2=" << q2_over_m2
                  << " q2=" << q2;
#endif
        std::cout << " pi=" << pi_sum << "\n";
      }
    }

    std::cout << "[run_nll_fit][zero-pi] events with pi<=0: "
              << n_zero_pi << " / " << events_in_window.size() << "\n";
#endif

  std::cout << "==================== Fit Result ====================\n";
  std::cout << "status   = " << res.status << "\n";
  std::cout << "nll_min  = " << res.nll_min << "\n";

  if (res.yields_hat.size() >= 3) {
    std::cout << "N_sig_hat = " << res.yields_hat[0] << "\n";
    std::cout << "N_rmd_hat = " << res.yields_hat[1] << "\n";
    std::cout << "N_acc_hat = " << res.yields_hat[2] << "\n";
  }

  if (res.yields_err.size() >= 3) {
    std::cout << "err_sig   = " << res.yields_err[0] << "\n";
    std::cout << "err_rmd   = " << res.yields_err[1] << "\n";
    std::cout << "err_acc   = " << res.yields_err[2] << "\n";
  }

  std::cout << "====================================================\n";
  return 0;
}
