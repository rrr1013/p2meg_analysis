#include <iostream>
#include <iomanip>

#include "p2meg/SignalPdf.h"
#include "p2meg/AnalysisWindow.h"
#include "p2meg/DetectorResolution.h"
#include "p2meg/Constants.h"

static void PrintPoint(const char* label,
                       double Ee, double Eg, double t, double th,
                       const AnalysisWindow4D& win,
                       const DetectorResolutionConst& res) {
  const double p = SignalPdf(Ee, Eg, t, th, win, res);
  std::cout << std::setw(12) << label
            << "  Ee=" << std::setw(7) << Ee
            << "  Eg=" << std::setw(7) << Eg
            << "  t="  << std::setw(7) << t
            << "  th=" << std::setw(9) << th
            << "  p="  << std::setprecision(12) << p
            << "\n";
}

int main() {
  // 既存の解析窓・分解能（名前を analysis_window / detres に変更済みの想定）
  const AnalysisWindow4D& win = analysis_window;
  const DetectorResolutionConst& res = detres;

  const double Ee0 = 0.5 * kMassesPDG.m_mu;
  const double Eg0 = 0.5 * kMassesPDG.m_mu;
  const double t0  = res.t_mean;
  const double th0 = pi;

  std::cout << "=== SignalPdf quick test ===\n";
  std::cout << "Ee0=" << Ee0 << " Eg0=" << Eg0 << " t0=" << t0 << " th0=" << th0 << "\n";
  std::cout << "win: Ee[" << win.Ee_min << "," << win.Ee_max << "]"
            << " Eg[" << win.Eg_min << "," << win.Eg_max << "]"
            << " t["  << win.t_min  << "," << win.t_max  << "]"
            << " th[" << win.theta_min << "," << win.theta_max << "]\n";
  std::cout << "res: sEe=" << res.sigma_Ee
            << " sEg=" << res.sigma_Eg
            << " st="  << res.sigma_t
            << " sth=" << res.sigma_theta
            << " t_mean=" << res.t_mean << "\n\n";

  // 1) 中心付近（窓に入っていれば最大級）
  PrintPoint("center", Ee0, Eg0, t0, th0, win, res);

  // 2) theta: δ=pi-theta の half-normal（折り返し）なので pi から離れるほど下がるはず
  PrintPoint("theta-1s", Ee0, Eg0, t0, pi - 1.0 * res.sigma_theta, win, res);
  PrintPoint("theta-2s", Ee0, Eg0, t0, pi - 2.0 * res.sigma_theta, win, res);

  // 3) Ee/Eg をずらす
  PrintPoint("Ee+1s", Ee0 + 1.0 * res.sigma_Ee, Eg0, t0, th0, win, res);
  PrintPoint("Eg+1s", Ee0, Eg0 + 1.0 * res.sigma_Eg, t0, th0, win, res);

  // 4) t をずらす
  PrintPoint("t+1s", Ee0, Eg0, t0 + 1.0 * res.sigma_t, th0, win, res);

  // 5) 窓外（0）
  PrintPoint("out_Ee", win.Ee_max + 1.0, Eg0, t0, th0, win, res);
  PrintPoint("out_th", Ee0, Eg0, t0, 0.0, win, res);

  return 0;
}
