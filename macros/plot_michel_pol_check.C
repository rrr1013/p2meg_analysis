// macros/plot_michel_pol_check.C
//
// Michel 偏極測定（別実験）チェック用プロット
//  - 4データセット（入替え運用）から cross-ratio 非対称 A_data(E) を作る
//  - 理論テンプレート K(E) を作る（Michel+エネルギー応答の畳み込み）
//  - A_data(E) と A_th(E)=P_hat*K(E) を同一図で比較
//  - 併せて K(E) と pull(E) も表示
//
// 実行例（ROOT）:
//   root -l
//   .L src/MichelEData.cc+
//   .L src/MichelPolTemplate.cc+
//   .L src/MichelPolFit.cc+        // (必須ではないが、依存関係があるなら先にロード)
//   .L macros/plot_michel_pol_check.C+
//   PlotMichelPolCheck("A_plus.dat","A_minus.dat","B_plus.dat","B_minus.dat");
//
// 注意:
//  - ACLiC は別共有ライブラリになるので、必要な .cc を先に .L ...+ でロードする。
//  - 本マクロは Ee-only フォーマット（1列Ee[MeV]）を想定。

#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TAxis.h"

#include "p2meg/MichelEData.h"
#include "p2meg/MichelPolConfig.h"
#include "p2meg/MichelPolFit.h"

// MichelPolTemplate.cc にある（ヘッダを作らない方針のため extern 宣言）
extern std::vector<double> BuildMichelPolKi(const MichelPolConfig& cfg);

static inline bool IsFinite(double x) { return std::isfinite(x); }

static inline int BinIndex(double Ee, double Emin, double Emax, int nbins) {
  if (!(Emax > Emin) || !(nbins > 0)) return -1;
  if (!(Ee >= Emin) || !(Ee < Emax)) return -1;
  const double binw = (Emax - Emin) / static_cast<double>(nbins);
  const int idx = static_cast<int>((Ee - Emin) / binw);
  if (idx < 0 || idx >= nbins) return -1;
  return idx;
}

static inline double BinCenter(int i, double Emin, double Emax, int nbins) {
  const double binw = (Emax - Emin) / static_cast<double>(nbins);
  return Emin + (i + 0.5) * binw;
}

static void FillCounts(const std::vector<MichelEEvent>& evs,
                       const MichelPolConfig& cfg,
                       std::vector<double>& counts) {
  for (const auto& ev : evs) {
    const int ib = BinIndex(ev.Ee, cfg.Ee_min, cfg.Ee_max, cfg.nbins_Ee);
    if (ib >= 0) counts[ib] += 1.0;
  }
}

// cross-ratio 非対称と誤差（src/MichelPolFit.cc と同じ近似）
static bool CrossRatioAsymAndErr(double a, double b, double c, double d,
                                 int min_counts_each,
                                 double& A, double& sigmaA) {
  if (!(a >= 0.0 && b >= 0.0 && c >= 0.0 && d >= 0.0)) return false;
  if (a < min_counts_each || b < min_counts_each || c < min_counts_each || d < min_counts_each) {
    return false;
  }

  const double ab = a * b;
  const double cd = c * d;
  if (!(ab > 0.0 && cd > 0.0) || !IsFinite(ab) || !IsFinite(cd)) return false;

  const double r = std::sqrt(ab / cd);
  if (!(r > 0.0) || !IsFinite(r)) return false;

  A = (r - 1.0) / (r + 1.0);
  if (!IsFinite(A)) return false;

  // Var(ln r) ≈ (1/4)(1/a+1/b+1/c+1/d),  u=(1/2)ln r,  A=tanh(u)
  const double var_ln_r = 0.25 * (1.0/a + 1.0/b + 1.0/c + 1.0/d);
  const double var_u    = 0.25 * var_ln_r;
  const double one_minus_A2 = 1.0 - A*A;
  const double var_A = (one_minus_A2 * one_minus_A2) * var_u;

  if (!(var_A > 0.0) || !IsFinite(var_A)) return false;
  sigmaA = std::sqrt(var_A);
  return (sigmaA > 0.0 && IsFinite(sigmaA));
}

static inline bool InFitRange(double Ecen, const MichelPolConfig& cfg) {
  return (Ecen >= cfg.fit_Ee_min && Ecen <= cfg.fit_Ee_max);
}

void PlotMichelPolCheck(const char* path_A_plus,
                        const char* path_A_minus,
                        const char* path_B_plus,
                        const char* path_B_minus) {
  gStyle->SetOptStat(0);

  const auto& cfg = michel_pol_config;
  const int nb = cfg.nbins_Ee;
  if (!(nb > 0) || !(cfg.Ee_max > cfg.Ee_min)) {
    std::cout << "[PlotMichelPolCheck] invalid config.\n";
    return;
  }

  // ---- 読み込み ----
  long long r1=0,s1=0, r2=0,s2=0, r3=0,s3=0, r4=0,s4=0;
  const auto evAp = ReadMichelEData(path_A_plus,  &r1, &s1);
  const auto evAm = ReadMichelEData(path_A_minus, &r2, &s2);
  const auto evBp = ReadMichelEData(path_B_plus,  &r3, &s3);
  const auto evBm = ReadMichelEData(path_B_minus, &r4, &s4);

  std::cout << "=== Input summary ===\n";
  std::cout << "A_plus : read="  << r1 << " skipped=" << s1 << "\n";
  std::cout << "A_minus: read="  << r2 << " skipped=" << s2 << "\n";
  std::cout << "B_plus : read="  << r3 << " skipped=" << s3 << "\n";
  std::cout << "B_minus: read="  << r4 << " skipped=" << s4 << "\n";
  std::cout << "total  : read="  << (r1+r2+r3+r4) << " skipped=" << (s1+s2+s3+s4) << "\n\n";

  // ---- ヒスト（カウント配列） ----
  std::vector<double> nAp(nb,0.0), nAm(nb,0.0), nBp(nb,0.0), nBm(nb,0.0);
  FillCounts(evAp, cfg, nAp);
  FillCounts(evAm, cfg, nAm);
  FillCounts(evBp, cfg, nBp);
  FillCounts(evBm, cfg, nBm);

  // ---- 理論テンプレ K_i ----
  const std::vector<double> K = BuildMichelPolKi(cfg);
  if (static_cast<int>(K.size()) != nb) {
    std::cout << "[PlotMichelPolCheck] BuildMichelPolKi returned wrong size.\n";
    return;
  }

  // ---- A_data, sigma, K_i を作り、解析解で P_hat ----
  struct Pt { double E; double A; double s; double K; int ibin; };
  std::vector<Pt> pts_all;
  std::vector<Pt> pts_fit;

  double S1 = 0.0;
  double S2 = 0.0;

  for (int i = 0; i < nb; ++i) {
    const double Ecen = BinCenter(i, cfg.Ee_min, cfg.Ee_max, cfg.nbins_Ee);

    // cross-ratio: a=nA+(i), b=nB+(i), c=nA-(i), d=nB-(i)
    const double a = nAp[i];
    const double b = nBp[i];
    const double c = nAm[i];
    const double d = nBm[i];

    double A=0.0, sA=0.0;
    if (!CrossRatioAsymAndErr(a,b,c,d,cfg.min_counts_each,A,sA)) continue;

    const double Ki = K[i];
    if (!IsFinite(Ki)) continue;

    Pt p{Ecen, A, sA, Ki, i};
    pts_all.push_back(p);

    if (InFitRange(Ecen, cfg)) {
      const double w = 1.0 / (sA*sA);
      S1 += Ki * A * w;
      S2 += Ki * Ki * w;
      pts_fit.push_back(p);
    }
  }

  if (!(S2 > 0.0) || pts_fit.size() < 2) {
    std::cout << "[PlotMichelPolCheck] not enough bins for fit.\n";
    return;
  }

  const double P_hat = S1 / S2;
  const double errP  = 1.0 / std::sqrt(S2);

  double chi2 = 0.0;
  for (const auto& p : pts_fit) {
    const double pull = (p.A - P_hat * p.K) / p.s;
    chi2 += pull * pull;
  }
  const int ndf = static_cast<int>(pts_fit.size()) - 1;

  std::cout << "=== Quick fit (macro) ===\n";
  std::cout << std::fixed << std::setprecision(6);
  std::cout << "P_mu_hat = " << P_hat << "  +/- " << errP << "\n";
  if (ndf > 0) std::cout << "chi2/ndf = " << chi2 << " / " << ndf << " = " << (chi2 / ndf) << "\n";
  std::cout << "bins used (fit range) = " << pts_fit.size() << "\n\n";

  // ---- A_data と A_th のプロット ----
  TCanvas* c1 = new TCanvas("c_michel_pol_asym", "Michel pol: asymmetry", 900, 650);

  TGraphErrors* gA = new TGraphErrors();
  gA->SetTitle("Michel polarization check;E_{e} [MeV];A_{data}(E)");
  gA->SetMarkerStyle(20);
  gA->SetMarkerSize(0.9);

  for (int ip = 0; ip < (int)pts_all.size(); ++ip) {
    const auto& p = pts_all[ip];
    gA->SetPoint(ip, p.E, p.A);
    gA->SetPointError(ip, 0.0, p.s);
  }

  // 理論: A_th(E)=P_hat*K_i（同じ点に打つ）
  TGraph* gTh = new TGraph();
  gTh->SetLineWidth(2);

  for (int ip = 0; ip < (int)pts_all.size(); ++ip) {
    const auto& p = pts_all[ip];
    gTh->SetPoint(ip, p.E, P_hat * p.K);
  }

  gA->Draw("AP");
  gA->GetYaxis()->SetTitleOffset(1.2);
  gTh->Draw("L SAME");

  // フィット範囲の目印
  TLine* l1 = new TLine(cfg.fit_Ee_min, gA->GetYaxis()->GetXmin(), cfg.fit_Ee_min, gA->GetYaxis()->GetXmax());
  TLine* l2 = new TLine(cfg.fit_Ee_max, gA->GetYaxis()->GetXmin(), cfg.fit_Ee_max, gA->GetYaxis()->GetXmax());
  l1->SetLineStyle(2);
  l2->SetLineStyle(2);
  // 注意: TGraph の軸範囲取得が面倒なので、描画後に pad の範囲で引く
  const double y1 = gPad->GetUymin();
  const double y2 = gPad->GetUymax();
  l1->SetY1(y1); l1->SetY2(y2);
  l2->SetY1(y1); l2->SetY2(y2);
  l1->Draw("SAME");
  l2->Draw("SAME");

  TLegend* leg = new TLegend(0.58, 0.74, 0.88, 0.88);
  leg->AddEntry(gA, "A_{data} (cross-ratio)", "pe");
  leg->AddEntry(gTh, "A_{th}=#hat{P}_{#mu} K(E)", "l");
  leg->AddEntry((TObject*)0, Form("#hat{P}_{#mu}=%.4f #pm %.4f", P_hat, errP), "");
  if (ndf > 0) leg->AddEntry((TObject*)0, Form("#chi^{2}/ndf=%.2f/%d", chi2, ndf), "");
  leg->Draw();

  c1->Update();

  // ---- K(E) のプロット ----
  TCanvas* c2 = new TCanvas("c_michel_pol_K", "Michel pol: K(E)", 900, 650);

  TGraph* gK = new TGraph();
  gK->SetTitle("Sensitivity curve K(E);E_{e} [MeV];K(E)=|cos#theta| V/U");
  gK->SetLineWidth(2);

  int ipk = 0;
  for (int i = 0; i < nb; ++i) {
    const double Ecen = BinCenter(i, cfg.Ee_min, cfg.Ee_max, cfg.nbins_Ee);
    const double Ki = K[i];
    if (!IsFinite(Ki)) continue;
    gK->SetPoint(ipk++, Ecen, Ki);
  }

  gK->Draw("AL");
  c2->Update();

  // ---- pull(E) のプロット（フィット範囲内だけ）----
  TCanvas* c3 = new TCanvas("c_michel_pol_pull", "Michel pol: pulls", 900, 650);

  TGraph* gP = new TGraph();
  gP->SetTitle("Pulls in fit range;E_{e} [MeV];(A_{data}-A_{th})/#sigma_{A}");

  int ipp = 0;
  for (const auto& p : pts_fit) {
    const double pull = (p.A - P_hat * p.K) / p.s;
    gP->SetPoint(ipp++, p.E, pull);
  }
  gP->SetMarkerStyle(20);
  gP->SetMarkerSize(0.9);
  gP->Draw("AP");

  TLine* l0 = new TLine(cfg.fit_Ee_min, 0.0, cfg.fit_Ee_max, 0.0);
  l0->SetLineStyle(2);
  l0->Draw("SAME");

  c3->Update();

  std::cout << "[PlotMichelPolCheck] canvases: c_michel_pol_asym, c_michel_pol_K, c_michel_pol_pull\n";
}
