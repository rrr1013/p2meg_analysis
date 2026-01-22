// src/HistUtils.cc
#include "p2meg/HistUtils.h"

#include <cmath>
#include <vector>

#include "p2meg/MathUtils.h"

int Hist_AxisBracketUniform(const TAxis& ax, double x, int& i0, int& i1, double& f) {
  const int n = ax.GetNbins();
  if (n < 2) return 1;

  const double xmin = ax.GetXmin();
  const double xmax = ax.GetXmax();
  const double dx = (xmax - xmin) / n;
  if (!(dx > 0.0) || !Math_IsFinite(dx)) return 2;

  // ビン中心基準の連続座標（1番ビン中心が 1.0）
  double i_float = (x - xmin) / dx + 0.5;

  // i0 と i1=i0+1 が必要なので i0 は 1..n-1
  if (i_float < 1.0) i_float = 1.0;
  const double eps = 1e-12;
  const double max_ifloat = static_cast<double>(n) - eps; // < n
  if (i_float > max_ifloat) i_float = max_ifloat;

  i0 = static_cast<int>(std::floor(i_float));
  if (i0 < 1) i0 = 1;
  if (i0 > n - 1) i0 = n - 1;

  i1 = i0 + 1;
  f = i_float - static_cast<double>(i0);
  if (f < 0.0) f = 0.0;
  if (f > 1.0) f = 1.0;

  return 0;
}

double Hist_InterpEeEg4(const THnD& h, double Ee, double Eg, int bin_phi_e, int bin_phi_g) {
  const TAxis* axE = h.GetAxis(0);
  const TAxis* axG = h.GetAxis(1);

  int i0, i1, j0, j1;
  double fE, fG;

  if (Hist_AxisBracketUniform(*axE, Ee, i0, i1, fE) != 0) return 0.0;
  if (Hist_AxisBracketUniform(*axG, Eg, j0, j1, fG) != 0) return 0.0;

  std::vector<int> idx(4, 1);

  auto get = [&](int ie, int ig) -> double {
    idx[0] = ie;
    idx[1] = ig;
    idx[2] = bin_phi_e;
    idx[3] = bin_phi_g;
    const Long64_t bin = h.GetBin(idx.data());
    const double v = h.GetBinContent(bin);
    if (v > 0.0 && Math_IsFinite(v)) return v;
    return 0.0;
  };

  const double v00 = get(i0, j0);
  const double v10 = get(i1, j0);
  const double v01 = get(i0, j1);
  const double v11 = get(i1, j1);

  const double a0 = (1.0 - fE) * v00 + fE * v10;
  const double a1 = (1.0 - fE) * v01 + fE * v11;
  const double v  = (1.0 - fG) * a0  + fG * a1;

  return (v > 0.0 && Math_IsFinite(v)) ? v : 0.0;
}

double Hist_SumAllBins4(const THnD& h) {
  const int n0 = h.GetAxis(0)->GetNbins();
  const int n1 = h.GetAxis(1)->GetNbins();
  const int n2 = h.GetAxis(2)->GetNbins();
  const int n3 = h.GetAxis(3)->GetNbins();

  std::vector<int> idx(4, 1);
  double sum = 0.0;

  for (int i0 = 1; i0 <= n0; ++i0) {
    idx[0] = i0;
    for (int i1 = 1; i1 <= n1; ++i1) {
      idx[1] = i1;
      for (int i2 = 1; i2 <= n2; ++i2) {
        idx[2] = i2;
        for (int i3 = 1; i3 <= n3; ++i3) {
          idx[3] = i3;
          const Long64_t bin = h.GetBin(idx.data());
          sum += h.GetBinContent(bin);
        }
      }
    }
  }
  return sum;
}
