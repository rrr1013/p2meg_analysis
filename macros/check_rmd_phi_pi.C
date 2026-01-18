// macros/check_rmd_phi_pi.C
// 使い方例:
//   root -l -q 'macros/check_rmd_phi_pi.C("rmd_grid.root","rmd")'
//   root -l -q 'macros/check_rmd_phi_pi.C("build/rmd_grid.root","rmd_grid")'

#include "TFile.h"
#include "THn.h"
#include "TAxis.h"
#include "TMath.h"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>

static double SumRange4(THnD* h,
                        int pe0, int pe1,  // [pe0, pe1]  (0..nPe+1)
                        int pg0, int pg1)  // [pg0, pg1]
{
  const int nEe = h->GetAxis(0)->GetNbins();
  const int nEg = h->GetAxis(1)->GetNbins();

  double sum = 0.0;
  int idx[4];

  for (int iEe = 1; iEe <= nEe; ++iEe) {
    idx[0] = iEe;
    for (int iEg = 1; iEg <= nEg; ++iEg) {
      idx[1] = iEg;
      for (int iPe = pe0; iPe <= pe1; ++iPe) {
        idx[2] = iPe;
        for (int iPg = pg0; iPg <= pg1; ++iPg) {
          idx[3] = iPg;
          const Long64_t bin = h->GetBin(idx);
          const double c = h->GetBinContent(bin);
          if (c > 0.0 && std::isfinite(c)) sum += c;
        }
      }
    }
  }
  return sum;
}

void check_rmd_phi_pi(const char* filepath = "rmd_grid.root",
                      const char* key      = "rmd")
{
  TFile f(filepath, "READ");
  if (f.IsZombie()) {
    std::cerr << "[check_rmd_phi_pi] cannot open: " << filepath << "\n";
    return;
  }

  THnD* h = dynamic_cast<THnD*>(f.Get(key));
  if (!h) {
    std::cerr << "[check_rmd_phi_pi] THnD not found: key=" << key << "\n";
    return;
  }
  if (h->GetNdimensions() != 4) {
    std::cerr << "[check_rmd_phi_pi] THnD dim != 4\n";
    return;
  }

  const TAxis* axPe = h->GetAxis(2);
  const TAxis* axPg = h->GetAxis(3);
  const int nPe = axPe->GetNbins();
  const int nPg = axPg->GetNbins();

  const double pi = TMath::Pi();
  const double pi_minus = std::nextafter(pi, 0.0);

  std::cout << std::setprecision(17);

  std::cout << "=== Axis info ===\n";
  std::cout << "phi_e: nbins=" << nPe
            << " xmin=" << axPe->GetXmin()
            << " xmax=" << axPe->GetXmax()
            << " overflow_index=" << (nPe + 1) << "\n";
  std::cout << "phi_g: nbins=" << nPg
            << " xmin=" << axPg->GetXmin()
            << " xmax=" << axPg->GetXmax()
            << " overflow_index=" << (nPg + 1) << "\n";

  std::cout << "\n=== FindBin check (phi_e axis) ===\n";
  std::cout << "pi       =" << pi << "  FindBin(pi)      =" << axPe->FindBin(pi) << "\n";
  std::cout << "pi_minus =" << pi_minus << "  FindBin(pi_minus)=" << axPe->FindBin(pi_minus) << "\n";
  std::cout << "last in-range bin edges: ["
            << axPe->GetBinLowEdge(nPe) << ", " << axPe->GetBinUpEdge(nPe) << ")\n";

  std::cout << "\n=== Mass check (raw bin contents, before density normalization is irrelevant here) ===\n";

  // in-range mass (あなたの正規化ループと同じ：phi も 1..n のみ)
  const double mass_inrange =
      SumRange4(h, 1, nPe, 1, nPg);

  // overflow に落ちている質量（phi_e overflow を含む）
  const double mass_pe_over =
      SumRange4(h, nPe + 1, nPe + 1, 1, nPg);

  // overflow に落ちている質量（phi_g overflow を含む）
  const double mass_pg_over =
      SumRange4(h, 1, nPe, nPg + 1, nPg + 1);

  // 両方 overflow（phi_e overflow かつ phi_g overflow）
  const double mass_both_over =
      SumRange4(h, nPe + 1, nPe + 1, nPg + 1, nPg + 1);

  // 「最後の in-range ビン」に入っている質量（phi_e 最終ビン、phi_g は全 in-range）
  const double mass_pe_last =
      SumRange4(h, nPe, nPe, 1, nPg);

  // 「最後の in-range ビン」に入っている質量（phi_g 最終ビン、phi_e は全 in-range）
  const double mass_pg_last =
      SumRange4(h, 1, nPe, nPg, nPg);

  std::cout << "mass_inrange      = " << mass_inrange << "\n";
  std::cout << "mass_pe_last      = " << mass_pe_last << "   (phi_e bin = last in-range)\n";
  std::cout << "mass_pg_last      = " << mass_pg_last << "   (phi_g bin = last in-range)\n";
  std::cout << "mass_pe_overflow  = " << mass_pe_over << "   (phi_e bin = overflow)\n";
  std::cout << "mass_pg_overflow  = " << mass_pg_over << "   (phi_g bin = overflow)\n";
  std::cout << "mass_both_overflow= " << mass_both_over << "   (both overflow)\n";

  std::cout << "\nInterpretation:\n";
  std::cout << " - FindBin(pi) が (nbins+1) なら、phi=pi は overflow 扱いです。\n";
  std::cout << " - mass_pe_overflow / mass_pg_overflow が 0 でなく、mass_*_last が小さいなら、phi=pi 由来が overflow に落ちています。\n";
}
