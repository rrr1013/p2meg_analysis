// scripts/inspect_rmd_grid.cc
//
// 使い方:
//   ./build/inspect_rmd_grid data/pdf_cache/rmd_grid.root rmd_grid
//
// rmd_grid.root に保存された 4D 格子PDF（key）と
// メタ情報（key+"_meta", key+"_d_min", key+"_seed", key+"_N_theta",
//           key+"_N_phi_e", key+"_N_phi_g", key+"_phi_e_min/max", key+"_phi_g_min/max",
//           key+"_P_mu"）を表示する。
// 格子の軸は (Ee, Eg, phi_detector_e, phi_detector_g) を想定。

#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include "TFile.h"
#include "TNamed.h"
#include "TParameter.h"
#include "THn.h"
#include "TAxis.h"
#include "TArrayD.h"

// 密度格子の体積積分（4D）と簡易統計を計算する
static void SumDensityStats4D(const THnD& h, double& sum, long& n_negative, long& n_nonfinite) {
  const TAxis* ax0 = h.GetAxis(0);
  const TAxis* ax1 = h.GetAxis(1);
  const TAxis* ax2 = h.GetAxis(2);
  const TAxis* ax3 = h.GetAxis(3);

  sum = 0.0;
  n_negative = 0;
  n_nonfinite = 0;

  if (!ax0 || !ax1 || !ax2 || !ax3) return;

  const int n0 = ax0->GetNbins();
  const int n1 = ax1->GetNbins();
  const int n2 = ax2->GetNbins();
  const int n3 = ax3->GetNbins();

  std::vector<int> idx(4, 1);
  for (int i0 = 1; i0 <= n0; ++i0) {
    idx[0] = i0;
    const double w0 = ax0->GetBinWidth(i0);
    for (int i1 = 1; i1 <= n1; ++i1) {
      idx[1] = i1;
      const double w1 = ax1->GetBinWidth(i1);
      for (int i2 = 1; i2 <= n2; ++i2) {
        idx[2] = i2;
        const double w2 = ax2->GetBinWidth(i2);
        for (int i3 = 1; i3 <= n3; ++i3) {
          idx[3] = i3;
          const double w3 = ax3->GetBinWidth(i3);
          const double vol = w0 * w1 * w2 * w3;
          if (!(vol > 0.0)) continue;

          const Long64_t bin = h.GetBin(idx.data());
          const double v = h.GetBinContent(bin);
          if (!std::isfinite(v)) {
            ++n_nonfinite;
            continue;
          }
          if (v < 0.0) ++n_negative;
          sum += v * vol;
        }
      }
    }
  }
}

// THnD の軸情報を出力する
static void PrintAxisInfo(const TAxis* ax, const char* label) {
  if (!ax || !label) return;

  std::cout << "[inspect_rmd_grid] axis " << label
            << ": nbins=" << ax->GetNbins()
            << " range=[" << ax->GetXmin() << "," << ax->GetXmax() << "]"
            << " title=\"" << ax->GetTitle() << "\"\n";

  const TArrayD* edges = ax->GetXbins();
  if (edges && edges->GetSize() > 0) {
    std::cout << "[inspect_rmd_grid] axis " << label << " edges (" << edges->GetSize() << "): ";
    for (int i = 0; i < edges->GetSize(); ++i) {
      std::cout << std::setprecision(6) << edges->At(i);
      if (i + 1 < edges->GetSize()) std::cout << ", ";
    }
    std::cout << "\n";
  } else {
    std::cout << "[inspect_rmd_grid] axis " << label
              << " bin_width=" << ax->GetBinWidth(1) << "\n";
  }
}

// 4D 格子の概要を出力する
static void PrintGridSummary(const THnD& h, const char* label) {
  std::cout << "[inspect_rmd_grid] " << label << " hist: name=" << h.GetName()
            << " title=\"" << h.GetTitle() << "\"\n";
  std::cout << "[inspect_rmd_grid] " << label << " ndim=" << h.GetNdimensions()
            << " entries=" << h.GetEntries() << "\n";

  if (h.GetNdimensions() == 4) {
    PrintAxisInfo(h.GetAxis(0), "0(Ee)");
    PrintAxisInfo(h.GetAxis(1), "1(Eg)");
    PrintAxisInfo(h.GetAxis(2), "2(phi_detector_e)");
    PrintAxisInfo(h.GetAxis(3), "3(phi_detector_g)");
    double mass = 0.0;
    long n_negative = 0;
    long n_nonfinite = 0;
    SumDensityStats4D(h, mass, n_negative, n_nonfinite);
    std::cout << "[inspect_rmd_grid] " << label
              << " density integral (sum v*vol)=" << std::setprecision(10)
              << mass << "\n";
    std::cout << "[inspect_rmd_grid] " << label
              << " negative bins=" << n_negative
              << " nonfinite bins=" << n_nonfinite << "\n";
  } else {
    std::cout << "[inspect_rmd_grid] " << label
              << " unexpected ndim (expected 4)\n";
  }
}

int main(int argc, char** argv)
{
  const char* filepath = (argc >= 2) ? argv[1] : "data/pdf_cache/rmd_grid.root";
  const char* key = (argc >= 3) ? argv[2] : "rmd_grid";

  TFile f(filepath, "READ");
  if (f.IsZombie()) {
    std::cerr << "[inspect_rmd_grid] cannot open: " << filepath << "\n";
    return 1;
  }

  std::cout << "[inspect_rmd_grid] file: " << filepath << "\n";
  std::cout << "[inspect_rmd_grid] key: " << key << "\n";

  const std::string meta_name = std::string(key) + "_meta";
  TNamed* meta = nullptr;
  f.GetObject(meta_name.c_str(), meta);
  if (meta) {
    std::cout << "[inspect_rmd_grid] meta:\n";
    std::cout << meta->GetTitle() << "\n";
  } else {
    std::cout << "[inspect_rmd_grid] meta not found: " << meta_name << "\n";
  }

  const std::string dmin_name = std::string(key) + "_d_min";
  const std::string seed_name = std::string(key) + "_seed";
  const std::string ntheta_name = std::string(key) + "_N_theta";
  const std::string nphi_e_name = std::string(key) + "_N_phi_e";
  const std::string nphi_g_name = std::string(key) + "_N_phi_g";
  const std::string phi_e_min_name = std::string(key) + "_phi_e_min";
  const std::string phi_e_max_name = std::string(key) + "_phi_e_max";
  const std::string phi_g_min_name = std::string(key) + "_phi_g_min";
  const std::string phi_g_max_name = std::string(key) + "_phi_g_max";
  const std::string pmu_name = std::string(key) + "_P_mu";

  TParameter<double>* par_dmin = nullptr;
  f.GetObject(dmin_name.c_str(), par_dmin);
  if (par_dmin) {
    std::cout << "[inspect_rmd_grid] d_min=" << par_dmin->GetVal() << "\n";
  } else {
    std::cout << "[inspect_rmd_grid] d_min not found: " << dmin_name << "\n";
  }

  TParameter<Long64_t>* par_seed = nullptr;
  f.GetObject(seed_name.c_str(), par_seed);
  if (par_seed) {
    std::cout << "[inspect_rmd_grid] seed=" << par_seed->GetVal() << "\n";
  } else {
    std::cout << "[inspect_rmd_grid] seed not found: " << seed_name << "\n";
  }

  TParameter<int>* par_ntheta = nullptr;
  f.GetObject(ntheta_name.c_str(), par_ntheta);
  if (par_ntheta) {
    std::cout << "[inspect_rmd_grid] N_theta=" << par_ntheta->GetVal() << "\n";
  } else {
    std::cout << "[inspect_rmd_grid] N_theta not found: " << ntheta_name << "\n";
  }

  TParameter<int>* par_nphi_e = nullptr;
  f.GetObject(nphi_e_name.c_str(), par_nphi_e);
  if (par_nphi_e) {
    std::cout << "[inspect_rmd_grid] N_phi_e=" << par_nphi_e->GetVal() << "\n";
  } else {
    std::cout << "[inspect_rmd_grid] N_phi_e not found: " << nphi_e_name << "\n";
  }

  TParameter<int>* par_nphi_g = nullptr;
  f.GetObject(nphi_g_name.c_str(), par_nphi_g);
  if (par_nphi_g) {
    std::cout << "[inspect_rmd_grid] N_phi_g=" << par_nphi_g->GetVal() << "\n";
  } else {
    std::cout << "[inspect_rmd_grid] N_phi_g not found: " << nphi_g_name << "\n";
  }

  TParameter<double>* par_phi_e_min = nullptr;
  TParameter<double>* par_phi_e_max = nullptr;
  f.GetObject(phi_e_min_name.c_str(), par_phi_e_min);
  f.GetObject(phi_e_max_name.c_str(), par_phi_e_max);
  if (par_phi_e_min && par_phi_e_max) {
    std::cout << "[inspect_rmd_grid] phi_e range=[" << par_phi_e_min->GetVal()
              << "," << par_phi_e_max->GetVal() << "]\n";
  } else {
    std::cout << "[inspect_rmd_grid] phi_e range not found: " << phi_e_min_name
              << ", " << phi_e_max_name << "\n";
  }

  TParameter<double>* par_phi_g_min = nullptr;
  TParameter<double>* par_phi_g_max = nullptr;
  f.GetObject(phi_g_min_name.c_str(), par_phi_g_min);
  f.GetObject(phi_g_max_name.c_str(), par_phi_g_max);
  if (par_phi_g_min && par_phi_g_max) {
    std::cout << "[inspect_rmd_grid] phi_g range=[" << par_phi_g_min->GetVal()
              << "," << par_phi_g_max->GetVal() << "]\n";
  } else {
    std::cout << "[inspect_rmd_grid] phi_g range not found: " << phi_g_min_name
              << ", " << phi_g_max_name << "\n";
  }

  TParameter<double>* par_pmu = nullptr;
  f.GetObject(pmu_name.c_str(), par_pmu);
  if (par_pmu) {
    std::cout << "[inspect_rmd_grid] P_mu=" << par_pmu->GetVal() << "\n";
  } else {
    std::cout << "[inspect_rmd_grid] P_mu not found: " << pmu_name << "\n";
  }

  TObject* obj = f.Get(key);
  if (!obj) {
    std::cout << "[inspect_rmd_grid] grid not found: " << key << "\n";
  } else {
    THnD* h = dynamic_cast<THnD*>(obj);
    if (!h) {
      std::cout << "[inspect_rmd_grid] object is not THnD: " << key << "\n";
    } else {
      PrintGridSummary(*h, "grid");
    }
  }

  f.Close();
  return 0;
}
