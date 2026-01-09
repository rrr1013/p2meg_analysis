// scripts/inspect_rmd_grid.cc
//
// 使い方:
//   g++ -std=c++17 scripts/inspect_rmd_grid.cc $(root-config --cflags --libs) -o build/inspect_rmd_grid
//   ./build/inspect_rmd_grid data/pdf_cache/rmd_grid.root rmd_grid
//
// rmd_grid.root に保存されたメタ情報を表示する。

#include <iostream>
#include <string>

#include "TFile.h"
#include "TNamed.h"
#include "TParameter.h"

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

  f.Close();
  return 0;
}
