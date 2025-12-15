#include <iostream>
#include "p2meg/MakeRMDGridPdf.h"

int main() {
  const char* out = "data/pdf_cache/rmd_grid.root";
  const char* key = "rmd_grid";

  const int rc = MakeRMDGridPdf(out, key);
  std::cout << "MakeRMDGridPdf returned " << rc << std::endl;
  return rc;
}
