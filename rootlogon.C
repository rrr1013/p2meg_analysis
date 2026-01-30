// rootlogon.C
// p2MEG: ACLiC 生成物を build/root_aclic にまとめる

#include "TSystem.h"

{
    const char* build_dir = "build/root_aclic";
    gSystem->mkdir(build_dir, true);
    gSystem->SetBuildDir(build_dir, true);
}
