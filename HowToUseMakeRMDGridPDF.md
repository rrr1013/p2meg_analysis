## RMD PDF の仕組み（MakeRMDGridPdf / RMDGridPdf）

### 目的
RMD の理論式 `RMD_d3B_dEe_dEg_dcos(Ee, Eg, cosθ)` と検出器分解能（Ee, Eg, θ の独立ガウシアン）を用いて、解析窓内で正規化された RMD の PDF を評価できるようにする。

尤度解析の際、毎回畳み込み計算を行うのは計算量が多すぎるため、事前にモンテカルロ法でPDFを作っておき、それを参照するという形式で計算する。

この実装では計算量と運用性のため、4変数 PDF を次のように分離して扱う。

\[
P(E_e,E_\gamma,t,\theta)=p_3(E_e,E_\gamma,\theta)\times p_t(t)
\]

- \(p_3(E_e,E_\gamma,\theta)\): Ee, Eg, θ を変数とする 3D の格子 PDF（オフラインで ROOT ファイルに保存）
- \(p_t(t)\): 時間差 t の PDF（窓内正規化ガウシアン）。4DのPDFを作るときに解析的に掛ける（格子には含めない）
- 解析窓（Ee, Eg, t, θ）と分解能（σEe, σEg, σt, σθ, t_mean）は `include/p2meg/...` にある設定から参照される想定

### 何が生成されるか
`MakeRMDGridPdf` は以下を作る。

- 解析窓内で RMD 理論式をサンプリングし、Ee, Eg, θ を分解能でスメアして 3D ヒスト（THn）に充填
- 3D ヒストを窓内で正規化して「3D PDF（密度）」として保存
- 生成条件（ビニング、乱数 seed、d_min、解析窓、分解能など）をメタ情報として同じ ROOT ファイルに保存

出力例:
- ファイル: `data/pdf_cache/rmd_grid.root`
- キー: `rmd_grid`（3D 格子 PDF）
- メタ: `rmd_grid_meta`（生成条件の記録）

### 注意点：　GitHub に ROOT ファイルが入らない点
`data/pdf_cache/rmd_grid.root` は生成物であり、`.gitignore` により GitHub に push されない。
そのため、リポジトリを clone した共同実験者は各自の環境で必ず一度 `MakeRMDGridPdf` を実行して格子 ROOT を作る必要がある。
格子 ROOT が無いと `RMDGridPdf_Load(...)` は失敗し、`RMDGridPdf(...)` は評価できない。

### 解析時の評価
解析時は以下の流れで PDF を評価する。

1. `RMDGridPdf_Load(filepath, key)` で 3D 格子 PDF（Ee, Eg, θ）をメモリにロード
2. イベントごとに `RMDGridPdf(Ee, Eg, t, theta)` を呼び、PDF 値を得る

内部では:
- 3D 格子から \(p_3(E_e,E_\gamma,\theta)\) を補間して取得
- \(p_t(t)\) を時間窓内正規化ガウシアンとして解析的に計算（t_mean を考慮）
- \(p_3 \times p_t\) を返す
- 解析窓外や未ロード等の場合は 0 を返す

### 正規化
- 3D 格子は解析窓内で
\[
\int dE_e\,dE_\gamma\,d\theta\;p_3(E_e,E_\gamma,\theta)=1
\]
となるように正規化されている。
- 時間 PDF は解析窓内で
\[
\int_{t_{\min}}^{t_{\max}} p_t(t)\,dt=1
\]
となるように正規化されている。
- よって 4D PDF も解析窓内で正規化される。

### 生成（MakeRMDGridPdf）の実行例
```cpp
#include <iostream>
#include <filesystem>
#include "p2meg/MakeRMDGridPdf.h"

int main() {
  std::filesystem::create_directories("data/pdf_cache");
  const char* out = "data/pdf_cache/rmd_grid.root";
  const char* key = "rmd_grid";

  const int rc = MakeRMDGridPdf(out, key);
  std::cout << "MakeRMDGridPdf returned " << rc << std::endl;
  return rc;
}
```

### 解析側（RMDGridPdf）の使用例
```cpp
#include <iostream>
#include "p2meg/RMDGridPdf.h"

int main() {
  const char* in  = "data/pdf_cache/rmd_grid.root";
  const char* key = "rmd_grid";

  if (!RMDGridPdf_Load(in, key)) {
    std::cerr << "RMDGridPdf_Load failed" << std::endl;
    return 1;
  }

  const double Ee = 50.0;   // [MeV]
  const double Eg = 48.0;   // [MeV]
  const double t  = 0.1;    // [ns]
  const double th = 3.05;   // [rad]

  const double p = RMDGridPdf(Ee, Eg, t, th);
  std::cout << "RMDGridPdf = " << p << std::endl;
  return 0;
}
```