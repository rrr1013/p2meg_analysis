# Functions

このファイルは p2meg_analysis の関数・クラスの「使い方」をまとめたものです。

# Contribution

新しい関数・クラスを追加した場合は、この`FUNCTIONS.md` に「使い方セクション」と Index の1行を追加してください。

## AIに「関数の使い方」を書かせるプロンプト

関数をAIで作った場合は以下をそのままAIに貼り付けてください。

---





### プロンプト本文

今実装した関数について、`p2meg_analysis/FUNCTIONS.md` に追記する「使い方」だけを書いてください。実装の変更提案や改善案は不要です。  

返答は コードブロック1つだけにしてください。コードブロックの外に文章を一切書かないでください（前置き・説明・箇条書き・空行の追加も禁止）。

そのコードブロックは 言語をmdにし、開始と終了は4つのバッククォートで囲ってください（````md 〜 ````）。

コードブロックの中身は、私が `FUNCTIONS.md` にそのまま貼り付けられる Markdown断片のみにしてください。

また、使用例はコードブロックの中で通常どおり
```cpp
...
```
の fenced code block を使ってください。

次の項目だけをこの順番で書いてください：
1. 関数名（またはクラス名）
2. ヘッダー名（`include/p2meg/...` のパス）
3. 目的
4. シグネチャ
5. 入力（引数ごと：意味を1行ずつ）
6. 出力（戻り値：意味を1〜2文）
7. 使用例（簡単なC++例。10〜20行程度で、コピペで動く形）

出力フォーマットは必ず次に従ってください：

```md
### <関数名>
- Header: `include/p2meg/<...>.h`
- 目的: ...

- シグネチャ

- 入力:
  - `<arg1>`: ...
  - `<arg2>`: ...

- 出力:
  - 戻り値: ...

- 使用例:
```cpp
#include <iostream>
#include "p2meg/<...>.h"

int main() {
  // ...
  return 0;
}
```






# Function list


### Michel_d2Shape_dE_dCosTheta
- Header: `include/p2meg/MichelSpectrum.h`
- 目的: ミシェル崩壊の二重微分分布に比例する非正規化 shape（全体定数は省略）を返します。

- シグネチャ
  - `double Michel_d2Shape_dE_dCosTheta(double Ee, double costh, double P_mu, const MichelParams& mp = kMichelSM, const ParticleMasses& ms = kMassesPDG);`

- 入力:
  - `Ee`: 陽電子エネルギー（MeV）
  - `costh`: $ \cos\theta $（$ \theta $ は偏極軸と陽電子運動量のなす角）
  - `P_mu`: ミューオン偏極度（角度項の係数。一般に [-1, 1] を想定）
  - `mp`: ミシェルパラメータ（`rho, eta, xi, delta`）。省略時は `kMichelSM`
  - `ms`: 粒子質量（`m_mu, m_e`）。省略時は `kMassesPDG`

- 出力:
  - 戻り値: 非正規化の shape 値（任意単位）。定義域外（`costh` が [-1,1] 外、または `x` が (0,1) 外）では 0 を返します。

- 使用例:
```cpp
#include <iostream>
#include "p2meg/MichelSpectrum.h"
#include "p2meg/Constants.h"

int main() {
  const double Emax = Michel_Emax(kMassesPDG);
  const double Ee = 0.7 * Emax;
  const double P_mu = 1.0;

  const double w_p = Michel_d2Shape_dE_dCosTheta(Ee, +1.0, P_mu);
  const double w_0 = Michel_d2Shape_dE_dCosTheta(Ee,  0.0, P_mu);
  const double w_m = Michel_d2Shape_dE_dCosTheta(Ee, -1.0, P_mu);

  std::cout << "w(+1) = " << w_p << "\n";
  std::cout << "w( 0) = " << w_0 << "\n";
  std::cout << "w(-1) = " << w_m << "\n";
  return 0;
}
```

### Michel_dShape_dE
- Header: `include/p2meg/MichelSpectrum.h`
- 目的: ミシェル崩壊のエネルギースペクトルに比例する角度積分後の非正規化 shape（全体定数は省略）を返します。

- シグネチャ
  - `double Michel_dShape_dE(double Ee, const MichelParams& mp = kMichelSM, const ParticleMasses& ms = kMassesPDG);`

- 入力:
  - `Ee`: 陽電子エネルギー（MeV）
  - `mp`: ミシェルパラメータ（`rho, eta, xi, delta`）。省略時は `kMichelSM`
  - `ms`: 粒子質量（`m_mu, m_e`）。省略時は `kMassesPDG`

- 出力:
  - 戻り値: 角度積分後の非正規化 shape 値（任意単位）。`x` が (0,1) 外では 0 を返します。

- 使用例:
```cpp
#include <iostream>
#include "p2meg/MichelSpectrum.h"
#include "p2meg/Constants.h"

int main() {
  const double Emax = Michel_Emax(kMassesPDG);
  const double Ee1 = 0.3 * Emax;
  const double Ee2 = 0.7 * Emax;

  std::cout << "dShape/dE(Ee1) = " << Michel_dShape_dE(Ee1) << "\n";
  std::cout << "dShape/dE(Ee2) = " << Michel_dShape_dE(Ee2) << "\n";
  return 0;
}
```

### RMD_d3B_dEe_dEg_dcos
- Header: `include/p2meg/RMDSpectrum.h`
- 目的: 停止ミューオン静止系における RMD の三重微分分岐比（shape の核）$d^3B/(dE_e\,dE_\gamma\,d\cos\theta_{e\gamma})$ をエネルギー変数（MeV）で返す。

- シグネチャ
  - `double RMD_d3B_dEe_dEg_dcos(double Ee_MeV, double Eg_MeV, double cosTheta, double d_min = 1e-6);`

- 入力:
  - `Ee_MeV`: 陽電子エネルギー $E_e$ [MeV]。
  - `Eg_MeV`: ガンマ線エネルギー $E_\gamma$ [MeV]。soft photon 発散があるため、呼び出し側で必ず $E_\gamma>E_{\gamma,\min}$ のカットを入れる。
  - `cosTheta`: $\cos\theta_{e\gamma}$（範囲は $[-1,1]$）。
  - `d_min`: 数値安定化のための下限（`RMD_d3B_dxdy_dcos` と同じ）。

- 出力:
  - 戻り値: 三重微分分岐比（核）$d^3B/(dE_e\,dE_\gamma\,d\cos\theta_{e\gamma})$（単位は MeV$^{-2}$）。運動学的に許されない領域では 0 を返す。

- 使用例:
```cpp
#include <iostream>
#include "p2meg/RMDSpectrum.h"
#include "p2meg/Constants.h"

int main() {
  const double Ee = 50.0;   // [MeV]
  const double Eg = 40.0;   // [MeV]（この点は角度の許容範囲が狭い）
  const double c  = -0.99;  // ほぼ反平行

  const double wE = RMD_d3B_dEe_dEg_dcos(Ee, Eg, c);
  std::cout << "d3B/dEe dEg dcos = " << wE << " [MeV^-2]" << std::endl;

  // 物理的に許されない領域では 0 が返る
  const double wE_bad = RMD_d3B_dEe_dEg_dcos(Ee, Eg, -0.5);
  std::cout << "bad point = " << wE_bad << " [MeV^-2]" << std::endl;

  return 0;
}
```

### RMD_d6B_dEe_dEg_dOmegae_dOmegag
- Header: `include/p2meg/RMDSpectrum.h`
- 目的: 停止ミューオン静止系における RMD の完全な微分分岐比（偏極込み）に対応する核 $d^6B/(dE_e\,dE_\gamma\,d\Omega_e\,d\Omega_\gamma)$ をエネルギー変数（MeV）で返す。

- シグネチャ
  - `double RMD_d6B_dEe_dEg_dOmegae_dOmegag(double Ee_MeV, double Eg_MeV, double cosThetaEG, double cosThetaE, double cosThetaG, double Pmu, double d_min = 1e-6);`

- 入力:
  - `Ee_MeV`: 陽電子エネルギー $E_e$ [MeV]。
  - `Eg_MeV`: ガンマ線エネルギー $E_\gamma$ [MeV]。soft photon 発散があるため、呼び出し側で必ず $E_\gamma>E_{\gamma,\min}$ のカットを入れる。
  - `cosThetaEG`: $\cos\theta_{e\gamma}=\hat{p}_e\cdot\hat{k}$（範囲は $[-1,1]$）。
  - `cosThetaE`: $\cos\theta_e=\hat{P}\cdot\hat{p}_e$（範囲は $[-1,1]$）。
  - `cosThetaG`: $\cos\theta_\gamma=\hat{P}\cdot\hat{k}$（範囲は $[-1,1]$）。
  - `Pmu`: 偏極度（スカラー、符号込み）。
  - `d_min`: 数値安定化のための下限（`RMD_d6B_dxdy_dOmegae_dOmegag` と同じ）。

- 出力:
  - 戻り値: 完全式（偏極込み）の核 $d^6B/(dE_e\,dE_\gamma\,d\Omega_e\,d\Omega_\gamma)$（単位は MeV$^{-2}$）。運動学的に許されない領域では 0 を返す。

- 使用例:
```cpp
#include <iostream>
#include "p2meg/RMDSpectrum.h"
#include "p2meg/Constants.h"

int main() {
  const double Ee = 50.0;   // [MeV]
  const double Eg = 10.0;   // [MeV]

  const double cosEG = -0.8; // p̂e · k̂
  const double cosE  = -0.2; // P̂ · p̂e
  const double cosG  =  0.6; // P̂ · k̂
  const double Pmu   = -1.0; // 理想偏極（例）

  const double wE = RMD_d6B_dEe_dEg_dOmegae_dOmegag(Ee, Eg, cosEG, cosE, cosG, Pmu);
  std::cout << "d6B/dEe dEg dOmegae dOmegag = " << wE << " [MeV^-2]" << std::endl;
  return 0;
}
```

### MakeRMDGridPdf
- Header: `include/p2meg/MakeRMDGridPdf.h`
- 目的: RMD 理論式（`RMD_d3B_dEe_dEg_dcos`）と検出器分解能（独立ガウシアン）を用いて、畳み込み済みの 3D 格子 PDF（Ee, Eg, θ）を生成し、ROOT ファイルに保存する。

- シグネチャ
  - `int MakeRMDGridPdf(const char* out_filepath, const char* key);`

- 入力:
  - `out_filepath`: 出力先 ROOT ファイルパス（例：`"data/pdf_cache/rmd_grid.root"`）。
  - `key`: ROOT ファイル中に保存するオブジェクトのキー名（例：`"rmd_grid"`）。

- 出力:
  - 戻り値: 成功時 0、失敗時は非0を返す。成功時、指定ファイルに 3D 格子 PDF（`key`）とメタ情報（`<key>_meta`）などを保存する。

- 使用例:
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

### RMDGridPdf_Load
- Header: `include/p2meg/RMDGridPdf.h`
- 目的: オフラインで生成した RMD 3D 格子 PDF（Ee, Eg, θ）を ROOT ファイルから読み込み、`RMDGridPdf(...)` で評価できる状態に初期化する。

- シグネチャ
  - `bool RMDGridPdf_Load(const char* filepath, const char* key);`

- 入力:
  - `filepath`: 入力 ROOT ファイルパス（例：`"data/pdf_cache/rmd_grid.root"`）。
  - `key`: ROOT ファイル中の格子 PDF のキー名（例：`"rmd_grid"`）。

- 出力:
  - 戻り値: 読み込みと初期化に成功したら `true`、失敗したら `false` を返す。

- 使用例:
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
  std::cout << "loaded" << std::endl;
  return 0;
}
```

### RMDGridPdf_IsLoaded
- Header: `include/p2meg/RMDGridPdf.h`
- 目的: RMD 格子 PDF がロード済みかどうかを返す（解析コード側の安全チェック用）。

- シグネチャ
  - `bool RMDGridPdf_IsLoaded();`

- 入力:
  - （なし）

- 出力:
  - 戻り値: ロード済みなら `true`、未ロードなら `false` を返す。

- 使用例:
```cpp
#include <iostream>
#include "p2meg/RMDGridPdf.h"

int main() {
  std::cout << "loaded? " << RMDGridPdf_IsLoaded() << std::endl;
  return 0;
}
```

### RMDGridPdf
- Header: `include/p2meg/RMDGridPdf.h`
- 目的: 観測値 (Ee, Eg, t, θ) に対して RMD の PDF 値を返す。3D 格子から p3(Ee,Eg,θ) を補間し、時間因子 p(t)（窓内正規化ガウシアン）を解析的に掛けて最終 PDF を評価する。

- シグネチャ
  - `double RMDGridPdf(double Ee, double Eg, double t, double theta);`

- 入力:
  - `Ee`: 陽電子エネルギー Ee [MeV]。
  - `Eg`: ガンマ線エネルギー Eg [MeV]。
  - `t`: 到達時間差 Δt [ns]。
  - `theta`: e と γ のなす角 θ [rad]（範囲は 0<θ<π を想定）。

- 出力:
  - 戻り値: 解析窓内なら PDF 値（密度）を返す。解析窓外、未ロード、数値的に不正な場合は 0 を返す。

- 使用例:
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

  const double Ee = 50.0;    // [MeV]
  const double Eg = 48.0;    // [MeV]
  const double t  = 0.1;     // [ns]
  const double th = 3.05;    // [rad]（ほぼ反平行）

  const double p = RMDGridPdf(Ee, Eg, t, th);
  std::cout << "RMDGridPdf = " << p << std::endl;
  return 0;
}
```

### SignalPdf
- Header: `include/p2meg/SignalPdf.h`
- 目的: 停止ミューオンの信号（μ+→e+γ）に対する解析的な 4D PDF を返す。Ee, Eg, t は解析窓内で正規化したトランケート正規分布、角度は δ=π−θ（δ≥0）で折り返した half-normal を解析窓に対応する δ 範囲で正規化して用いる。

- シグネチャ
  - `double SignalPdf(double Ee, double Eg, double t, double theta, const AnalysisWindow4D& win, const DetectorResolutionConst& res, const ParticleMasses& ms = kMassesPDG);`

- 入力:
  - `Ee`: 陽電子エネルギー Ee [MeV]（解析窓 `win.Ee_min..win.Ee_max` を想定）
  - `Eg`: ガンマ線エネルギー Eg [MeV]（解析窓 `win.Eg_min..win.Eg_max` を想定）
  - `t`: 到達時間差 Δt [ns]（解析窓 `win.t_min..win.t_max` を想定）
  - `theta`: e と γ のなす角 θ [rad]（0≤θ≤π を想定。角度PDFは δ=π−θ を half-normal にして用いる）
  - `win`: 解析窓（Ee, Eg, t, theta の各範囲）
  - `res`: 分解能パラメータ（`sigma_Ee`, `sigma_Eg`, `sigma_t`, `sigma_theta`, `t_mean`）
  - `ms`: 粒子質量（`m_mu`, `m_e`）。省略時は `kMassesPDG`（信号真値 Ee0=Eg0=m_mu/2 に使用）

- 出力:
  - 戻り値: 解析窓内での PDF 密度 p(Ee,Eg,t,theta) を返す。解析窓外、theta が [0,π] 外、分解能が不正、正規化定数が不正などの場合は 0 を返す。

- 使用例:
```cpp
#include <iostream>
#include "p2meg/SignalPdf.h"
#include "p2meg/AnalysisWindow.h"
#include "p2meg/DetectorResolution.h"
#include "p2meg/Constants.h"

int main() {
  const double Ee0 = 0.5 * kMassesPDG.m_mu;
  const double Eg0 = 0.5 * kMassesPDG.m_mu;
  const double t0  = detres.t_mean;
  const double th0 = pi;

  const double p0  = SignalPdf(Ee0, Eg0, t0, th0, analysis_window, detres);
  const double pth = SignalPdf(Ee0, Eg0, t0, pi - detres.sigma_theta, analysis_window, detres);

  std::cout << "SignalPdf(center) = " << p0 << "\n";
  std::cout << "SignalPdf(theta-1sigma) = " << pth << "\n";
  return 0;
}
```
