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

### Michel_Emax
- Header: `include/p2meg/MichelSpectrum.h`
- 目的: 粒子質量からミシェル崩壊の端点エネルギー $ E_{\max} $ を返します（実装では $ E_{\max} = m_\mu/2 $）。

- シグネチャ
  - `double Michel_Emax(const ParticleMasses& ms);`

- 入力:
  - `ms`: 粒子質量（`ms.m_mu` と `ms.m_e` を持つ。単位は MeV）

- 出力:
  - 戻り値: 端点エネルギー（MeV）。実装では `ms.m_mu/2.0` を返します。

- 使用例:
```cpp
#include <iostream>
#include "p2meg/MichelSpectrum.h"
#include "p2meg/Constants.h"

int main() {
  const ParticleMasses ms = kMassesPDG;
  const double Emax = Michel_Emax(ms);
  std::cout << "Emax = " << Emax << " MeV\n";
  return 0;
}
```

### Michel_x
- Header: `include/p2meg/MichelSpectrum.h`
- 目的: 陽電子エネルギーから無次元変数 $ x = E_e/E_{\max} $ を計算します。

- シグネチャ
  - `double Michel_x(double Ee, const ParticleMasses& ms);`

- 入力:
  - `Ee`: 陽電子エネルギー（MeV）
  - `ms`: 粒子質量（端点計算に `ms.m_mu` を使用。単位は MeV）

- 出力:
  - 戻り値: 無次元変数 `x`（`Ee / Michel_Emax(ms)`）。

- 使用例:
```cpp
#include <iostream>
#include "p2meg/MichelSpectrum.h"
#include "p2meg/Constants.h"

int main() {
  const ParticleMasses ms = kMassesPDG;
  const double Emax = Michel_Emax(ms);
  const double Ee = 0.7 * Emax;
  std::cout << "x = " << Michel_x(Ee, ms) << "\n";
  return 0;
}
```

### Michel_x0
- Header: `include/p2meg/MichelSpectrum.h`
- 目的: ミシェル分布の $ \eta $ 項に現れる係数 $ x_0 = m_e/E_{\max} $ を計算します。

- シグネチャ
  - `double Michel_x0(const ParticleMasses& ms);`

- 入力:
  - `ms`: 粒子質量（`ms.m_e` と `ms.m_mu` を使用。単位は MeV）

- 出力:
  - 戻り値: 無次元係数 `x0`（`ms.m_e / Michel_Emax(ms)`）。

- 使用例:
```cpp
#include <iostream>
#include "p2meg/MichelSpectrum.h"
#include "p2meg/Constants.h"

int main() {
  const ParticleMasses ms = kMassesPDG;
  std::cout << "x0 = " << Michel_x0(ms) << "\n";
  return 0;
}
```

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

### RMD_d3B_dxdy_dcos
- Header: `include/p2meg/RMDSpectrum.h`
- 目的: 停止ミューオン静止系における RMD（$\mu^+\to e^+\nu\bar{\nu}\gamma$）の三重微分分岐比（shape の核）$d^3B/(dx\,dy\,d\cos\theta_{e\gamma})$ を返す。

- シグネチャ
  - `double RMD_d3B_dxdy_dcos(double x, double y, double cosTheta, double d_min = 1e-6);`

- 入力:
  - `x`: $x=2E_e/m_\mu$（無次元）。
  - `y`: $y=2E_\gamma/m_\mu$（無次元）。soft photon 発散があるため、呼び出し側で必ず $y>y_{\min}$（$E_\gamma>E_{\gamma,\min}$）のカットを入れる。
  - `cosTheta`: $\cos\theta_{e\gamma}$（$e$ と $\gamma$ のなす角の余弦、範囲は $[-1,1]$）。
  - `d_min`: 数値安定化のための下限。$d=1-\beta\cos\theta_{e\gamma}$ が `d_min` 未満なら `d_min` に置き換える（物理カットではない）。

- 出力:
  - 戻り値: 三重微分分岐比（核）$d^3B/(dx\,dy\,d\cos\theta_{e\gamma})$（無次元）。運動学的に許されない領域では 0 を返す。

- 使用例:
```cpp
#include <iostream>
#include "p2meg/RMDSpectrum.h"
#include "p2meg/Constants.h"

int main() {
  const double mmu = kMassesPDG.m_mu;

  const double Ee = 50.0;   // [MeV]
  const double Eg = 10.0;   // [MeV]（必ず下限を入れる）
  const double c  = -0.8;   // cos(theta_eγ)

  const double x = 2.0 * Ee / mmu;
  const double y = 2.0 * Eg / mmu;

  const double w = RMD_d3B_dxdy_dcos(x, y, c);
  std::cout << "d3B/dx dy dcos = " << w << std::endl;
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
