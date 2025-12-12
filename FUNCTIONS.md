# Functions

このファイルは **p2meg_analysis の関数・クラスの「使い方」**をまとめたものです。

# Contribution

新しい関数・クラスを追加した場合は、この`FUNCTIONS.md` に「使い方セクション」と Index の1行を追加してください。

## AIに「関数の使い方」を書かせるプロンプト

関数をAIで作った場合は以下をそのままAIに貼り付けてください。

---

### プロンプト本文

今実装した関数について、`p2meg_analysis/FUNCTIONS.md` に追記する「使い方」だけを書いてください。実装の変更提案や改善案は不要です。  

返答は コードブロック1つだけにしてください。コードブロックの外に文章を一切書かないでください（前置き・説明・箇条書き・空行の追加も禁止）。

そのコードブロックは 言語をmdにし、開始と終了は4つのバッククォートで囲ってください（````md 〜 ````）。

コードブロックの中身は、私が `FUNCTIONS.md` にそのまま貼り付けられる **Markdown断片のみ**にしてください。

また、使用例はコードブロックの中で通常どおり
```cpp
...
```
の fenced code block を使ってください。

次の項目だけをこの順番で書いてください：
1. 関数名（またはクラス名）
2. ヘッダー名（`include/p2meg/...` のパス）
3. 目的
4. 入力（引数ごと：意味を1行ずつ）
5. 出力（戻り値：意味を1〜2文）
6. 使用例（簡単なC++例。10〜20行程度で、コピペで動く形）

**出力フォーマットは必ず次に従ってください：**

```md
### <関数名>
- Header: `include/p2meg/<...>.h`
- 目的: ...

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

# Function list

### Michel_Emax
- Header: `include/p2meg/MichelSpectrum.h`
- 目的: 粒子質量からミシェル崩壊の端点エネルギー [$ E_{\max}] を返します（実装では [$ E_{\max} = m_\mu/2]）。

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
- 目的: 陽電子エネルギーから無次元変数 [$ x = E_e/E_{\max}] を計算します。

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
- 目的: ミシェル分布の [$ \eta] 項に現れる係数 [$ x_0 = m_e/E_{\max}] を計算します。

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
- 目的: ミシェル崩壊の二重微分分布に比例する「非正規化 shape」[$ d^2\Gamma/(dE_e\,d\cos\theta)] を返します（全体定数は省略）。

- 入力:
  - `Ee`: 陽電子エネルギー（MeV）
  - `costh`: [$ \cos\theta]（[$ \theta] は偏極軸と陽電子運動量のなす角）
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

  const double w_p = Michel_d2Shape_dE_dCosTheta(Ee, +1.0, 1.0);
  const double w_0 = Michel_d2Shape_dE_dCosTheta(Ee,  0.0, 1.0);
  const double w_m = Michel_d2Shape_dE_dCosTheta(Ee, -1.0, 1.0);

  std::cout << "w(costh=+1) = " << w_p << "\n";
  std::cout << "w(costh= 0) = " << w_0 << "\n";
  std::cout << "w(costh=-1) = " << w_m << "\n";
  return 0;
}
```

### Michel_dShape_dE
- Header: `include/p2meg/MichelSpectrum.h`
- 目的: ミシェル崩壊のエネルギースペクトルに比例する「角度積分後の非正規化 shape」[$ d\Gamma/dE_e] を返します（全体定数は省略）。

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

  const double w1 = Michel_dShape_dE(Ee1);
  const double w2 = Michel_dShape_dE(Ee2);

  std::cout << "dShape/dE(Ee1) = " << w1 << "\n";
  std::cout << "dShape/dE(Ee2) = " << w2 << "\n";
  return 0;
}
```
