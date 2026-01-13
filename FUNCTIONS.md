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


### MakeRMDGridPdf
- Header: `include/p2meg/MakeRMDGridPdf.h`
- 目的: 停止ミューオン静止系における RMD の偏極込み核 `RMD_d6B_dEe_dEg_dOmegae_dOmegag` と検出器分解能（Ee, Eg の独立ガウシアン）を用いて、角度離散化（`N_theta`）込みの 4D 格子 PDF（Ee, Eg, cos_detector_e, cos_detector_g）を 2枝（cosΔφ=+1/-1）で生成し、ROOT ファイルに保存します（時間 t は評価側で解析的に掛ける）。

- シグネチャ
  - `int MakeRMDGridPdf(const char* out_filepath, const char* key);`

- 入力:
  - `out_filepath`: 出力先 ROOT ファイルパス（例：`"data/pdf_cache/rmd_grid.root"`）。
  - `key`: ROOT ファイル中に保存する格子 PDF のキー名（例：`"rmd_grid"`）。

- 出力:
  - 戻り値: 成功時 0、失敗時は非0を返す。成功時、指定ファイルに 4D 格子 PDF（`key+"_p"`, `key+"_m"`）とメタ情報（`<key>_meta`, `<key>_N_theta` など）を保存する。


### RMDGridPdf_Load
- Header: `include/p2meg/RMDGridPdf.h`
- 目的: オフラインで生成した RMD 4D 格子 PDF（Ee, Eg, cos_detector_e, cos_detector_g）を ROOT ファイルから読み込み、`RMDGridPdf(...)` で評価できる状態に初期化します。内部で `key+"_p"`（cosΔφ=+1）と `key+"_m"`（cosΔφ=-1）の2枝を同時にロードします。

- シグネチャ
  - `bool RMDGridPdf_Load(const char* filepath, const char* key);`

- 入力:
  - `filepath`: 入力 ROOT ファイルパス（例：`"data/pdf_cache/rmd_grid.root"`）
  - `key`: 格子 PDF のベースキー名（例：`"rmd_grid"`）。実際には `key+"_p"`, `key+"_m"` を読み込みます。

- 出力:
  - 戻り値: 2枝（`_p`, `_m`）のロードと内部クローン生成に成功したら `true`、失敗したら `false` を返します。

### RMDGridPdf_IsLoaded
- Header: `include/p2meg/RMDGridPdf.h`
- 目的: RMD 格子 PDF がロード済みかどうか（2枝とも揃っているか）を返します。解析コード側の安全チェック用です。

- シグネチャ
  - `bool RMDGridPdf_IsLoaded();`

- 入力:
  - （なし）

- 出力:
  - 戻り値: ロード済み（2枝とも利用可能）なら `true`、未ロードなら `false` を返します。



### RMDGridPdf
- Header: `include/p2meg/RMDGridPdf.h`
- 目的: 観測値 (Ee, Eg, t, theta, cos_detector_e, cos_detector_g) に対して RMD の PDF 値を返します。ROOT からロードした 4D 格子（Ee,Eg,cos_e,cos_g）を用いて評価し、時間因子 p(t)（窓内正規化ガウシアン）を解析的に掛けて最終 PDF を計算します。`theta` は生データ由来の離散化値を用い、cosΔφ=+1/-1 のどちらの枝かをハード判定して適切な格子（`_p` または `_m`）を選びます。

- シグネチャ
  - `double RMDGridPdf(double Ee, double Eg, double t, double theta, double cos_detector_e, double cos_detector_g);`

- 入力:
  - `Ee`: 陽電子エネルギー Ee [MeV]（解析窓 `analysis_window.Ee_min..Ee_max` を想定）
  - `Eg`: ガンマ線エネルギー Eg [MeV]（解析窓 `analysis_window.Eg_min..Eg_max` を想定）
  - `t`: 到達時間差 Δt [ns]（解析窓 `analysis_window.t_min..t_max` を想定）
  - `theta`: e と γ のなす角 θ [rad]（生データの θ を `detres.N_theta` 格子に最近傍丸めした値を想定）
  - `cos_detector_e`: 偏極軸と e 側検出器代表方向の内積（無次元、[-1,1]）。内部で最近傍の格子点に落として評価します。
  - `cos_detector_g`: 偏極軸と γ 側検出器代表方向の内積（無次元、[-1,1]）。内部で最近傍の格子点に落として評価します。

- 出力:
  - 戻り値: 解析窓内なら PDF 密度 p(Ee,Eg,t,cos_e,cos_g) を返します。窓外、未ロード、不正入力、格子評価が不正な場合は 0 を返します。


### SignalPdf
- Header: `include/p2meg/SignalPdf.h`
- 目的: 停止ミューオンの信号（μ+→e+γ）に対する解析的な 4D PDF を返す。Ee, Eg, t は解析窓内で正規化したトランケート正規分布、角度は `N_theta` 格子に丸めた離散角で扱い、理想化により θ=π のみに重みを持たせる。

- シグネチャ
  - `double SignalPdf(double Ee, double Eg, double t, double theta, const AnalysisWindow4D& win, const DetectorResolutionConst& res, const ParticleMasses& ms = kMassesPDG);`

- 入力:
  - `Ee`: 陽電子エネルギー Ee [MeV]（解析窓 `win.Ee_min..win.Ee_max` を想定）
  - `Eg`: ガンマ線エネルギー Eg [MeV]（解析窓 `win.Eg_min..win.Eg_max` を想定）
  - `t`: 到達時間差 Δt [ns]（解析窓 `win.t_min..win.t_max` を想定）
  - `theta`: e と γ のなす角 θ [rad]（0≤θ≤π を想定。内部で `N_theta` 格子に最近傍丸めし、θ=π のみに重み）
  - `win`: 解析窓（Ee, Eg, t, theta の各範囲）
  - `res`: 分解能パラメータ（`sigma_Ee`, `sigma_Eg`, `sigma_t`, `N_theta`, `t_mean`）
  - `ms`: 粒子質量（`m_mu`, `m_e`）。省略時は `kMassesPDG`（信号真値 Ee0=Eg0=m_mu/2 に使用）

- 出力:
  - 戻り値: 解析窓内での PDF 密度 p(Ee,Eg,t,theta) を返す。解析窓外、theta が [0,π] 外、分解能が不正、正規化定数が不正などの場合は 0 を返す。


### Event
- Header: `include/p2meg/Event.h`
- 目的: 尤度計算・PDF評価・入出力を疎結合にするための、1事象の観測量をまとめたデータ構造を提供します。

- シグネチャ
```cpp
struct Event {
  double Ee;
  double Eg;
  double t;
  double theta;
  double cos_detector_e;
  double cos_detector_g;
};
```

- 入力:
  - `Ee`: 陽電子エネルギー $E_e$ [MeV]
  - `Eg`: ガンマ線エネルギー $E_\gamma$ [MeV]
  - `t`: 到達時間差 $\Delta t$ [ns]
  - `theta`: e と $\gamma$ のなす角 $\theta$ [rad]
  - `cos_detector_e`: ミューオン偏極軸 $\hat{P}$ と e+ 側検出器代表方向 $\hat{d}_e$ の内積 $\hat{P}\cdot\hat{d}_e$（無次元、[-1,1]）
  - `cos_detector_g`: ミューオン偏極軸 $\hat{P}$ と $\gamma$ 側検出器代表方向 $\hat{d}_\gamma$ の内積 $\hat{P}\cdot\hat{d}_\gamma$（無次元、[-1,1]）

- 出力:
  - 戻り値: （なし）

### PdfComponent
- Header: `include/p2meg/Likelihood.h`
- 目的: 拡張尤度の混合モデルにおける1成分（例: sig, rmd, acc）の PDF 評価関数とその設定（ctx）をまとめた構造体を提供します。

- シグネチャ
```cpp
typedef double (*PdfEval)(const Event& ev, const void* ctx);

struct PdfComponent {
  const char* name;
  PdfEval eval;
  const void* ctx;
};
```

- 入力:
  - `name`: 成分名（例: `"sig"`, `"rmd"`。デバッグ用）
  - `eval`: PDF評価関数（`p_k(ev)` を返す関数ポインタ）
  - `ctx`: `eval` に渡す任意の設定ポインタ（不要なら `nullptr`）

- 出力:
  - 戻り値: （なし）

### ConstraintNLL
- Header: `include/p2meg/Likelihood.h`
- 目的: 事象数パラメータ（`yields`）に対して NLL に加算する制約項を返します。制約が不要な場合は 0 を返す実装として使用します。

- シグネチャ
```cpp
double ConstraintNLL(const std::vector<double>& yields);
```

- 入力:
  - `yields`: 解析窓内の期待事象数の配列（`{N_sig, N_rmd, (N_acc, ...)}`）。並びは `components` と同順

- 出力:
  - 戻り値: 制約項として NLL に加算する値（制約なしなら 0）

### NLL
- Header: `include/p2meg/Likelihood.h`
- 目的: 拡張尤度に基づく負の対数尤度（NLL）を計算します。PDF成分 `p_k(x)`（解析窓内で正規化済み）と期待事象数 `N_k`（`yields`）を用い、ポアソン項と混合項、および `ConstraintNLL` を加算します。

- シグネチャ
```cpp
double NLL(
  const std::vector<Event>& events,
  const std::vector<PdfComponent>& components,
  const std::vector<double>& yields
);
```

- 入力:
  - `events`: 解析窓内のイベント配列（各要素は `Event`）
  - `components`: PDF成分の配列（各要素は `PdfComponent`）
  - `yields`: 期待事象数の配列（`{N_sig, N_rmd, (N_acc, ...)}`）。`components` と同順

- 出力:
  - 戻り値: NLL 値（`(Σ_k N_k) - Σ_i log(Σ_k N_k p_k(x_i)) + ConstraintNLL(yields)` に対応）

### SignalPdfContext
- Header: `include/p2meg/PdfWrappers.h`
- 目的: `SignalPdf` を `PdfComponent` から評価できるようにするための設定（解析窓・分解能・質量）をまとめたコンテキストを提供します。

- シグネチャ
```cpp
struct SignalPdfContext {
  AnalysisWindow4D win;
  DetectorResolutionConst res;
  ParticleMasses ms;
};
```

- 入力:
  - `win`: 解析窓（`AnalysisWindow4D`）
  - `res`: 分解能パラメータ（`DetectorResolutionConst`）
  - `ms`: 粒子質量（`ParticleMasses`。信号真値 $E_{e0}=E_{\gamma0}=m_\mu/2$ に使用）

- 出力:
  - 戻り値: （なし）

### SignalPdfEval
- Header: `include/p2meg/PdfWrappers.h`
- 目的: `Event` を入力として `SignalPdf` を評価し、信号成分の PDF 密度を返します（`PdfEval` 互換）。

- シグネチャ
```cpp
double SignalPdfEval(const Event& ev, const void* ctx);
```

- 入力:
  - `ev`: 観測イベント（`Event`）
  - `ctx`: `SignalPdfContext` へのポインタ（`win,res,ms` を保持）

- 出力:
  - 戻り値: 信号 PDF 密度 $p_{\mathrm{sig}}(x)$（解析窓内で正規化済み）。窓外・不正入力などは 0

### MakeSignalComponent
- Header: `include/p2meg/PdfWrappers.h`
- 目的: 信号 PDF を尤度計算で扱える `PdfComponent` として生成します。

- シグネチャ
```cpp
PdfComponent MakeSignalComponent(const SignalPdfContext* ctx);
```

- 入力:
  - `ctx`: `SignalPdfContext` へのポインタ（呼び出し側で生存管理）

- 出力:
  - 戻り値: 信号成分の `PdfComponent`（`name="sig"`, `eval=&SignalPdfEval`, `ctx=ctx`）

### RMDGridPdfEval
- Header: `include/p2meg/PdfWrappers.h`
- 目的: `Event` を入力として `RMDGridPdf` を評価し、RMD成分の PDF 密度を返します（`PdfEval` 互換）。`Event` の `(Ee, Eg, t, theta)` に加えて `(cos_detector_e, cos_detector_g)` を渡し、角度離散化込みの評価を行います。

- シグネチャ
```cpp
double RMDGridPdfEval(const Event& ev, const void* ctx);
```

- 入力:
  - `ev`: 観測イベント（`Event`）。`ev.Ee`, `ev.Eg`, `ev.t`, `ev.theta`, `ev.cos_detector_e`, `ev.cos_detector_g` を使用します。
  - `ctx`: 未使用（`nullptr` を想定）。

- 出力:
  - 戻り値: RMD PDF 密度 $p_{\mathrm{rmd}}(x)$（解析窓内で正規化済み）。窓外・未ロード・不正入力などは 0 を返します。

### MakeRMDComponent
- Header: `include/p2meg/PdfWrappers.h`
- 目的: RMD PDF を尤度計算で扱える `PdfComponent` として生成します。

- シグネチャ
```cpp
PdfComponent MakeRMDComponent();
```

- 入力:
  - （なし）

- 出力:
  - 戻り値: RMD成分の `PdfComponent`（`name="rmd"`, `eval=&RMDGridPdfEval`, `ctx=nullptr`）

### FitConfig
- Header: `include/p2meg/NLLFit.h`
- 目的: `FitNLL` に渡す最小化設定（初期値・反復回数・収束判定）をまとめた構造体を提供します。

- シグネチャ
```cpp
struct FitConfig {
  std::vector<double> start_yields;
  int max_calls;
  double tol;
};
```

- 入力:
  - `start_yields`: 期待事象数の初期値配列（`components` と同順）
  - `max_calls`: 最大評価回数（実装側で解釈）
  - `tol`: 収束判定の許容値（実装側で解釈）

- 出力:
  - 戻り値: （なし）

### FitResult
- Header: `include/p2meg/NLLFit.h`
- 目的: `FitNLL` の結果（推定値・誤差・最小NLL・ステータス）を格納します。

- シグネチャ
```cpp
struct FitResult {
  int status;
  std::vector<double> yields_hat;
  std::vector<double> yields_err;
  double nll_min;
};
```

- 入力:
  - `status`: フィットの成否コード（実装側で定義）
  - `yields_hat`: 推定された期待事象数（`components` と同順）
  - `yields_err`: 推定誤差（取れる場合のみ。取れない場合は空でもよい）
  - `nll_min`: 最小 NLL 値

- 出力:
  - 戻り値: （なし）

### FitNLL
- Header: `include/p2meg/NLLFit.h`
- 目的: `NLL(...)` を最小化して、期待事象数（`yields`）の最尤推定値を求めます。尤度計算と最小化過程を分離するための入口関数です。

- シグネチャ
```cpp
FitResult FitNLL(
  const std::vector<Event>& events,
  const std::vector<PdfComponent>& components,
  const FitConfig& cfg
);
```

- 入力:
  - `events`: 解析窓内のイベント配列（`Event`）
  - `components`: PDF成分の配列（`PdfComponent`）。並びが `yields` の意味を決める
  - `cfg`: 最小化設定（`FitConfig`）

- 出力:
  - 戻り値: フィット結果（`FitResult`）
