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

### Math_IsFinite
- Header: `include/p2meg/MathUtils.h`
- 目的: 数値が有限かどうかを判定します（NaN/inf を除外するガード）。

- シグネチャ
  - `static inline bool Math_IsFinite(double x);`

- 入力:
  - `x`: 判定対象の実数。

- 出力:
  - 戻り値: 有限なら `true`、NaN または ±inf なら `false` を返します。

### Math_Clamp
- Header: `include/p2meg/MathUtils.h`
- 目的: 値を指定範囲 [lo, hi] にクリップします（物理カットではない数値ガード）。

- シグネチャ
  - `static inline double Math_Clamp(double x, double lo, double hi);`

- 入力:
  - `x`: 入力値。
  - `lo`: 下限。
  - `hi`: 上限。

- 出力:
  - 戻り値: `x` を [lo, hi] に収めた値を返します。

### Math_GetNTheta
- Header: `include/p2meg/MathUtils.h`
- 目的: 分解能設定 `N_theta` を安全に int として取得し、最小値 1 を保証します。

- シグネチャ
  - `static inline int Math_GetNTheta(const DetectorResolutionConst& res);`

- 入力:
  - `res`: 検出器分解能設定（`DetectorResolutionConst`）。

- 出力:
  - 戻り値: `res.N_theta` を丸めた整数値（下限 1）。不正値でも 1 を返します。

### Angle_ClipPhi0Pi
- Header: `include/p2meg/AngleUtils.h`
- 目的: 角度 φ を [0, π] にクリップします（不正入力は 0 に落とす）。

- シグネチャ
  - `static inline double Angle_ClipPhi0Pi(double phi);`

- 入力:
  - `phi`: 角度 φ [rad]。

- 出力:
  - 戻り値: φ を [0, π] に収めた値を返します。

### Angle_DiscretizePhi
- Header: `include/p2meg/AngleUtils.h`
- 目的: 角度 φ を [0, π] にクリップしたうえで、離散格子点 φ_i=iπ/N_theta に丸めます。

- シグネチャ
  - `static inline double Angle_DiscretizePhi(double phi, int N_theta);`

- 入力:
  - `phi`: 角度 φ [rad]。
  - `N_theta`: 角度分割数（`N_theta>=1` を想定）。

- 出力:
  - 戻り値: 離散格子点に丸めた φ を返します。不正値の場合は 0 を返します。

### Angle_ThetaFromPhiStrict
- Header: `include/p2meg/AngleUtils.h`
- 目的: φ_e, φ_g から e-γ の相対角 θ=|φ_e-φ_g| を求めます（範囲外は不正扱い）。

- シグネチャ
  - `static inline double Angle_ThetaFromPhiStrict(double phi_e, double phi_g);`

- 入力:
  - `phi_e`: e 側検出器の角度 φ_e [rad]（0..π を想定）。
  - `phi_g`: γ 側検出器の角度 φ_g [rad]（0..π を想定）。

- 出力:
  - 戻り値: 有効なら θ=|φ_e-φ_g| を返します。不正入力では -1 を返します。

### Angle_ThetaFromPhiClipped
- Header: `include/p2meg/AngleUtils.h`
- 目的: φ_e, φ_g を [0, π] に収めてから相対角 θ=|φ_e-φ_g| を求めます。

- シグネチャ
  - `static inline double Angle_ThetaFromPhiClipped(double phi_e, double phi_g);`

- 入力:
  - `phi_e`: e 側検出器の角度 φ_e [rad]。
  - `phi_g`: γ 側検出器の角度 φ_g [rad]。

- 出力:
  - 戻り値: クリップ後の θ を返します。不正入力では 0 を返します。

### AnalysisWindow_In4D
- Header: `include/p2meg/AnalysisWindowUtils.h`
- 目的: 解析窓 (Ee, Eg, t, theta) の範囲内かどうかを判定します。

- シグネチャ
  - `static inline bool AnalysisWindow_In4D(const AnalysisWindow4D& win, double Ee, double Eg, double t, double theta);`

- 入力:
  - `win`: 解析窓設定（`AnalysisWindow4D`）。
  - `Ee`: 陽電子エネルギー Ee [MeV]。
  - `Eg`: ガンマ線エネルギー Eg [MeV]。
  - `t`: 到達時間差 Δt [ns]。
  - `theta`: e-γ 相対角 θ [rad]。

- 出力:
  - 戻り値: 解析窓内なら `true`、外なら `false`。

### AnalysisWindow_In3D
- Header: `include/p2meg/AnalysisWindowUtils.h`
- 目的: 解析窓 (Ee, Eg, theta) の範囲内かどうかを判定します。

- シグネチャ
  - `static inline bool AnalysisWindow_In3D(const AnalysisWindow4D& win, double Ee, double Eg, double theta);`

- 入力:
  - `win`: 解析窓設定（`AnalysisWindow4D`）。
  - `Ee`: 陽電子エネルギー Ee [MeV]。
  - `Eg`: ガンマ線エネルギー Eg [MeV]。
  - `theta`: e-γ 相対角 θ [rad]。

- 出力:
  - 戻り値: 解析窓内なら `true`、外なら `false`。

### AnalysisWindow_InTime
- Header: `include/p2meg/AnalysisWindowUtils.h`
- 目的: 解析窓の時間範囲に入っているかどうかを判定します。

- シグネチャ
  - `static inline bool AnalysisWindow_InTime(const AnalysisWindow4D& win, double t);`

- 入力:
  - `win`: 解析窓設定（`AnalysisWindow4D`）。
  - `t`: 到達時間差 Δt [ns]。

- 出力:
  - 戻り値: 解析窓内なら `true`、外なら `false`。

### AnalysisWindow_InTimeSideband
- Header: `include/p2meg/AnalysisWindowUtils.h`
- 目的: 全時間範囲 [t_all_min, t_all_max] の中で、解析窓の時間範囲 [win.t_min, win.t_max] の外側（TSB: timing sideband）に入っているかを判定します。

- シグネチャ
  - `static inline bool AnalysisWindow_InTimeSideband(const AnalysisWindow4D& win, double t, double t_all_min, double t_all_max);`

- 入力:
  - `win`: 解析窓設定（`AnalysisWindow4D`。ここでは `t_min`, `t_max` を使用）。
  - `t`: 到達時間差 Δt [ns]。
  - `t_all_min`: 全時間範囲の下限 [ns]。
  - `t_all_max`: 全時間範囲の上限 [ns]。

- 出力:
  - 戻り値: 全時間範囲内かつ解析窓の外側（TSB）なら `true`、それ以外（全時間範囲外、解析窓内、不正入力など）は `false` を返します。

### AnalysisWindow_InTSB
- Header: `include/p2meg/AnalysisWindowUtils.h`
- 目的: (Ee, Eg, theta) が解析窓内で、かつ t が全時間範囲内で解析窓外（TSB）であるかを判定します（TSB イベント選別用）。

- シグネチャ
  - `static inline bool AnalysisWindow_InTSB(const AnalysisWindow4D& win, double Ee, double Eg, double theta, double t, double t_all_min, double t_all_max);`

- 入力:
  - `win`: 解析窓設定（`AnalysisWindow4D`）。
  - `Ee`: 陽電子エネルギー Ee [MeV]。
  - `Eg`: ガンマ線エネルギー Eg [MeV]。
  - `theta`: e-γ 相対角 θ [rad]。
  - `t`: 到達時間差 Δt [ns]。
  - `t_all_min`: 全時間範囲の下限 [ns]。
  - `t_all_max`: 全時間範囲の上限 [ns]。

- 出力:
  - 戻り値: (Ee,Eg,theta) が窓内かつ t が TSB なら `true`、それ以外は `false` を返します。

### AnalysisWindow_TimeSidebandWidth
- Header: `include/p2meg/AnalysisWindowUtils.h`
- 目的: 全時間範囲 [t_all_min, t_all_max] に対する TSB（解析窓外側）の合計幅（左 + 右）を返します（TSB→解析窓へのスケール係数計算などに使用）。

- シグネチャ
  - `static inline double AnalysisWindow_TimeSidebandWidth(const AnalysisWindow4D& win, double t_all_min, double t_all_max);`

- 入力:
  - `win`: 解析窓設定（`AnalysisWindow4D`。ここでは `t_min`, `t_max` を使用）。
  - `t_all_min`: 全時間範囲の下限 [ns]。
  - `t_all_max`: 全時間範囲の上限 [ns]。

- 出力:
  - 戻り値: TSB の合計幅 [ns] を返します。不正入力や幅が 0 の場合は 0 を返します。

### Hist_AxisBracketUniform
- Header: `include/p2meg/HistUtils.h`
- 目的: 等間隔ビン軸に対して、補間用の隣接ビンと係数 (i0, i1, f) を求めます。

- シグネチャ
  - `int Hist_AxisBracketUniform(const TAxis& ax, double x, int& i0, int& i1, double& f);`

- 入力:
  - `ax`: ROOT の軸（等間隔ビン前提）。
  - `x`: 対象の座標。
  - `i0`: 出力（下側のビン番号 1..n-1）。
  - `i1`: 出力（上側のビン番号 i0+1）。
  - `f`: 出力（補間係数、0..1）。

- 出力:
  - 戻り値: 成功時 0、失敗時は非0 を返します。

### Hist_InterpEeEg4
- Header: `include/p2meg/HistUtils.h`
- 目的: 4D THnD の Ee/Eg を 2D（4点）補間し、phi は固定ビンで評価します。

- シグネチャ
  - `double Hist_InterpEeEg4(const THnD& h, double Ee, double Eg, int bin_phi_e, int bin_phi_g);`

- 入力:
  - `h`: 4D ヒストグラム（Ee, Eg, phi_e, phi_g）。
  - `Ee`: 陽電子エネルギー Ee [MeV]。
  - `Eg`: ガンマ線エネルギー Eg [MeV]。
  - `bin_phi_e`: phi_e のビン番号。
  - `bin_phi_g`: phi_g のビン番号。

- 出力:
  - 戻り値: 補間値（不正/負値は 0 を返します）。

### Hist_SumAllBins4
- Header: `include/p2meg/HistUtils.h`
- 目的: 4D THnD の全ビン内容の総和を返します。

- シグネチャ
  - `double Hist_SumAllBins4(const THnD& h);`

- 入力:
  - `h`: 4D ヒストグラム。

- 出力:
  - 戻り値: 全ビンの総和を返します。


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
- 目的: 停止ミューオン静止系における RMD の偏極込み核 `RMD_d6B_dEe_dEg_dOmegae_dOmegag` と検出器分解能（`energy_response_shape_e/g` による Ee/Eg スメア）を用いて、角度離散化（`N_theta`）込みの 4D 格子 PDF（Ee, Eg, phi_detector_e, phi_detector_g）を生成し、ROOT ファイルに保存します（時間 t は評価側で解析的に掛ける）。真値サンプル窓は解析窓を energy_response の 0.1 倍点まで広げ、物理領域でクリップしたものを使います。

- シグネチャ
  - `int MakeRMDGridPdf(const char* out_filepath, const char* key);`
  - `int MakeRMDGridPdfWithTruthWindow(const char* out_filepath, const char* key, const AnalysisWindow4D& truth_win);`

- 入力:
  - `out_filepath`: 出力先 ROOT ファイルパス（例：`"data/pdf_cache/rmd_grid.root"`）。
  - `key`: ROOT ファイル中に保存する格子 PDF のキー名（例：`"rmd_grid"`）。
  - `truth_win`: 解析窓として扱う Ee/Eg を指定し、`energy_response_shape_e/g` の 0.1 倍点まで広げた真値サンプル窓を内部で作る（非物理領域は除外、Eg_min>0 が必要）。t, theta は未使用。単位は MeV。

- 出力:
  - 戻り値: 成功時 0、失敗時は非0を返す。成功時、指定ファイルに 4D 格子 PDF（`key`）とメタ情報（`<key>_meta`, `<key>_N_theta` など）を保存する。


### RMDGridPdf_Load
- Header: `include/p2meg/RMDGridPdf.h`
- 目的: オフラインで生成した RMD 4D 格子 PDF（Ee, Eg, phi_detector_e, phi_detector_g）を ROOT ファイルから読み込み、`RMDGridPdf(...)` で評価できる状態に初期化します。

- シグネチャ
  - `bool RMDGridPdf_Load(const char* filepath, const char* key);`

- 入力:
  - `filepath`: 入力 ROOT ファイルパス（例：`"data/pdf_cache/rmd_grid.root"`）
  - `key`: 格子 PDF のキー名（例：`"rmd_grid"`）。

- 出力:
  - 戻り値: ロードと内部クローン生成に成功したら `true`、失敗したら `false` を返します。

### RMDGridPdf_IsLoaded
- Header: `include/p2meg/RMDGridPdf.h`
- 目的: RMD 格子 PDF がロード済みかどうかを返します。解析コード側の安全チェック用です。

- シグネチャ
  - `bool RMDGridPdf_IsLoaded();`

- 入力:
  - （なし）

- 出力:
  - 戻り値: ロード済みなら `true`、未ロードなら `false` を返します。



### RMDGridPdf
- Header: `include/p2meg/RMDGridPdf.h`
- 目的: 観測値 (Ee, Eg, t, phi_detector_e, phi_detector_g) に対して RMD の PDF 値を返します。ROOT からロードした 4D 格子（Ee,Eg,phi_e,phi_g）を用いて評価し、時間因子 p(t)（窓内正規化ガウシアン）を解析的に掛けて最終 PDF を計算します。`theta_eg=|phi_e-phi_g|` を作って解析窓カットを行います。

- シグネチャ
  - `double RMDGridPdf(double Ee, double Eg, double t, double phi_detector_e, double phi_detector_g);`

- 入力:
  - `Ee`: 陽電子エネルギー Ee [MeV]（解析窓 `analysis_window.Ee_min..Ee_max` を想定）
  - `Eg`: ガンマ線エネルギー Eg [MeV]（解析窓 `analysis_window.Eg_min..Eg_max` を想定）
  - `t`: 到達時間差 Δt [ns]（解析窓 `analysis_window.t_min..t_max` を想定）
  - `phi_detector_e`: 偏極軸と e 側検出器代表方向の角度 φ_e [rad]（0..π を想定。評価時は 0..π にクリップ）
  - `phi_detector_g`: 偏極軸と γ 側検出器代表方向の角度 φ_g [rad]（0..π を想定。評価時は 0..π にクリップ）

- 出力:
  - 戻り値: 解析窓内なら PDF 密度 p(Ee,Eg,t,phi_e,phi_g) を返します。窓外、未ロード、不正入力、格子評価が不正な場合は 0 を返します。

### MakeACCGridPdf
- Header: `include/p2meg/MakeACCGridPdf.h`
- 目的: ACC (accidental) 成分の TSB（時間が全時間範囲内で解析窓外）イベントから、解析窓内の 4D 格子 PDF（Ee, Eg, phi_detector_e, phi_detector_g）を作成して ROOT に保存します。phi は N_theta の離散点に丸め、正規化は「Σ_{phi_e,phi_g} ∫ dEe dEg p4 = 1」となるよう Ee/Eg の bin 幅のみを測度に入れます。時間は評価側で解析窓内一様として掛けます。

- シグネチャ
  - `int MakeACCGridPdf(const std::vector<Event>& events, const char* out_filepath, const char* key);`

- 入力:
  - `events`: 入力イベント配列（各 `Event` に `Ee, Eg, t, phi_detector_e, phi_detector_g` を格納）
  - `out_filepath`: 出力先 ROOT ファイルパス（例：`"data/pdf_cache/acc_grid.root"`）
  - `key`: ROOT ファイル中に保存する格子 PDF のキー名（例：`"acc_grid"`）

- 出力:
  - 戻り値: 成功時 0、失敗時は非0を返す。成功時、指定ファイルに 4D 格子 PDF（`key`）とメタ情報（`<key>_meta`, `<key>_N_theta` など）を保存する。

### ACCGridPdf_Load
- Header: `include/p2meg/ACCGridPdf.h`
- 目的: オフラインで生成した ACC 4D 格子 PDF（Ee, Eg, phi_detector_e, phi_detector_g）を ROOT ファイルから読み込み、`ACCGridPdf(...)` で評価できる状態に初期化します。

- シグネチャ
  - `bool ACCGridPdf_Load(const char* filepath, const char* key);`

- 入力:
  - `filepath`: 入力 ROOT ファイルパス（例：`"data/pdf_cache/acc_grid.root"`）
  - `key`: 格子 PDF のキー名（例：`"acc_grid"`）。

- 出力:
  - 戻り値: ロードと内部クローン生成に成功したら `true`、失敗したら `false` を返します。

### ACCGridPdf_IsLoaded
- Header: `include/p2meg/ACCGridPdf.h`
- 目的: ACC 格子 PDF がロード済みかどうかを返します。解析コード側の安全チェック用です。

- シグネチャ
  - `bool ACCGridPdf_IsLoaded();`

- 入力:
  - （なし）

- 出力:
  - 戻り値: ロード済みなら `true`、未ロードなら `false` を返します。

### ACCGridPdf
- Header: `include/p2meg/ACCGridPdf.h`
- 目的: 観測値 (Ee, Eg, t, phi_detector_e, phi_detector_g) に対して ACC の PDF 値を返します。ROOT からロードした 4D 格子（Ee,Eg,phi_e,phi_g）を用いて評価し、時間因子は解析窓内一様として解析的に掛けます。`theta_eg=|phi_e-phi_g|` を作って解析窓カットを行います。

- シグネチャ
  - `double ACCGridPdf(double Ee, double Eg, double t, double phi_detector_e, double phi_detector_g);`

- 入力:
  - `Ee`: 陽電子エネルギー Ee [MeV]（解析窓 `analysis_window.Ee_min..Ee_max` を想定）
  - `Eg`: ガンマ線エネルギー Eg [MeV]（解析窓 `analysis_window.Eg_min..Eg_max` を想定）
  - `t`: 到達時間差 Δt [ns]（解析窓 `analysis_window.t_min..t_max` を想定）
  - `phi_detector_e`: 偏極軸と e 側検出器代表方向の角度 φ_e [rad]（0..π を想定。評価時は 0..π にクリップして離散化）
  - `phi_detector_g`: 偏極軸と γ 側検出器代表方向の角度 φ_g [rad]（0..π を想定。評価時は 0..π にクリップして離散化）

- 出力:
  - 戻り値: 解析窓内なら PDF 密度 p(Ee,Eg,t,phi_e,phi_g) を返します。窓外、未ロード、不正入力、格子評価が不正な場合は 0 を返します。


### SignalPdf
- Header: `include/p2meg/SignalPdf.h`
- 目的: 停止ミューオンの信号（μ+→e+γ）に対する解析的な 5D PDF を返す。Ee, Eg は `energy_response_shape_e/g` を解析窓内で正規化したエネルギー応答、t は解析窓内で正規化したトランケート正規分布、角度は `N_theta` 格子に丸めた離散角で扱い、理想化により θ=π のみに重みを持たせる。角度は `phi_detector_e/g`（0..π）から `theta_eg=|phi_e-phi_g|` を作って評価する（phi ベースの角度評価）。

- シグネチャ
  - `double SignalPdf(double Ee, double Eg, double t, double phi_detector_e, double phi_detector_g, const AnalysisWindow4D& win, const DetectorResolutionConst& res, const ParticleMasses& ms = kMassesPDG);`

- 入力:
  - `Ee`: 陽電子エネルギー Ee [MeV]（解析窓 `win.Ee_min..win.Ee_max` を想定）
  - `Eg`: ガンマ線エネルギー Eg [MeV]（解析窓 `win.Eg_min..win.Eg_max` を想定）
  - `t`: 到達時間差 Δt [ns]（解析窓 `win.t_min..win.t_max` を想定）
  - `phi_detector_e`: 偏極軸と e 側検出器代表方向の角度 φ_e [rad]（0..π を想定）
  - `phi_detector_g`: 偏極軸と γ 側検出器代表方向の角度 φ_g [rad]（0..π を想定）
  - `win`: 解析窓（Ee, Eg, t, theta の各範囲）
  - `res`: 分解能パラメータ（`sigma_t`, `N_theta`, `t_mean` を使用。Ee/Eg の応答 shape は `DetectorResolution.h` の `energy_response_shape_e/g` を参照）
  - `ms`: 粒子質量（`m_mu`, `m_e`）。省略時は `kMassesPDG`（信号真値 Ee0=Eg0=m_mu/2 に使用）

- 出力:
  - 戻り値: 解析窓内での PDF 密度 p(Ee,Eg,t,phi_e,phi_g) を返す。解析窓外（theta は phi から作った角度を格子に丸めた値で判定）、phi が [0,π] 外、分解能が不正、正規化定数が不正などの場合は 0 を返す。


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
- 目的: `Event` を入力として `RMDGridPdf` を評価し、RMD成分の PDF 密度を返します（`PdfEval` 互換）。`Event` の `(Ee, Eg, t, phi_detector_e, phi_detector_g)` を渡して評価します。

- シグネチャ
```cpp
double RMDGridPdfEval(const Event& ev, const void* ctx);
```

- 入力:
  - `ev`: 観測イベント（`Event`）。`ev.Ee`, `ev.Eg`, `ev.t`, `ev.phi_detector_e`, `ev.phi_detector_g` を使用します。
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

### ACCGridPdfEval
- Header: `include/p2meg/PdfWrappers.h`
- 目的: `Event` を入力として `ACCGridPdf` を評価し、ACC成分の PDF 密度を返します（`PdfEval` 互換）。`Event` の `(Ee, Eg, t, phi_detector_e, phi_detector_g)` を渡して評価します。

- シグネチャ
```cpp
double ACCGridPdfEval(const Event& ev, const void* ctx);
```

- 入力:
  - `ev`: 観測イベント（`Event`）。`ev.Ee`, `ev.Eg`, `ev.t`, `ev.phi_detector_e`, `ev.phi_detector_g` を使用します。
  - `ctx`: 未使用（`nullptr` を想定）。

- 出力:
  - 戻り値: ACC PDF 密度 $p_{\mathrm{acc}}(x)$（解析窓内で正規化済み）。窓外・未ロード・不正入力などは 0 を返します。

### MakeACCComponent
- Header: `include/p2meg/PdfWrappers.h`
- 目的: ACC PDF を尤度計算で扱える `PdfComponent` として生成します。

- シグネチャ
```cpp
PdfComponent MakeACCComponent();
```

- 入力:
  - （なし）

- 出力:
  - 戻り値: ACC成分の `PdfComponent`（`name="acc"`, `eval=&ACCGridPdfEval`, `ctx=nullptr`）

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

### FitNLLFixedSignal
- Header: `include/p2meg/UpperLimit.h`
- 目的: signal yield `N_sig` を固定し、他の背景 yield を profile fit で再最適化します。μ→eγ upper limit 計算で各仮説値 `s=N_sig` に対する条件付き最尤点を求めるために使います。

- シグネチャ
```cpp
FitResult FitNLLFixedSignal(
    const std::vector<Event>& events,
    const std::vector<PdfComponent>& components,
    const FitConfig& cfg,
    double N_sig_fixed
);
```

- 入力:
  - `events`: 解析窓内のイベント配列（`Event`）。
  - `components`: PDF成分の配列（`PdfComponent`）。先頭成分を signal (`N_sig`) とみなします。
  - `cfg`: 最小化設定（`FitConfig`）。`start_yields` は `components` と同じ長さが必要です。
  - `N_sig_fixed`: 固定する signal yield 仮説値 `s`。単位は期待事象数。

- 出力:
  - 戻り値: 条件付き最尤フィット結果（`FitResult`）。失敗時は `status!=0` を返します。

### GenerateToyDatasetFromModel
- Header: `include/p2meg/UpperLimit.h`
- 目的: 最終解析で使う PDF 成分から直接 toy pseudo-experiment を 1 本生成します。各成分の事象数には Poisson fluctuation を入れます。

- シグネチャ
```cpp
bool GenerateToyDatasetFromModel(
    const std::vector<PdfComponent>& components,
    const std::vector<double>& mean_yields,
    const ToyGeneratorConfig& cfg,
    unsigned long long toy_index,
    std::vector<Event>& out_events,
    std::vector<double>* out_generated_yields = nullptr
);
```

- 入力:
  - `components`: toy 生成に使う PDF 成分配列（`PdfComponent`）。
  - `mean_yields`: 各成分の平均期待事象数。`components` と同じ順序・同じ長さで与えます。
  - `cfg`: toy 生成設定（`ToyGeneratorConfig`）。seed や棄却法の `pmax` 推定条件を含みます。
  - `toy_index`: 同じ seed で複数 toy を区別するための通し番号。
  - `out_events`: 生成した toy イベント列を返す出力先。
  - `out_generated_yields`: 各成分で実際に Poisson 生成された事象数を受け取る任意出力。

- 出力:
  - 戻り値: toy 生成に成功したら `true`、入力不正やサンプリング失敗なら `false`。

### EvaluateUpperLimitPoint
- Header: `include/p2meg/UpperLimit.h`
- 目的: 固定 signal 仮説値 `s` 1点に対して、`q(s)=-2ln lambda(s)` の観測値と toy 分布を比較し、90% C.L. などの受容判定を返します。

- シグネチャ
```cpp
UpperLimitPointResult EvaluateUpperLimitPoint(
    const std::vector<Event>& events,
    const std::vector<PdfComponent>& components,
    const UpperLimitPointConfig& cfg
);
```

- 入力:
  - `events`: 実データのイベント列（`Event`）。
  - `components`: 尤度に使う PDF 成分配列（`PdfComponent`）。
  - `cfg`: 仮説値 `N_sig_test`、toy 本数、信頼水準、fit 設定をまとめた点ごとの設定。

- 出力:
  - 戻り値: その `s` 点における `q_obs(s)`、toy からの `p_value`、受容判定、実データの free/profile fit を含む結果構造体（`UpperLimitPointResult`）。

### EvaluateUpperLimitScan
- Header: `include/p2meg/UpperLimit.h`
- 目的: signal 仮説値 `s` を複数点走査し、受容された最大の `N_sig` を `N_sig^90` として返します。`N_mu_eff` が与えられていれば `BR_90` へも変換します。

- シグネチャ
```cpp
UpperLimitScanResult EvaluateUpperLimitScan(
    const std::vector<Event>& events,
    const std::vector<PdfComponent>& components,
    const UpperLimitScanConfig& cfg
);
```

- 入力:
  - `events`: 実データのイベント列（`Event`）。
  - `components`: 尤度評価と toy 生成に使う PDF 成分配列（`PdfComponent`）。
  - `cfg`: `N_sig` の走査点、toy 本数、`N_mu_eff`、fit 設定、乱数設定をまとめた scan 設定。

- 出力:
  - 戻り値: 通常 fit、各走査点の結果、`N_sig_90`、`BR_90` をまとめた結果構造体（`UpperLimitScanResult`）。

### AccTimeFitResult
- Header: `include/p2meg/AccTimeFit.h`
- 目的: ACC 時間形状 fit の結果（パラメータと適合度）をまとめます。

- シグネチャ
```cpp
struct AccTimeFitResult {
  int fit_status;
  double chi2;
  int ndf;
  double chi2_ndf;
  double A;
  double sigma;
  double C;
};
```

- 入力:
  - `fit_status`: ROOT fit の終了コード。
  - `chi2`: カイ二乗値。
  - `ndf`: 自由度。
  - `chi2_ndf`: `chi2/ndf`。
  - `A`: ガウス成分の振幅。
  - `sigma`: ガウス幅 [ns]。
  - `C`: pedestal 項の高さ。

- 出力:
  - 戻り値: （なし）

### AccTimeFit_DensityCoreValue
- Header: `include/p2meg/AccTimeFit.h`
- 目的: ACC 時間形状モデル `A * exp(-t^2 / (2 sigma^2)) + C` の密度値を直接返します。

- シグネチャ
```cpp
static inline double AccTimeFit_DensityCoreValue(double t, double A, double sigma, double C);
```

- 入力:
  - `t`: 時間差 [ns]。
  - `A`: ガウス成分の振幅。
  - `sigma`: ガウス幅 [ns]。
  - `C`: pedestal 項。

- 出力:
  - 戻り値: モデル密度値。`sigma<=0` や非有限値など不正入力では `0` を返します。

### AccTimeFit_Run
- Header: `include/p2meg/AccTimeFit.h`
- 目的: `TH1D` に対して、共通仕様の ACC 時間形状 fit を実行します。blind/reject 窓を除外した fit を macro と本解析コードで共通化するための入口です。

- シグネチャ
```cpp
static inline int AccTimeFit_Run(TH1D& h,
                                 double t_reject_min,
                                 double t_reject_max,
                                 double sigma_min,
                                 double sigma_max,
                                 AccTimeFitResult& out);
```

- 入力:
  - `h`: fit 対象の時間ヒストグラム。
  - `t_reject_min`: fit から除外する時間窓の下限 [ns]。
  - `t_reject_max`: fit から除外する時間窓の上限 [ns]。
  - `sigma_min`: `sigma` の下限制約 [ns]。
  - `sigma_max`: `sigma` の上限制約 [ns]。
  - `out`: fit 結果の書き込み先。

- 出力:
  - 戻り値: 成功なら `0`、空ヒストや不正 fit パラメータなどでは非 0 を返します。

### AccTimeFit_FillDensityHistogramFromResult
- Header: `include/p2meg/AccTimeFit.h`
- 目的: `AccTimeFitResult` のパラメータから、ビン平均密度のヒストグラムを作ります。

- シグネチャ
```cpp
static inline int AccTimeFit_FillDensityHistogramFromResult(TH1D& h_out,
                                                            const AccTimeFitResult& fit);
```

- 入力:
  - `h_out`: 出力先ヒストグラム。各ビンに密度値が書き込まれます。
  - `fit`: 使う fit 結果。

- 出力:
  - 戻り値: 成功なら `0`、`fit` が不正なら非 0 を返します。

### ToyGeneratorConfig
- Header: `include/p2meg/UpperLimit.h`
- 目的: toy dataset 生成の乱数・棄却法設定をまとめます。

- シグネチャ
```cpp
struct ToyGeneratorConfig {
    unsigned long long seed;
    int pmax_scan_trials;
    double pmax_safety;
    double pmax_update;
    int event_pool_size_per_component;
};
```

- 入力:
  - `seed`: 乱数 seed。
  - `pmax_scan_trials`: 棄却法で `pmax` を見積もる試行数。
  - `pmax_safety`: 見積もった `pmax` に掛ける安全係数。
  - `pmax_update`: 生成中に `pmax` を更新するときの倍率。
  - `event_pool_size_per_component`: 成分ごとに事前生成する event pool の大きさ。

- 出力:
  - 戻り値: （なし）

### UpperLimitPointConfig
- Header: `include/p2meg/UpperLimit.h`
- 目的: 固定 `N_sig` 1 点の upper limit 評価設定をまとめます。

- シグネチャ
```cpp
struct UpperLimitPointConfig {
    double N_sig_test;
    int n_toys;
    double cl;
    FitConfig free_fit_cfg;
    FitConfig prof_fit_cfg;
    ToyGeneratorConfig toy_cfg;
};
```

- 入力:
  - `N_sig_test`: 固定する signal yield 仮説値。
  - `n_toys`: 生成する toy 本数。
  - `cl`: 信頼水準（例: `0.90`）。
  - `free_fit_cfg`: 実データ free fit 用設定。
  - `prof_fit_cfg`: 固定 `N_sig` profile fit 用設定。
  - `toy_cfg`: toy 生成設定。

- 出力:
  - 戻り値: （なし）

### UpperLimitPointResult
- Header: `include/p2meg/UpperLimit.h`
- 目的: 固定 `N_sig` 1 点での受容判定結果をまとめます。

- シグネチャ
```cpp
struct UpperLimitPointResult {
    double N_sig_test;
    double q_obs;
    double p_value;
    double acceptance_threshold;
    bool accepted;
    int n_toys_requested;
    int n_toys_valid;
    FitResult fit_free_obs;
    FitResult fit_prof_obs;
};
```

- 入力:
  - `N_sig_test`: 評価した仮説値。
  - `q_obs`: 実データの検定統計量。
  - `p_value`: toy 分布に対する右側 p-value。
  - `acceptance_threshold`: 受容しきい値（通常 `1-CL`）。
  - `accepted`: 受容判定。
  - `n_toys_requested`: 要求 toy 数。
  - `n_toys_valid`: 実際に有効だった toy 数。
  - `fit_free_obs`: 実データ free fit 結果。
  - `fit_prof_obs`: 実データ profile fit 結果。

- 出力:
  - 戻り値: （なし）

### UpperLimitScanConfig
- Header: `include/p2meg/UpperLimit.h`
- 目的: 複数 `N_sig` 点の走査設定をまとめます。

- シグネチャ
```cpp
struct UpperLimitScanConfig {
    std::vector<double> N_sig_scan;
    int n_toys_per_point;
    double cl;
    double N_mu_eff;
    FitConfig free_fit_cfg;
    FitConfig prof_fit_cfg;
    ToyGeneratorConfig toy_cfg;
};
```

- 入力:
  - `N_sig_scan`: 走査する signal yield 値の配列。
  - `n_toys_per_point`: 各点の toy 本数。
  - `cl`: 信頼水準。
  - `N_mu_eff`: `BR = N_sig / N_mu_eff` 変換に使う有効停止ミューオン数。
  - `free_fit_cfg`: 実データ free fit 用設定。
  - `prof_fit_cfg`: 固定 `N_sig` profile fit 用設定。
  - `toy_cfg`: toy 生成設定。

- 出力:
  - 戻り値: （なし）

### UpperLimitScanResult
- Header: `include/p2meg/UpperLimit.h`
- 目的: `N_sig` 走査全体の結果をまとめます。

- シグネチャ
```cpp
struct UpperLimitScanResult {
    FitResult fit_free_obs;
    std::vector<UpperLimitPointResult> points;
    double N_sig_90;
    double BR_90;
    double N_mu_eff;
};
```

- 入力:
  - `fit_free_obs`: 実データ通常 fit の結果。
  - `points`: 各走査点の結果。
  - `N_sig_90`: 90% C.L. の signal upper limit。
  - `BR_90`: 90% C.L. の branching ratio upper limit。
  - `N_mu_eff`: 変換に使った有効停止ミューオン数。

- 出力:
  - 戻り値: （なし）

### NormalizationUncertaintyConfig
- Header: `include/p2meg/UpperLimit.h`
- 目的: `N_mu_eff` の公称値と不確かさをまとめ、BR 上限での normalisation uncertainty を指定します。

- シグネチャ
```cpp
struct NormalizationUncertaintyConfig {
    double N_mu_eff_nom;
    double N_mu_eff_sigma;
};
```

- 入力:
  - `N_mu_eff_nom`: 有効停止ミューオン数の公称値。
  - `N_mu_eff_sigma`: 有効停止ミューオン数の絶対誤差。

- 出力:
  - 戻り値: （なし）

### UpperLimitBRPointConfig
- Header: `include/p2meg/UpperLimit.h`
- 目的: 固定 BR 1 点の upper limit 評価設定をまとめます。

- シグネチャ
```cpp
struct UpperLimitBRPointConfig {
    double BR_test;
    int n_toys;
    double cl;
    FitConfig free_fit_cfg;
    FitConfig prof_fit_cfg;
    ToyGeneratorConfig toy_cfg;
    NormalizationUncertaintyConfig norm_cfg;
};
```

- 入力:
  - `BR_test`: 固定する分岐比仮説値。
  - `n_toys`: toy 本数。
  - `cl`: 信頼水準。
  - `free_fit_cfg`: 実データ free fit 用設定。
  - `prof_fit_cfg`: 固定 BR profile fit 用設定。
  - `toy_cfg`: toy 生成設定。
  - `norm_cfg`: `N_mu_eff` の公称値と不確かさ。

- 出力:
  - 戻り値: （なし）

### UpperLimitBRPointResult
- Header: `include/p2meg/UpperLimit.h`
- 目的: 固定 BR 1 点での受容判定結果をまとめます。

- シグネチャ
```cpp
struct UpperLimitBRPointResult {
    double BR_test;
    double N_sig_test_nominal;
    double q_obs;
    double p_value;
    double acceptance_threshold;
    bool accepted;
    int n_toys_requested;
    int n_toys_valid;
    FitResult fit_free_obs;
    FitResult fit_prof_obs;
};
```

- 入力:
  - `BR_test`: 評価した分岐比仮説値。
  - `N_sig_test_nominal`: 公称 `N_mu_eff` に対応する signal yield。
  - `q_obs`: 実データの検定統計量。
  - `p_value`: toy 分布に対する右側 p-value。
  - `acceptance_threshold`: 受容しきい値。
  - `accepted`: 受容判定。
  - `n_toys_requested`: 要求 toy 数。
  - `n_toys_valid`: 実際に有効だった toy 数。
  - `fit_free_obs`: 実データ free fit 結果。
  - `fit_prof_obs`: 実データ profile fit 結果。

- 出力:
  - 戻り値: （なし）

### UpperLimitBRScanConfig
- Header: `include/p2meg/UpperLimit.h`
- 目的: BR 走査による upper limit 設定をまとめます。

- シグネチャ
```cpp
struct UpperLimitBRScanConfig {
    std::vector<double> BR_scan;
    int n_toys_per_point;
    double cl;
    FitConfig free_fit_cfg;
    FitConfig prof_fit_cfg;
    ToyGeneratorConfig toy_cfg;
    NormalizationUncertaintyConfig norm_cfg;
};
```

- 入力:
  - `BR_scan`: 走査する分岐比の配列。
  - `n_toys_per_point`: 各点の toy 本数。
  - `cl`: 信頼水準。
  - `free_fit_cfg`: 実データ free fit 用設定。
  - `prof_fit_cfg`: 固定 BR profile fit 用設定。
  - `toy_cfg`: toy 生成設定。
  - `norm_cfg`: `N_mu_eff` の公称値と不確かさ。

- 出力:
  - 戻り値: （なし）

### UpperLimitBRScanResult
- Header: `include/p2meg/UpperLimit.h`
- 目的: BR 走査全体の結果をまとめます。

- シグネチャ
```cpp
struct UpperLimitBRScanResult {
    FitResult fit_free_obs;
    std::vector<UpperLimitBRPointResult> points;
    double BR_90;
    double N_sig_90_nominal;
    NormalizationUncertaintyConfig norm_cfg;
};
```

- 入力:
  - `fit_free_obs`: 実データ通常 fit の結果。
  - `points`: 各 BR 点の結果。
  - `BR_90`: 90% C.L. の branching ratio upper limit。
  - `N_sig_90_nominal`: 公称 `N_mu_eff` に対する signal upper limit。
  - `norm_cfg`: 使った normalisation uncertainty 設定。

- 出力:
  - 戻り値: （なし）

### ProfileLikelihoodQPoint
- Header: `include/p2meg/UpperLimit.h`
- 目的: 1 つの固定 `N_sig` 点で評価した profile-likelihood 統計量をまとめます。

- シグネチャ
```cpp
struct ProfileLikelihoodQPoint {
    double N_sig_test;
    double q_value;
    FitResult fit_prof;
};
```

- 入力:
  - `N_sig_test`: 固定した signal yield 仮説値。
  - `q_value`: `q = -2 ln lambda` の値。
  - `fit_prof`: その点の profile fit 結果。

- 出力:
  - 戻り値: （なし）

### ProfileLikelihoodQScanResult
- Header: `include/p2meg/UpperLimit.h`
- 目的: 同一 dataset に対して複数の固定 `N_sig` 点を評価した結果をまとめます。

- シグネチャ
```cpp
struct ProfileLikelihoodQScanResult {
    FitResult fit_free;
    std::vector<ProfileLikelihoodQPoint> points;
};
```

- 入力:
  - `fit_free`: 同一 dataset に対する free fit 結果。
  - `points`: 各固定仮説値に対する評価結果。

- 出力:
  - 戻り値: （なし）

### EvaluateProfileLikelihoodQ
- Header: `include/p2meg/UpperLimit.h`
- 目的: 固定 `N_sig` 仮説値に対して `q(s) = -2 ln lambda(s)` を計算します。free fit と profile fit の両方を返します。

- シグネチャ
```cpp
double EvaluateProfileLikelihoodQ(
    const std::vector<Event>& events,
    const std::vector<PdfComponent>& components,
    const FitConfig& free_fit_cfg,
    const FitConfig& prof_fit_cfg,
    double N_sig_fixed,
    FitResult& fit_free_out,
    FitResult& fit_prof_out
);
```

- 入力:
  - `events`: 評価対象 dataset。
  - `components`: PDF 成分配列。
  - `free_fit_cfg`: 通常 fit の設定。
  - `prof_fit_cfg`: 固定 `N_sig` profile fit の設定。
  - `N_sig_fixed`: 固定する signal yield 仮説値。
  - `fit_free_out`: free fit 結果の出力先。
  - `fit_prof_out`: profile fit 結果の出力先。

- 出力:
  - 戻り値: `q(s)` の値。数値誤差で負になった場合は `0` に丸められます。

### EvaluateProfileLikelihoodQScan
- Header: `include/p2meg/UpperLimit.h`
- 目的: 同一 dataset に対して複数の固定 `N_sig` 仮説値をまとめて評価します。free fit は 1 回だけ行い、各点で再利用します。

- シグネチャ
```cpp
ProfileLikelihoodQScanResult EvaluateProfileLikelihoodQScan(
    const std::vector<Event>& events,
    const std::vector<PdfComponent>& components,
    const FitConfig& free_fit_cfg,
    const FitConfig& prof_fit_cfg,
    const std::vector<double>& N_sig_scan
);
```

- 入力:
  - `events`: 評価対象 dataset。
  - `components`: PDF 成分配列。
  - `free_fit_cfg`: 通常 fit の設定。
  - `prof_fit_cfg`: 固定 `N_sig` profile fit の設定。
  - `N_sig_scan`: まとめて評価する signal yield 仮説値の配列。

- 出力:
  - 戻り値: free fit と各走査点の `q` をまとめた結果（`ProfileLikelihoodQScanResult`）。

### EvaluateUpperLimitBRPoint
- Header: `include/p2meg/UpperLimit.h`
- 目的: 固定 BR 1 点に対して、normalisation uncertainty を含む toy MC 受容判定を返します。

- シグネチャ
```cpp
UpperLimitBRPointResult EvaluateUpperLimitBRPoint(
    const std::vector<Event>& events,
    const std::vector<PdfComponent>& components,
    const UpperLimitBRPointConfig& cfg
);
```

- 入力:
  - `events`: 実データのイベント列。
  - `components`: 尤度評価と toy 生成に使う PDF 成分配列。
  - `cfg`: BR 仮説値、toy 本数、信頼水準、fit 設定、normalisation uncertainty をまとめた設定。

- 出力:
  - 戻り値: 固定 BR 1 点での `q_obs`、toy 由来 `p_value`、受容判定、fit 結果をまとめた構造体（`UpperLimitBRPointResult`）。

### EvaluateUpperLimitBRScan
- Header: `include/p2meg/UpperLimit.h`
- 目的: BR 仮説値を複数点走査し、受容された最大の BR を `BR_90` として返します。

- シグネチャ
```cpp
UpperLimitBRScanResult EvaluateUpperLimitBRScan(
    const std::vector<Event>& events,
    const std::vector<PdfComponent>& components,
    const UpperLimitBRScanConfig& cfg
);
```

- 入力:
  - `events`: 実データのイベント列。
  - `components`: 尤度評価と toy 生成に使う PDF 成分配列。
  - `cfg`: BR 走査点、toy 本数、信頼水準、fit 設定、normalisation uncertainty をまとめた設定。

- 出力:
  - 戻り値: 通常 fit、各 BR 点の結果、`BR_90`、公称 `N_sig_90` をまとめた構造体（`UpperLimitBRScanResult`）。
