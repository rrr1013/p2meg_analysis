# Repository Structure Guide

このメモは、今後このリポジトリを拡張するときに「どこへ何を置くか」を決めるための基準です。

## 基本原則

1. 物理ロジックは `include/p2meg/` と `src/` に寄せる
2. 生波形から解析入力を作る executable は `datashaping/`
3. 端末から再実行する解析入口は `scripts/`
4. 図作成と ROOT macro は `macros/`
5. 一回限りの切り分けと quicklook は `check/`
6. 生データや生成データは `data/` に閉じ込める
7. 手書き文書と生成物は `doc/` で区別する

## 現在の repo で困りやすい点

- thesis 固有コードと再利用コードが同じ階層にある
- `doc/` に生成物が混ざりやすい
- `scripts/build_all.sh` と `datashaping/` のビルド入口が分かれていた
- README が不足していて、入口がコード名からしか分からない

## これからの置き方

### `include/p2meg/` と `src/`

再利用する関数・クラスのみを置く。

### `datashaping/`

以下のような「再構成段階」の executable のみを置く。

- 波形分解
- baseline 推定
- PS / NaI パルス抽出
- run ごとの最終 5D 入力作成

### `scripts/`

後輩が端末から実行する正式入口を置く。

- grid PDF 生成
- NLL fit
- upper limit / sensitivity
- mockdata 生成
- 要約用 CLI

### `macros/`

ROOT で図を描くもの、見た目を確認するものを置く。

### `check/`

一回限りの切り分け、run 固有の調査、quicklook を置く。
ここにあるものは「正式フローではない」と考える。

### `doc/`

文書は次の考えで残す。

- 残す: `.md`, `.tex`, 要約 `.txt`
- 原則残さない: `.aux`, `.dvi`, `.log`
- 大量 scan の数値結果は、必要なら別の results 保管先も検討する

## 後輩向けのおすすめ開始順

1. `README.md`
2. `data/README.md`
3. `datashaping/README.md`
4. `scripts/README.md`
5. 必要に応じて `FUNCTIONS.md`
