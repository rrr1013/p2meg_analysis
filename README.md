# p2MEG Analysis Repository

このリポジトリは、課題研究 P2 の p2MEG 実験に使った解析・再構成コードをまとめたものです。

対象読者:
- 実験を引き継ぐ後輩
- 再解析や条件変更を行う人
- final likelihood / Michel 解析の入力を作り直したい人

## まず読むもの

- 全体方針と運用ルール: [AGENTS.md](./AGENTS.md)
- 公開 API の一覧: [FUNCTIONS.md](./FUNCTIONS.md)
- ディレクトリ整理方針: [doc/REPO_STRUCTURE.md](./doc/REPO_STRUCTURE.md)
- 文書・結果の置き方: [doc/README.md](./doc/README.md)

## 最短の入口

1. ビルド

```bash
./scripts/build_all.sh
```

2. 生波形から run ごとの再構成出力を作る

```bash
./build/run_datashaping_pipeline data/rawdata/run8
```

3. `final_module_pulses_runM.txt` から最終 5D 入力を作る

```bash
./build/final_pair_observables --input-dir data/shapeddata/run8 --run 8
```

4. 生成した 5D データを使って PDF 生成・fit・upper limit を回す

```bash
./build/run_nll_fit data/finaldata/step4f_run1to20_eg_doubleonly_allpatterns_500ns.txt
./build/run_upper_limit_br_toymc ...
```

## 現行の主な流れ

### 1. 再構成

生波形 `wave_<CH>_runM.txt` を入力し、[datashaping](./datashaping/README.md) のパイプラインで次を作ります。

- `final_module_pulses_runM.txt`
- `unmatched_ps_runM.csv`
- `step3f_runM_{pp,pg,gg}.txt`
- `step4f_runM_eg.txt`

### 2. 最終解析

最終 5D 入力 `(Ee, Eg, t, phi_detector_e, phi_detector_g)` を使い、[scripts](./scripts/README.md) の CLI と [macros](./macros/README.md) の ROOT macro で以下を行います。

- RMD / ACC / signal PDF の生成
- NLL fit
- Feldman-Cousins / toy MC upper limit
- thesis 用の図作成

### 3. 健全性確認と一回限りの調査

[check](./check/README.md) には quicklook と run 固有の調査 macro を置いてあります。
ここは「正式パイプライン」ではなく、波形確認や切り分け用です。

## ディレクトリの見方

- `include/p2meg/`: 公開ヘッダ
- `src/`: ライブラリ実装
- `datashaping/`: 生波形から解析入力を作る実行ファイル群
- `scripts/`: CLI 実行ファイル群
- `macros/`: ROOT macro
- `check/`: quicklook / 探索的調査
- `data/`: ローカルデータ置き場（git 管理しない）
- `doc/`: 手書き文書・整理メモ・結果

## 後輩向けの運用ルール

- まず `README.md` と `doc/REPO_STRUCTURE.md` を読む
- 新しいコードを追加するときは、どの段階のものかを先に決めてから置き場所を選ぶ
- run 固有・卒論固有の一回限りコードは、正式入口に混ぜず `check` や `doc` 側で明示する
- 生成物（aux, log, dvi, 実行ファイル, 生データ）は原則コミットしない

## 注意

- `Event.t` は解析コード内部では ns を想定しています
- `final_pair_observables` が出す `step4f_runM_eg.txt` の `t` は秒なので、後段へ渡す前に単位を確認してください
- このリポジトリには自動テスト基盤はまだありません。最低限、ビルドと代表コマンドの再実行確認をしてください
