# macros

このディレクトリは ROOT macro 用です。

## 何を置くか

- 解析結果の可視化
- final figure 生成
- exploratory な分布確認

## 何を置かないか

- 生波形から解析入力を作る処理
- 長時間の大量 toy を回す主入口
- 汎用ライブラリとして再利用したい実装

## 現在の使い分け

- `plot_*_thesis.C` や `*_final.C` は thesis / final figure 寄り
- `fit_*`, `predict_*` は解析補助
- 一回限りの run 固有調査は本来 `check/` に寄せるのが望ましい

## 後輩向けメモ

- macro は「図を出す」「見た目を確認する」用途を優先する
- 本解析の正式な数値結果を出す入口は、できるだけ `scripts/` 側に置く
