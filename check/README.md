# check

このディレクトリは quicklook と探索的調査用です。

## 役割

- 生波形の目視確認
- ADC 積分やピークの健全性確認
- run 固有の切り分け
- selection の再評価

## 現在の中身

- `eventdisplay.cc`: quicklook の主入口
- `hist_*`, `study_*`, `report_*`, `reassess_*`: run 固有・問題切り分け寄り

## 後輩向けルール

- 正式パイプラインに組み込むものは `datashaping/` または `scripts/` に昇格する
- その場限りの検証や run 固有の解析はここに置く
- thesis 専用コードは README やコメントで「再利用前提ではない」ことを明示する
