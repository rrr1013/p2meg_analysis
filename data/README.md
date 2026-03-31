# data

このディレクトリはローカル専用です。原則として git 管理しません。

## 置くもの

- `rawdata/`: digitizer の生波形
- `shapeddata/`: run ごとの再構成中間物と `final_module_pulses_runM.txt`
- `finaldata/`: run 結合済みの最終 5D 入力
- `pdf_cache/`: RMD / ACC 格子 PDF
- `mockdata*`: toy / mock data
- `toy_cache/`: toy 計算の中間キャッシュ

## 基本方針

- 生データは再取得困難なので、別媒体にも保存する
- `finaldata/` に置くファイルは「後段の解析コードが直接読むもの」だけにする
- 一時ファイルやデバッグ出力は run ごとの作業ディレクトリに閉じ込める

## 後輩向けメモ

- まず `rawdata/` の命名をそろえる
- 次に `run_datashaping_pipeline` で run ごとの `shapeddata/` を作る
- 最後に必要な run を結合して `finaldata/` を作る
