# scripts

このディレクトリは、CLI 実行ファイルとして使う解析スクリプト群です。

## 主な分類

### 本解析の主要入口

- `run_nll_fit.cc`
- `run_upper_limit_toymc.cc`
- `run_upper_limit_br_toymc.cc`
- `run_sensitivity_br_toymc.cc`
- `make_acc_grid_pdf.cc`
- `make_rmd_grid_pdf_main.cc`

### 診断・検証

- `diagnose_fc_qtoy_steps.cc`
- `diagnose_nsig_bound_effect.cc`
- `benchmark_fc_parts.cc`
- `evaluate_rmd_aw_sm.cc`
- `inspect_rmd_grid.cc`

### mock / toy data 生成

- `make_acc_mockdata.cc`
- `make_pdfmix_mockdata.cc`
- `make_michel_E_mockdata.cc`
- `merge_mockdata.cc`
- `smear_finaldata.cc`

### 図化・要約

- `plot_acc_grid_pdf.cc`
- `plot_rmd_theory.cc`
- `plot_signal_rmd_pdfs.cc`
- `plot_pdf_shapes_final.cc`
- `plot_fc_toy_components.cc`
- `plot_aw_result_megstyle.cc`

## 運用ルール

- 後輩が「端末から再実行する入口」は原則ここに置く
- ROOT macro で十分なものは `macros/` に置く
- run 固有の一回限り調査は `check/` に寄せる
- ファイル名は `run_*`, `make_*`, `plot_*`, `diagnose_*` のように役割を明示する
