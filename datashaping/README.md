# datashaping

このディレクトリは、生波形から解析入力を作るための実行ファイル群です。

## 現行の正式フロー

1. `wave_to_pileup_csv.cpp`
2. `wave_to_baseline_template.cpp`
3. `run_pileup_residual_standalone.cpp`
4. `run_pssignal_standalone.cpp`
5. `build_module_summary.cpp`
6. `final_pair_observables.cpp`

通常はこれらを個別に叩かず、`run_datashaping_pipeline.cpp` を入口にします。

## 代表的な出力

- `final_module_pulses_runM.txt`
- `unmatched_ps_runM.csv`
- `step3f_runM_pp.txt`
- `step3f_runM_pg.txt`
- `step3f_runM_gg.txt`
- `step4f_runM_eg.txt`

## 単位の注意

- `final_module_pulses_runM.txt` の時刻列は ns
- `step4f_runM_eg.txt` の `t` は現状では秒
- 解析ライブラリの `Event.t` は ns を想定している

## 置き方のルール

- 生波形から final 5D 入力を作る executable はここに置く
- ライブラリとして再利用するものは `src/` と `include/p2meg/` に置く
- final figure 作成や一回限りの調査はここに置かない
