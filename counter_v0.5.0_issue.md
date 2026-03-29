# v0.5.0 issue (Counter)

## この repo での目的
`iwa-rnaseq-counter` 側では、v0.5.0 において sample sheet 契約を固定し、GUI 実行時の証跡出力を安定化する。
また、README と example を現行の 13 列契約に完全同期させる。

## スコープ
### 1. sample_sheet.csv 契約の最終固定
- 13 列（sample_id, r1_path, r2_path, layout, group, condition, replicate, batch, pair_id, display_name, color, exclude, note）の標準化。
- example と README の同期。

### 2. GUI 実行時の証跡出力（gui_input_status.csv）
- GUI で設定した解析名、入力パス、Index パス、スレッド数、strandedness 等を `inputs/gui_input_status.csv` として保存する。

### 3. 出力 Spec の提示
- `matrix.spec.json` および `execution-run.spec.json` が出力されることを README に明記する。

## タスクリスト
- [x] sample_sheet.csv 列定義の固定
- [x] example.sample.sheet.csv の更新
- [ ] app.py で gui_input_status.csv を出力する
- [ ] README.md の「出力されるもの」を最新化する
