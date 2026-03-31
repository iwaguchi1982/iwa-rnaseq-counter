# Developer_README

## この文書の目的

この文書は `iwa-rnaseq-counter` の開発者向け内部整理である。対象は、現状の状態、基準、契約、データフロー、および変更時に崩してはならない前提の共有にある。利用者向けの導線や説明は `README` に置き、この文書では内部契約と実装上の判断基準を扱う。

## このアプリの責務

`iwa-rnaseq-counter` の責務は、RNA-Seq の入力を受け取り、実行条件を正規化し、後段が安全に読める共通成果物へ変換することにある。主役は GUI と `sample_sheet.csv` だが、内部ではそれらを Spec-aware な中間成果物へ正規化して扱う。

このアプリは次を担当する。

- FASTQ と sample sheet を入口として受ける
- 実行条件を検証する
- quant 実行を backend job として起動する
- run artifact を標準構造で出力する
- 後段の reporter が読める matrix / metadata / annotation 契約へ寄せる

このアプリは reporter の責務まで侵食しない。解釈、比較、図示、レポート生成は reporter 側の責務である。

## 現状の状態

現状では、GUI / CLI の両方から run を作成できる。Runs Root、Job History、active job 復元、疑似タブによる Run/Input と Job History の分離が入っており、利用者向けの導線は一定程度整理済みである。

内部では、次の流れが成立している。

- `sample_sheet.csv` もしくは GUI 入力を受ける
- 実行設定を整理する
- backend 実行を起動する
- matrix / execution / supporting input を run artifact に出力する
- annotation を作れる場合は `results/feature_annotation.tsv` を生成する

v0.5 系の主題は、ここに reporter 接続を見据えた annotation / contract 固定を加えることにある。

## 入口契約

### `sample_sheet.csv`

`sample_sheet.csv` は利用者向けの正式入口である。内部では Spec に正規化するが、利用者に JSON を直接書かせない。

正式列契約は次の通りである。

必須列
- `sample_id`
- `r1_path`
- `layout`

条件付き必須
- `r2_path` ただし paired のとき

任意列
- `group`
- `condition`
- `replicate`
- `batch`
- `pair_id`
- `display_name`
- `color`
- `exclude`
- `note`

`layout=single` のとき `r2_path` は空でよい。`layout=paired` のとき `r1_path` と `r2_path` の両方が必要である。

### GUI 入力

GUI は利用者向けの主入口である。GUI 上の入力は、そのまま内部真実ではない。内部真実は sample table と Spec 群であり、GUI 入力はそこへ正規化される。

## 内部契約

### `ExecutionRunSpec`

`ExecutionRunSpec` は、何をどう実行したかの実行記録である。実行 backend、command、リソース、開始時刻、状態、出力先を保持する。ツール差はここへ寄せる。

### `MatrixSpec`

`MatrixSpec` は、後段が読む matrix の契約である。Reporter は backend 固有出力を直接読むのではなく、最終的にこの契約へ正規化された成果物を読む。

`MatrixSpec` は少なくとも次の意味を持つ。

- 何の matrix か
- どの feature space を持つか
- どの run から来たか
- どの annotation を参照するか

### `feature_annotation_path` のルール

`feature_annotation_path` は、reporter 表示用 annotation 契約ファイルへの参照として扱う。annotation 生成に成功した場合のみ `results/feature_annotation.tsv` へのパスを設定し、生成できなかった場合は `None` (JSON では `null`) とする。`tx2gene` は gene 集約用資源であり、`feature_annotation_path` の代用品として設定しない。`tx2gene` と `feature_annotation.tsv` は同一由来でも意味を区別し、前者は集約用、後者は reporter における `gene_symbol` / `display_label` 供給用 annotation 契約として扱う。空文字 `""` は使用しない。

### `feature_annotation_available` のルール

`feature_annotation_available` は実態と一致させる。annotation を生成できたときのみ `True` とし、失敗時は `False` とする。存在しない annotation を存在するものとして metadata に書かない。

## annotation 契約

### `tx2gene` の役割

`tx2gene` は transcript から gene への集約用資源である。集約時の参照には使うが、そのまま reporter 表示用 annotation 契約とみなしてはならない。

### `feature_annotation.tsv` の役割

`feature_annotation.tsv` は reporter 表示用 annotation 契約である。最低限、次の列を持つ。

- `feature_id`
- `gene_symbol`

`gene_symbol` が取れない場合は annotation を作らず、reporter 側で `feature_id` fallback を許容する。偽の annotation を作って意味を曖昧にしない。

### 標準配置

annotation を生成できた場合の標準配置は `results/feature_annotation.tsv` とする。GUI backend と CLI / runner は同じ helper でこの配置を決める。

## データフロー

### GUI / CLI から run artifact まで

1. 利用者が GUI または `sample_sheet.csv` で入力する
2. validator が入力条件を確認する
3. backend 実行用 command を組み立てる
4. run directory を標準構造で作る
5. matrix / execution / supporting input を出力する
6. `tx2gene` から annotation を作れる場合は `results/feature_annotation.tsv` を置く
7. `MatrixSpec.feature_annotation_path` に annotation 参照を設定する

### counter から reporter まで

counter は後段に渡す共通成果物を作る。reporter は backend 固有出力ではなく、matrix / execution / annotation 契約を読む。

## 実装上の基準

### 保つべきこと

- GUI と CLI で run artifact の意味がずれない
- annotation path 決定を GUI backend と runner で分岐させない
- `feature_annotation_path` に tx2gene を混ぜない
- 利用者向け README に内部契約を出しすぎない
- reporter の責務を counter に持ち込まない

### 避けるべきこと

- `feature_annotation.tsv` を作れないときに偽の path を設定すること
- GUI 側だけ annotation を特別扱いし、CLI とずれること
- matrix 契約より先に backend 固有出力へ依存すること
- 利用者向け README に内部ルールを過剰に書くこと

## 変更時の確認項目

- sample sheet 契約を崩していないか
- GUI と CLI の出力説明がずれていないか
- `results/feature_annotation.tsv` の生成条件が一貫しているか
- `feature_annotation_path` の意味が曖昧になっていないか
- reporter へ渡す成果物説明が README と一致しているか

## project_root 配下の必要ディレクトリと主要ファイル

以下は、開発時に意味を理解しておくべき最小構成である。

### `project_root/app.py`

Streamlit の GUI 入口である。Runs Root、Job History、疑似タブ、Run/Input 画面、結果表示の大枠を持つ。

### `project_root/cli.py`

CLI 入口である。GUI 非依存の実行導線を持つ。

### `project_root/ui/sections.py`

GUI の表示部品をまとめる。利用者向け表示を整理する層であり、実行契約の真実をここに持ち込まない。

### `project_root/src/iwa_rnaseq_counter/models/`

内部契約モデル群を置く。`assay.py`、`matrix.py`、`execution_run.py` が主であり、counter の Spec-aware な中心である。

### `project_root/src/iwa_rnaseq_counter/io/`

Spec / sample sheet の入出力を置く。読み書きの境界をここへ寄せる。

### `project_root/src/iwa_rnaseq_counter/pipeline/`

実行フロー本体を置く。`runner.py` が CLI 系、`gui_backend.py` が GUI backend 系の主要入口である。

### `project_root/src/iwa_rnaseq_counter/builders/`

GUI artifact export のような組み立て処理を置く。

### `project_root/src/iwa_rnaseq_counter/legacy/`

現行実装でまだ重要なロジックを含む。`annotation_helper.py`、`sample_parser.py`、`validators.py`、`salmon_runner.py`、`run_artifacts.py` などがある。段階的移行前提のため、ここを雑に壊さない。

### `project_root/README.md`

利用者向け README である。GUI / sample sheet / 出力物の意味を説明する場であり、内部契約の細部はここに出しすぎない。

### `project_root/Developer_README.md`

この文書である。現状の状態、基準、契約、データフローを開発者向けに整理する。
