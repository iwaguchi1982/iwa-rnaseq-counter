# iwa-rnaseq-counter

FASTQ から RNA-Seq の count データを作成する、**Wet ファースト**の解析アプリです。  
まずは Web アプリとして、**入力をそろえて、実行して、結果を保存する** ところまでを分かりやすく使えることを重視しています。

このアプリは、次のような方を想定しています。

- Wet 実験が主な研究者
- RNA-Seq 解析を自分でも少し進めたい人
- Dry 解析の専門家ではないが、count データまでは自分で作りたい人
- 解析担当者へ渡す前の結果をそろえたい人

---

## 1. このアプリでできること

- FASTQ から RNA-Seq の定量を行う
- gene / transcript レベルの count を出力する
- サンプルごとの結果を保存する
- 複数サンプルの結果をまとめて、次の解析で使いやすい形にする
- 実行時の設定や結果を保存して、あとで見直しやすくする

---

## 2. まずは Web アプリとして使う

このアプリは、まず Web アプリとして使うことを想定しています。

```bash
pixi run streamlit run app.py
```
> 以前の iwa_rnaseq_counter.py は app.py に名称変更されています。

### Web アプリでの基本的な流れ
1. サンプル認識モードを選択する（FASTQ直接指定 / CSV読込 / ディレクトリ自動認識）
2. サンプル情報や比較群、除外フラグ（Exclude）を整理する
3. `RUN START` で実行する
4. 実行中・完了・失敗の状態を Job History で確認する
5. count 結果やレポート連携用の成果物を保存する
6. 必要に応じて、次の解析アプリ（reporter）へ進む

> 単体サンプルでも複数サンプルでも、同じ `sample_sheet.csv` 形式で扱えるようにしています。
---
### 実行後の見方
実行した解析は、画面左の Job History から後で見直せます。

- 現在実行中の job
- 完了した job
- 失敗した job

を一覧で確認できます。

一覧は `Runs Root` で指定した保存先配下から読み込みます。  
通常は既定の `output/` をそのまま使えば大丈夫です。

保存先 root は切り替えられますが、各解析結果の内部構造は固定です。  
そのため、保存場所を分けても、同じ見方で結果を追えます。
---
### Web アプリで意識すればよいこと
まず利用時に気にすればよいのは次の点です。

- どの FASTQ を使うか
- single-end か paired-end か
- サンプル名をどう付けるか
- 何が出力されるか

内部のデータ形式や開発用の仕様は、使うだけなら意識しなくて大丈夫です。
---
## 3. CLI で使う場合
CLI でも実行できます。
ただし、このアプリでは 単体サンプルでも複数サンプルでも、入口は sample sheet.csv に統一する方針です。
基本コマンド
```bash
pixi run python cli.py run-batch \
  --sample-sheet path/to/sample_sheet.csv \
  --salmon-index path/to/salmon_index \
  --tx2gene path/to/tx2gene.tsv \
  --outdir output/run_001
```
この形式なら、

- 1サンプルだけのときも
- 複数サンプルまとめてのときも

同じ考え方で使えます。

## 4. sample sheet とは何か
sample sheet は、**どのサンプルを、どの FASTQ で解析するか** をまとめた CSV です。
このアプリでは、1行 = 1 assay の考え方で扱います。

例: example.sample.sheet.csv
```csv
sample_id,r1_path,r2_path,layout,group,condition,replicate,batch,pair_id,display_name,color,exclude,note
SP001,/data/fastq/SP001_R1.fastq.gz,/data/fastq/SP002_R2.fastq.gz,paired,case,baseline,1,batch1,,Case 1,#1f77b4,false,
SP002,/data/fastq/SP002_R1.fastq.gz,/data/fastq/SP002_R2.fastq.gz,paired,case,baseline,2,batch1,,Case 2,#1f77b4,false,
SP003,/data/fastq/SP003_R1.fastq.gz,/data/fastq/SP003_R2.fastq.gz,paired,control,baseline,1,batch1,,Control 1,#ff7f0e,false,
SP004,/data/fastq/SP004_R1.fastq.gz,,single,control,baseline,2,batch2,,Control 2,#ff7f0e,false,single-end example
```

## 5. sample_sheet.csv の列構成（v0.5.0 契約）
解析の入力となる CSV は以下の列で構成されます。  
人間向けの入口（GUI / CSV）はここで行われ、内部的には `AssaySpec` 等の内部形式に変換・正規化されます。

### 必須列 (Mandatory)
- `sample_id`: サンプル識別子。結果や表示の主キーとなります。
- `r1_path`: R1 FASTQ への絶対パスまたは相対パス。
- `layout`: `single` または `paired`。

### 条件付き必須 (Conditional)
- `r2_path`: `layout=paired` の場合に必須です。`single` の場合は空欄にしてください。

### 任意列 (Optional/Metadata)
比較グループや図示の制御に使用されます。指定がない場合は空文字などで補完されます。
- `group`, `condition`, `replicate`, `batch`, `pair_id`: 解析の群分け・要因指定用
- `display_name`: ユーザーインターフェースでの表示名（未指定時は `sample_id`）
- `color`: 表示色（16進数カラーコード等）
- `exclude`: `true` で解析から除外
- `note`: 自由記入欄

## 6. 出力物と Run Artifact（v0.5.0 契約）
解析が実行されると、出力先ディレクトリに以下の標準構造でファイルが生成されます。

### フォルダ構成
- `inputs/`: 実行に使用された入力設定（`run_config.json`, `sample_sheet.csv`, `gui_input_status.csv`）
- `work/`: 解析の中間生成物
- `results/`: 解析の最終成果物（reporter 連携用。下記参照）
- `specs/`: 解析仕様書（`matrix.spec.json`, `execution-run.spec.json`）
- `logs/`: 実行記録（`run.log`）
- `tmp/`: 一時ファイル

### 主要な出力ファイル
- `results/gene_numreads.csv`: 遺伝子ごとのカウント（メイン成果物）
- `results/feature_annotation.tsv`: **※重要** Reporter 用のアノテーションファイル。`tx2gene` にシンボル情報が含まれる場合のみ生成されます。
- `dataset_manifest.json`: Reporter が解析結果をロードするためのカタログファイル。
- `inputs/gui_input_status.csv`: GUI 実行時のパラメータや実行日時を記録したトレーサビリティ用ファイル。

> [!NOTE]
> **v0.5.0 契約の核**:
> - `feature_annotation.tsv` が存在しなくても run 自体は成立します。
> - Reporter はアノテーションがない場合、自動的に `feature_id` にフォールバックして動作します。

## 7. 「Suite 契約の完了」の定義
v0.5.0 において、以下の状態が満たされていることを「追従完了」として定義しています。
1. **実装と文書の一致**: `sample_parser.py` が受け取れる 13 列と `README` の説明、および `example.sample.sheet.csv` の構成が完全に一致していること。
2. **出力位置の固定**: `feature_annotation.tsv` が必ず `results/` 配下に配置され、その有無が `MatrixSpec` のメタデータと一貫していること。
3. **fallback の保証**: アノテーションが生成されないケースでも、`dataset_manifest.json` を通じて Reporter が安全に起動・表示できること。

Web アプリでは、これらの結果を Job History から後で見直せます。  
Job History は `Runs Root` で指定した保存先配下の解析結果を一覧表示します。  
通常は既定の `output/` をそのまま使う想定です。

## 8. 動作確認済み条件 (v0.5.0)
- OS: Linux (WSL2 / Ubuntu 24.04 等)
- Python: 3.12+
- Salmon: v1.10.1 (推奨)
  - v1.11.x 以降ではインデックスの versionInfo.json が厳格チェックされるため注意が必要です。
- データ: Saccharomyces cerevisiae (Yeast) R64-1-1 / Human (GRCh38)

## 9. 注意事項
- 現バージョンでは Salmon を使った解析のみをサポートしています
- 参照データとして、Salmon index と tx2gene を事前に準備してください
- **Suite 連携**: reporter で遺伝子名（gene_symbol）を表示したい場合は、`results/feature_annotation.tsv` が存在することを確認してください。これは `tx2gene` にシンボル情報が含まれていれば自動生成されます。

## 11. Spec Contract (v0.5.0 技術仕様)
Suite 間のデータ受け渡しに使用される `MatrixSpec` の正式な振る舞いです。

### feature_annotation_path の扱い
- **定義**: Reporter での表示に使用するアノテーションファイルへの参照です。
- **値の型**: `Optional[str]` (JSON 上では `string` または `null`)。
- **格納条件**:
  - アノテーション生成に成功した場合のみ、`results/feature_annotation.tsv` へのパス（絶対または相対）を格納します。
  - 生成されなかった場合は `null` (None) とします。空文字 `""` は使用しません。
- **禁止事項**: 内部集約用の `tx2gene` へのパスを代用品として流し込むことは禁止されています。

### feature_annotation_available の同期
- `MatrixSpec.metadata["feature_annotation_available"]` は、`feature_annotation_path` が有効なファイルを指している場合のみ `true` となります。
- Reporter はこのフラグを第一のヒントとしてアノテーションのロードを試みます。

## License

This repository is distributed under the Iwa Collections Non-Resale License 1.0.
Commercial resale of the software itself, or paid redistribution of derivative versions where the software is the primary value, is prohibited.

本リポジトリは Iwa Collections Non-Resale License 1.0 で公開しています。
ソフトウェア自体の有償販売、および本ソフトウェアが主たる価値となる派生物の有償再配布は禁止です。

