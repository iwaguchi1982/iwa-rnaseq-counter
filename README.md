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
4. count 結果やレポート連携用の manifest 等を保存する
5. 必要に応じて、次の解析アプリ（reporter）へ進む
> 単体サンプルでも複数サンプルでも、同じ `sample_sheet.csv` 形式で扱えるようにしています。

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
SP001,/data/fastq/SP001_R1.fastq.gz,/data/fastq/SP001_R2.fastq.gz,paired,case,baseline,1,batch1,,Case 1,#1f77b4,false,
SP002,/data/fastq/SP002_R1.fastq.gz,/data/fastq/SP002_R2.fastq.gz,paired,case,baseline,2,batch1,,Case 2,#1f77b4,false,
SP003,/data/fastq/SP003_R1.fastq.gz,/data/fastq/SP003_R2.fastq.gz,paired,control,baseline,1,batch1,,Control 1,#ff7f0e,false,
SP004,/data/fastq/SP004_R1.fastq.gz,,single,control,baseline,2,batch2,,Control 2,#ff7f0e,false,single-end example
```

## 5. sample sheet の各列
### 必須列
- sample_id
  サンプル名です。結果表の列名や表示名の元になります。
- r1_path
  R1 FASTQ のパスです。single-end でも使います。
- layout
  single または paired を指定します。
- paired-end のときに必要
  r2_path
- layout=paired のとき必須です。layout=single のときは空欄で構いません。
### 任意列
- group
  比較群の目安です。例: case, control
- condition
  条件名のメモです。例: baseline, treated
- replicate
  反復番号や区別用の情報です。
- batch
  バッチ情報です。
- pair_id
  対応のあるサンプルを扱いたいときの識別子です。
- display_name
  画面表示用の名前です。
- color
  表示色の指定です。
- exclude
  true の場合、その行を解析対象から外せるようにするための列です。
- note
  補足メモです。

## 6. single-end / paired-end の書き方
### paired-end の場合
- layout=paired
- r1_path, r2_path の両方を記入
### single-end の場合
- layout=single
- r1_path を記入
- r2_path は空欄

## 7. 出力されるもの

主に次のものが出力されます。

- gene / transcript の count 結果
- 実行結果の保存フォルダ
- `dataset_manifest.json`
- `inputs/auto_generated.sample_sheet.csv`（FASTQ 直読み込み時）
- `inputs/gui_input_status.csv`（GUI 実行時）
- 次の解析アプリに渡しやすい形にまとめた結果
- 実行記録や内部設定ファイル（内部連携用）


## 8. 動作確認済み条件 (v0.2.0)
- OS: Linux (WSL2 / Ubuntu 24.04 等)
- Python: 3.12+
- Salmon: v1.10.1 (推奨)
  - v1.11.x 以降ではインデックスの versionInfo.json が厳格チェックされるため注意が必要です。
- データ: Saccharomyces cerevisiae (Yeast) R64-1-1

## 9. 注意事項
- 現バージョンでは Salmon を使った解析のみをサポートしています
- 参照データとして、Salmon index と tx2gene を事前に準備してください
- インデックス作成時のバージョン違いによって動作差が出ることがあります

## 10. 開発者向けメモ
ここから下は、将来の自分向けの短いメモです。
利用だけが目的なら読まなくて大丈夫です。

### 10.1 現在の役割
iwa-rnaseq-counter は RNA-Seq Suite における count 作成の入口です。

### 10.2 現在の入口
人間向けの入口は sample_sheet.csv に統一する。
単体サンプルでも複数サンプルでも同じ入口を使う。

### 10.3 内部的な考え方
内部では assay / matrix / execution run の単位で扱うが、README の主役にはしない。

### 10.4 いま大事にしていること
- Wet First を崩さない
- README に Spec を出しすぎない
- 人間向け入力は CSV、内部契約は別管理
- 次段アプリへ安全につながる出力を作る

## License

This repository is distributed under the Iwa Collections Non-Resale License 1.0.
Commercial resale of the software itself, or paid redistribution of derivative versions where the software is the primary value, is prohibited.

本リポジトリは Iwa Collections Non-Resale License 1.0 で公開しています。
ソフトウェア自体の有償販売、および本ソフトウェアが主たる価値となる派生物の有償再配布は禁止です。

