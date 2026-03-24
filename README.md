# iwa-rnaseq-counter

FASTQ から transcript / gene quant を出力するための、wet ファーストな RNA-Seq 解析アプリです。

## 主な機能
1. **FASTQ 自動検出**: 入力ディレクトリから R1/R2, Single-end, Lane 等を自動認識。
2. **strandedness 自動推定**: Salmon probe によるライブラリタイプの推定。
3. **実行前バリデーション**: チェックリストによるミス防止。
4. **成果物自動集約**: Transcript 単位および Gene 単位の TPM 値を出力。
5. **再現性**: 実行時の設定（run_config.json）とサンプルシートを自動保存。

## 動作確認済み条件 (v0.1.1)
- **OS**: Linux (WSL2 / Ubuntu 24.04 等)
- **Python**: 3.12+
- **Salmon**: v1.10.1 (推奨)
  - v1.11.x 以降ではインデックスの `versionInfo.json` が厳格チェックされるため注意が必要です。
- **データ**: Saccharomyces cerevisiae (Yeast) R64-1-1

## I/O Contract (v0.1.0)

`iwa-rnaseq-counter` は、RNA-Seq の解析を実行し、共通仕様（Spec）およびレガシー形式で成果物を出力する **Producer** です。

### 入力 (AssaySpec)
解析の入力条件は `AssaySpec` で定義されます。
- `fastq_r1`, `fastq_r2` パス
- `strandedness` (Auto-detect / unstranded / forward / reverse)
- `reference_resources` (Salmon index, tx2gene)

### 出力 (MatrixSpec & ExecutionRunSpec)
解析結果は以下の Spec 形式でメタデータと共に管理されます。
- **MatrixSpec**: 出力された遺伝子・転写産物発現行列の定義とパス。
- **ExecutionRunSpec**: アプリ名、バージョン、実行ステータス（completed/failed）、ログパス。

### レガシー成果物 (dataset contract) について
- **現行 Streamlit 実装**: 解析実行時に `dataset_manifest.json` などのレガシー成果物セットを生成・更新します。
- **CLI / Pipeline 内部**: 将来的な統合を見据え、現在は `MatrixSpec` / `ExecutionRunSpec` の生成を先行導入しています。

## 使い方
1. `pixi shell` で環境に入ります。
2. `streamlit run iwa_rnaseq_counter.py` で起動します。
3. 画面上の案内に従って以下のパスを入力します：
   - 入力ディレクトリ (例: `runs/SRA518891`)
   - Salmon index パス
   - tx2gene パス
4. **RUN START** をクリックして解析を開始します。

## 注意事項
- 現バージョンでは Salmon 解析のみをサポートしています。
- インデックス作成時には transcript ID と gene ID のマッピングファイル（tx2gene.csv）を事前に準備してください。

## License
This repository is distributed under the **Iwa Collections Non-Resale License 1.0**.
Commercial resale of the software itself, or paid redistribution of derivative versions where the software is the primary value, is prohibited.
本リポジトリは **Iwa Collections Non-Resale License 1.0** で公開しています。  
ソフトウェア自体の有償販売、および本ソフトウェアが主たる価値となる派生物の有償再配布は禁止です。
