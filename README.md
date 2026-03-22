# iwa-rnaseq-counter

FASTQ から transcript / gene quant を出力するための、wet ファーストな RNA-Seq 解析アプリです。

## 動作確認済み条件 (v0.1.1)
- **OS**: Linux (WSL2 / Ubuntu 24.04 等)
- **Python**: 3.12+
- **Salmon**: v1.10.1 (推奨)
  - v1.11.x 以降ではインデックスの `versionInfo.json` が厳格チェックされるため注意が必要です。
- **データ**: Saccharomyces cerevisiae (Yeast) R64-1-1

<<<<<<< HEAD
- 初期はSalmon のみ対応
  - 最終的にはSTAR, HISAT2, kallisto迄は対応予定
- FASTQ 検出
- single-end / paired-end 推定
- lane 論理統合
- strandedness 自動推定の表示
- transcript / gene quant CSV の出力
- 実行ログ / 設定保存
=======
## 主な機能
1. **FASTQ 自動検出**: 入力ディレクトリから R1/R2, Single-end, Lane 等を自動認識。
2. **strandedness 自動推定**: Salmon probe によるライブラリタイプの推定。
3. **実行前バリデーション**: チェックリストによるミス防止。
4. **成果物自動集約**: Transcript 単位および Gene 単位の TPM 値を出力。
5. **再現性**: 実行時の設定（run_config.json）とサンプルシートを自動保存。
>>>>>>> e88e41a (feat: v0.1.5 multi-sample and sample-sheet execution baseline)

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
