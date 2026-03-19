# iwa-rnaseq-counter

wet の人が、まず FASTQ から transcript / gene quant までたどり着くための RNA-Seq 入口アプリです。

## v0.1.0 scope

- 初期はSalmon のみ対応
  - 最終的にはSTAR, HISAT2, kallisto迄は対応予定
- FASTQ 検出
- single-end / paired-end 推定
- lane 論理統合
- strandedness 自動推定の表示
- transcript / gene quant CSV の出力
- 実行ログ / 設定保存

## Tech stack

- Python 3.11+
- Streamlit
- pandas
- Pixi
- Salmon (external command)

## Project layout

```text
.
├── iwa_rnaseq_counter.py
├── pixi.toml
├── ui/
│   └── sections.py
├── src/
│   ├── config.py
│   ├── fastq_discovery.py
│   ├── gene_aggregator.py
│   ├── salmon_runner.py
│   ├── sample_parser.py
│   ├── strandedness.py
│   └── validators.py
└── docs/
```

## Quick start

```bash
pixi install
pixi run streamlit run iwa_rnaseq_counter.py
```

## Notes

- v0.1.0 では既存 Salmon index の利用のみを対象にします。
- index 作成、STAR + featureCounts、kallisto は将来拡張です。
- 現在のコードは GitHub に上げて育てていくための最小雛形です。

## License

Repository owner intends to use: **Iwa Collections Non-Resale License 1.0**.

The full license text has not been finalized in this scaffold yet. Replace `LICENSE` with the final text before public release.
