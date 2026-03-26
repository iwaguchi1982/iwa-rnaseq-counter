# Developer_README.md

## 目的

この文書は `iwa-rnaseq-counter` の開発者向けメモです。  
利用者向け README ではなく、将来の自分が設計意図と現在の実装位置を思い出すための短い整理です。

---

## 1. このアプリの役割

`iwa-rnaseq-counter` は RNA-Seq Suite における **count 作成の入口** です。

表向きには、Wet 利用者が

- FASTQ を用意する
- サンプル情報をそろえる
- 実行する
- count を得る

ためのアプリです。

裏側では、

- assay 単位の入力を正規化する
- count 結果を後段へ安全に渡せる形にする
- 実行記録を残す
- reporter へ接続できる形を作る

ことを担います。

---

## 2. Wet First の原則

このアプリは最初から **Wet 主体利用** を前提にしています。  
「最終的に Wet も使える」ではなく、最初から Wet の人が使えることが前提です。

そのため、

- README は GUI 説明を先頭に置く
- CLI も Wet な Dry が理解できる入口にする
- 内部 Spec を利用者向け説明の主役にしない
- JSON を人に直接書かせる方向へ寄せすぎない

ことを重視します。

---

## 3. 入力の考え方

### 3.1 人間向け入口
人間向けの入口は `sample_sheet.csv` に統一する方針です。

- 単体サンプルでも sample sheet
- 複数サンプルでも sample sheet
- 1行1assay

### 3.2 内部表現
内部では必要に応じて以下へ正規化します。

- AssaySpec
- MatrixSpec
- ExecutionRunSpec

つまり、

- 人間向け入力 = CSV
- 内部契約 = Spec

の二層構造です。

---

## 4. sample sheet の現在方針

現在の基本列は以下です。

- sample_id
- r1_path
- r2_path
- layout
- group
- condition
- replicate
- batch
- pair_id
- display_name
- color
- exclude
- note

### single / paired の扱い
- `layout=single` のとき `r2_path` は空欄
- `layout=paired` のとき `r1_path`, `r2_path` を両方記入

### 内部変換
- `sample_id` は人間向け名称として使う
- 必要に応じて内部では `specimen_id` 相当へ写像する

---

## 5. 現在の主要機能

### 5.1 Web UI
主役は Streamlit UI。  
Wet 利用者はまずこちらを使う前提。

### 5.2 CLI
CLI は補助経路。  
主として以下を想定する。

- `run-batch --sample-sheet ...`
- 将来必要なら内部確認用コマンドを別途維持

### 5.3 出力
主な出力は以下。

- count table
- matrix artifact
- 実行記録
- 次段 reporter に渡せる形式

---

## 6. reporter との関係

`iwa-rnaseq-counter` は、最終的に `iwa-rnaseq-reporter` が読める結果を作る入口です。

意識する接続線は次の通りです。

1. count を作る
2. 必要なら複数サンプル結果を統合する
3. reporter 側で読み込める形をそろえる

大事なのは、counter が reporter の責務まで侵食しないことです。

- counter は結果を作る
- reporter は結果を見て、比較準備を進める

---

## 7. いまの設計上の芯

### 7.1 README に Spec を出しすぎない
内部では Spec-aware でも、利用者向け README では出しすぎない。

### 7.2 人間向け入口は CSV
人に JSON を直接書かせない。  
JSON / Spec は内部正規化用。

### 7.3 Wet の人に Dry 手続きを押しつけない
- データ契約
- provenance
- execution record
- 内部構造

は裏側で支える。  
利用者に意識させない。

### 7.4 日本標準の Dry を前提にする
世界標準の専業 Dry ではなく、  
「少し PC に覚えのある Wet 寄りの人」が扱うことを前提にする。

---

## 8. v0.3.0 での達成事項と今後の拡張メモ

### v0.3.0 での達成
- `run-batch --sample-sheet` による入口の一本化
- `sample_sheet.csv` からの AssaySpec 群の生成と安定化
- `iwa-rnaseq-reporter` へ安全に繋がる unified matrix の出力と接続確認

### 今後の拡張メモ
- paired design の本格対応
- covariates の本格対応
- output 構造のさらなる統一

---

## 9. 自分向けの確認ポイント

迷ったら次を確認する。

1. このファイル
2. 利用者向け README
3. sample sheet の example
4. 実装コード
5. CONTRACTS / specs 系文書

優先順位は常に、

- Wet 利用しやすさ
- 次段へ安全につながること
- 内部契約の破綻を避けること

の順で考える。
