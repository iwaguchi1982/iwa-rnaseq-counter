from __future__ import annotations

from typing import Any

import pandas as pd
import streamlit as st

from src.config import get_strandedness_options


def render_app_header() -> None:
    st.title("iwa-rnaseq-counter")
    st.markdown(
        "FASTQ から transcript / gene quant を出力するための、"
        "wet ファーストな RNA-Seq 入口アプリです。"
    )


def render_analysis_section() -> dict[str, Any]:
    st.subheader("1. 解析情報")
    col1, col2 = st.columns(2)
    with col1:
        analysis_name = st.text_input("解析名", value=st.session_state.analysis_name)
        input_dir = st.text_input("入力ディレクトリ", value=st.session_state.input_dir)
    with col2:
        output_dir = st.text_input("出力ディレクトリ", value=st.session_state.output_dir)
        scan_requested = st.button("FASTQ をスキャン")
    return {
        "analysis_name": analysis_name,
        "input_dir": input_dir,
        "output_dir": output_dir,
        "scan_requested": scan_requested,
    }


def render_fastq_section(fastq_df: pd.DataFrame | None) -> None:
    st.subheader("2. 入力 FASTQ")
    if fastq_df is None or fastq_df.empty:
        st.warning("FASTQ はまだ検出されていません。")
        return
    st.dataframe(fastq_df, use_container_width=True)


def render_sample_section(sample_df: pd.DataFrame | None) -> dict[str, Any]:
    st.subheader("3. サンプル認識結果")
    if sample_df is None or sample_df.empty:
        st.info("FASTQ をスキャンするとサンプル認識結果が表示されます。")
        return {"edited_sample_df": None}

    st.info("lane は物理結合せず、論理的に同一サンプルとして扱います。")
    edited_sample_df = st.data_editor(sample_df, use_container_width=True, disabled=["lane_count", "file_count"])
    return {"edited_sample_df": edited_sample_df}


def render_reference_section() -> dict[str, str]:
    st.subheader("4. 参照設定")
    salmon_index_path = st.text_input("Salmon index パス", value=st.session_state.salmon_index_path)
    tx2gene_path = st.text_input("tx2gene パス", value=st.session_state.tx2gene_path)
    st.caption("STAR + featureCounts / kallisto は将来拡張です。v0.1.0 では Salmon のみ有効です。")
    return {
        "salmon_index_path": salmon_index_path,
        "tx2gene_path": tx2gene_path,
    }


def render_run_section(
    strandedness_result: dict[str, Any] | None = None,
    run_validation: dict[str, Any] | None = None,
) -> dict[str, Any]:
    st.subheader("5. 実行設定")
    options = get_strandedness_options()
    default_index = options.index(st.session_state.strandedness_mode) if st.session_state.strandedness_mode in options else 0
    strandedness_mode = st.selectbox("strandedness", options=options, index=default_index)

    if strandedness_result:
        st.info(
            f"推定結果: {strandedness_result.get('mode', 'unknown')} / "
            f"confidence: {strandedness_result.get('confidence', 'unknown')} / "
            f"reason: {strandedness_result.get('reason', '')}"
        )

    can_run = bool(run_validation and run_validation.get("is_valid", False))
    if run_validation:
        for error in run_validation.get("errors", []):
            st.error(error)
        for warning in run_validation.get("warnings", []):
            st.warning(warning)
        if can_run:
            st.success("実行可能です。")

    run_requested = st.button("RUN", disabled=not can_run)
    return {
        "strandedness_mode": strandedness_mode,
        "run_requested": run_requested,
    }


def render_result_section(
    run_status: str,
    output_files: list[dict[str, str]] | None = None,
    log_summary: str | None = None,
) -> None:
    st.subheader("6. 実行・結果")
    st.write(f"実行状態: {run_status}")
    if log_summary:
        st.text_area("ログ要約", value=log_summary, height=180)
    if output_files:
        st.dataframe(pd.DataFrame(output_files), use_container_width=True)
