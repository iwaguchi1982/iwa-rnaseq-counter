from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd
import streamlit as st

from iwa_rnaseq_counter.legacy.config import get_strandedness_options
from iwa_rnaseq_counter.legacy.qc import evaluate_sample_qc
from iwa_rnaseq_counter.legacy.sample_parser import apply_sample_table_edits, parse_sample_sheet


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
        
        # Discovery Mode
        mode_options = {"auto": "FASTQ 自動認識 (input_dir)", "csv": "サンプルシート (CSV) 読込", "fastq": "FASTQ 直接指定 (1サンプル)"}
        current_mode = st.session_state.get("discovery_mode", "auto")
        idx = list(mode_options.keys()).index(current_mode) if current_mode in mode_options else 0
        discovery_mode = st.radio(
            "サンプル認識モード",
            options=list(mode_options.keys()),
            format_func=lambda x: mode_options[x],
            index=idx,
            horizontal=True
        )

    csv_path = ""
    load_csv_requested = False
    scan_requested = False
    fastq_direct_input = {}
    
    if discovery_mode == "csv":
        c1, c2 = st.columns([3, 1])
        with c1:
            csv_path = st.text_input("CSV パス (絶対パスまたは絶対パス推奨)", value=st.session_state.get("last_csv_path", ""))
        with c2:
            st.write("##") # alignment
            load_csv_requested = st.button("LOAD CSV", use_container_width=True)
        st.markdown("""
        **CSV形式の仕様:**
        - **必須:** `sample_id`, `r1_path`, `layout`
        - **条件付き必須:** `r2_path` (paired-end の場合)
        - **任意:** `group`, `condition`, `replicate`, `batch`, `pair_id`, `note`, `display_name`, `color`, `exclude`
        
        *※ r1_path, r2_path は入力ディレクトリからの相対パスまたは絶対パス。*
        
        *※ `exclude`: 今後の解析対象除外フラグ。現時点では保存のみで、quant 実行制御には未使用。*
        """)
        with st.expander("CSV サンプル例を表示"):
            st.code("sample_id,r1_path,r2_path,layout,group,condition,replicate,batch,pair_id,note,display_name,color,exclude\nSRR518891,runs/SRA518891/SRR518891_1.fastq,runs/SRA518891/SRR518891_2.fastq,paired-end,control,control,1,batch1,,yeast control replicate 1,Control-1,#4C78A8,false\nSRR518892,runs/SRA518892/SRR518892_1.fastq,runs/SRA518892/SRR518892_2.fastq,paired-end,treated,treated,1,batch1,,yeast treated replicate 1,Treated-1,#F58518,false\n", language="csv")
    elif discovery_mode == "fastq":
        st.write("### FASTQ 直接指定 (1サンプル)")
        f_c1, f_c2 = st.columns(2)
        with f_c1:
            fastq_direct_input["sample_id"] = st.text_input("Sample ID", value=st.session_state.get("fastq_sample_id", "sample_01"))
            fastq_direct_input["layout"] = st.selectbox("Layout", ["paired-end", "single-end"], index=0 if st.session_state.get("fastq_layout", "paired-end") == "paired-end" else 1)
            fastq_direct_input["exclude"] = st.checkbox("Exclude (解析対象から外す)", value=st.session_state.get("fastq_exclude", False))
        with f_c2:
            fastq_direct_input["r1_path"] = st.text_input("R1 FASTQ パス", value=st.session_state.get("fastq_r1_path", ""))
            fastq_direct_input["r2_path"] = st.text_input("R2 FASTQ パス (paired-end のみ)", value=st.session_state.get("fastq_r2_path", ""))
        
        load_fastq_requested = st.button("Set Sample", use_container_width=False)
        fastq_direct_input["load_requested"] = load_fastq_requested
    else:
        scan_requested = st.button("FASTQ をスキャン")
        
        if st.session_state.get("fastq_scan_completed") and st.session_state.get("fastq_scan_input_dir") == input_dir:
            st.info(
                f"✅ **スキャン完了**: {st.session_state.get('fastq_scan_timestamp')}\n\n"
                f"**対象**: `{st.session_state.get('fastq_scan_input_dir')}`\n\n"
                f"**検出**: {st.session_state.get('fastq_scan_file_count')} 個の FASTQ ファイル"
            )

    # Sync to session state
    st.session_state.analysis_name = analysis_name
    st.session_state.input_dir = input_dir
    st.session_state.output_dir = output_dir
    st.session_state.discovery_mode = discovery_mode
    if csv_path: st.session_state.last_csv_path = csv_path
    if discovery_mode == "fastq":
        st.session_state.fastq_sample_id = fastq_direct_input["sample_id"]
        st.session_state.fastq_layout = fastq_direct_input["layout"]
        st.session_state.fastq_r1_path = fastq_direct_input["r1_path"]
        st.session_state.fastq_r2_path = fastq_direct_input["r2_path"]
        st.session_state.fastq_exclude = fastq_direct_input["exclude"]

    return {
        "analysis_name": analysis_name.strip(),
        "input_dir": input_dir.strip(),
        "output_dir": output_dir.strip(),
        "discovery_mode": discovery_mode,
        "csv_path": csv_path.strip(),
        "load_csv_requested": load_csv_requested,
        "scan_requested": scan_requested,
        "fastq_direct_input": fastq_direct_input,
    }


def render_fastq_section(fastq_df: pd.DataFrame | None) -> None:
    st.subheader("2. 入力 FASTQ")
    if fastq_df is None or fastq_df.empty:
        st.warning("FASTQ はまだ検出されていません。")
        return
    st.dataframe(fastq_df, use_container_width=True)


def render_sample_section(sample_df: pd.DataFrame) -> pd.DataFrame:
    st.subheader("3. サンプル確認・編集")
    if sample_df.empty:
        st.info("サンプルが検出されていません。input_dir を確認するか、CSV を読み込んでください。")
        return sample_df

    st.caption("必要に応じて Sample ID や Layout, Group を編集してください。")
    
    # セッション状態での編集保持
    edited_df = st.data_editor(
        sample_df,
        column_config={
            "sample_id": st.column_config.TextColumn("Sample ID", required=True),
            "group": st.column_config.TextColumn("Group", help="control, treated 等"),
            "condition": st.column_config.TextColumn("Condition", help="drugA, wild-type 等"),
            "replicate": st.column_config.TextColumn("Replicate", help="1, 2, R1 等"),
            "batch": st.column_config.TextColumn("Batch", help="batch1, plateA 等"),
            "pair_id": st.column_config.TextColumn("PairID", help="Pair 解析用 ID"),
            "note": st.column_config.TextColumn("Note", help="備考"),
            "display_name": st.column_config.TextColumn("Display Name", help="レポート用表示名"),
            "color": st.column_config.TextColumn("Color", help="プロット用色 (例: #FF0000, red)"),
            "exclude": st.column_config.CheckboxColumn("Exclude", help="今後の解析対象除外フラグ（現時点では保存のみ）"),
            "layout_final": st.column_config.SelectboxColumn(
                "Layout",
                options=["paired-end", "single-end"],
                required=True,
            ),
            "status": st.column_config.TextColumn("Status", disabled=True),
            "notes": st.column_config.TextColumn("Notes", disabled=True),
            "r1_paths": None,  # 非表示
            "r2_paths": None,
            "all_paths": None,
            "layout_predicted": None,
            "lane_count": None,
            "file_count": None,
            "input_source": None,
        },
        hide_index=True,
        use_container_width=True,
        key="sample_editor"
    )
    
    return apply_sample_table_edits(sample_df, edited_df)


def render_reference_section(quantifier: str) -> dict[str, Any]:
    st.subheader("4. 参照設定")

    q = str(quantifier or "salmon").strip().lower()

    quantifier_index_path = st.text_input(
        "Quantifier index パス",
        value=st.session_state.get("quantifier_index_path", ""),
        help="選択中 backend の index / prefix / index file を指定します。",
    )

    tx2gene_required = q in {"salmon", "kallisto"}
    annotation_gtf_required = q == "hisat2"

    if q == "salmon":
        st.caption("Salmon: quantifier index と tx2gene が必要です。")
    elif q == "star":
        st.caption("STAR: genome index が必要です。tx2gene は任意ですが、annotation 用に保持できます。")
    elif q == "hisat2":
        st.caption("HISAT2: index prefix と annotation GTF が必要です。")
    elif q == "kallisto":
        st.caption("kallisto: index file と tx2gene が必要です。")

    tx2gene_path = st.text_input(
        f"tx2gene パス{' *' if tx2gene_required else ''}",
        value=st.session_state.get("tx2gene_path", ""),
        help="transcript-level backend から gene-level へ集約するときに使用します。",
    )

    annotation_gtf_path = st.text_input(
        f"Annotation GTF パス{' *' if annotation_gtf_required else ''}",
        value=st.session_state.get("annotation_gtf_path", ""),
        help="HISAT2 gene-level counting で使用します。",
        disabled=not annotation_gtf_required,
    )

    st.write("---")
    st.write("💡 Strandedness (鎖特異性) の自動推定")

    if q == "salmon":
        st.caption("サンプル一覧と quantifier index パスが揃っている場合、データをサンプリングして strandedness を推定できます。")
        estimate_strandedness = st.button("Strandedness を推定", use_container_width=False)
    else:
        st.caption(f"v0.8.0 時点では strandedness 自動推定は {q if q != 'salmon' else ''} backend 未対応です。")
        estimate_strandedness = False

    return {
        "quantifier_index_path": quantifier_index_path.strip(),
        "tx2gene_path": tx2gene_path.strip(),
        "annotation_gtf_path": annotation_gtf_path.strip(),
        "estimate_strandedness": estimate_strandedness,
    }


def render_run_section(
    strandedness_result: dict[str, Any] | None = None,
    run_validation: dict[str, Any] | None = None,
) -> dict[str, Any]:
    st.subheader("5. 実行設定")
    col1, col2 = st.columns(2)
    with col1:
        options = get_strandedness_options()
        default_index = options.index(st.session_state.strandedness_mode) if st.session_state.strandedness_mode in options else 0
        strandedness_mode = st.selectbox("strandedness", options=options, index=default_index)

    with col2:
        threads = st.number_input("Quantifier スレッド数", min_value=1, max_value=64, value=st.session_state.threads)

    if strandedness_result:
        color = "green" if strandedness_result.get("confidence") == "high" else "orange"
        st.markdown(
            f"**推定結果:** :{color}[{strandedness_result.get('mode', 'unknown')}] "
            f"(信頼度: {strandedness_result.get('confidence', 'unknown')})\n\n"
            f"*{strandedness_result.get('reason', '')}*"
        )
        if strandedness_result.get("probe_cmd"):
            with st.expander("Probe 詳細 (Salmon)"):
                st.code(" ".join(strandedness_result["probe_cmd"]), language="bash")
                if strandedness_result.get("stderr"):
                    st.text_area("Probe Stderr (Detailed)", value=strandedness_result["stderr"], height=300, key="probe_stderr")

    can_run = bool(run_validation and run_validation.get("is_valid", False))
    
    with st.container(border=True):
        st.write("### 実行前チェック")
        if run_validation:
            checks = run_validation.get("checks", {})
            c1, c2, c3 = st.columns(3)
            with c1:
                st.write(f"{'✅' if checks.get('fastq_detected') else '❌'} 入力サンプル検出")
                st.write(f"{'✅' if checks.get('sample_structure') else '❌'} サンプル不整合なし")
            with c2:
                st.write(f"{'✅' if checks.get('quantifier_index') else '❌'} Quantifier index")
                st.write(f"{'✅' if checks.get('backend_references') else '❌'} Backend reference")
            with c3:
                st.write(f"{'✅' if checks.get('strandedness') else '❌'} strandedness")
                st.write(f"{'✅' if checks.get('output_dir') else '❌'} 出力先")
            
            st.write(f"{'✅' if checks.get('quantifier_binary') else '❌'} Binary 依存関係")

            if not can_run:
                st.error("実行前に上記の ❌ 項目を修正してください。")
                for error in run_validation.get("errors", []):
                    st.write(f"- {error}")
            else:
                st.success("すべての準備が整いました。")

        run_requested = st.button("RUN START", type="primary", disabled=not can_run, use_container_width=True)

    return {
        "strandedness_mode": strandedness_mode,
        "threads": threads,
        "run_requested": run_requested,
    }


def render_result_section(
    run_status: str,
    output_files: list[dict[str, str]] | None = None,
    log_summary: str | None = None,
    run_summary: dict[str, Any] | None = None,
) -> None:
    st.subheader("6. 実行・結果")
    st.caption("v0.1.7: Reporter 連携用出力標準化（manifest 生成）対応。")

    if run_status == "idle":
        st.info("RUN ボタンを押すと解析が開始されます。")
        return

    if run_status == "running":
        st.warning("解析実行中...")
        return

    if run_status == "success" and run_summary:
        sample_count = run_summary.get('sample_count', 0)
        save_path = run_summary.get('save_path', 'N/A')
        analysis_name = run_summary.get('analysis_name', 'N/A')
        st.success(f"解析 `{analysis_name}` が完了しました。{sample_count} サンプルの結果は `{save_path}` に保存されました。")

        # 1. Dashboard metrics
        m1, m2, m3, m4 = st.columns(4)
        m1.metric("サンプル数", f"{sample_count} samples")
        m2.metric("遺伝子（集約数）", f"{run_summary.get('gene_rows', 0):,}")
        m3.metric("転写産物数", f"{run_summary.get('transcript_rows', 0):,}")
        elapsed = run_summary.get('elapsed_seconds', 0)
        m4.metric("実行時間", f"{elapsed:.1f}s" if elapsed < 60 else f"{elapsed/60:.1f}m")

        # 2. Reference Info & QC Overview
        st.write("---")
        c_ref, c_tbl = st.columns([1, 2])
        
        with c_ref:
            with st.container(border=True):
                st.markdown("#### 参照設定・環境")
                st.write(f"**Quantifier:** `{run_summary.get('quantifier', 'salmon')}`")
                st.write(f"**Quantifier Version:** `{run_summary.get('quantifier_version', 'N/A')}`")
                
                stranded_info = run_summary.get('strandedness') or {}
                st.write(f"**Inferred Lib Type:** `{stranded_info.get('mode', 'N/A')}`")
                
                index_name = Path(
                    run_summary.get("quantifier_index_path")
                    or run_summary.get("salmon_index_path", "")
                ).name
                tx2gene_name = Path(run_summary.get("tx2gene_path", "")).name
                annotation_gtf_path = run_summary.get("annotation_gtf_path", "")
                annotation_gtf_name = Path(annotation_gtf_path).name if annotation_gtf_path else ""

                st.write(f"**Index:** `{index_name or 'N/A'}`")
                st.write(f"**tx2gene:** `{tx2gene_name or 'N/A'}`")
                if annotation_gtf_name:
                    st.write(f"**Annotation GTF:** `{annotation_gtf_name}`")

        with c_tbl:
            if "outputs" in run_summary:
                with st.expander("マッピング統計サマリ", expanded=True):
                    data = []
                    for o in run_summary["outputs"]:
                        if not o.get("is_success"): continue
                        
                        qc = evaluate_sample_qc(
                            sample_id=o["sample_id"],
                            num_processed=o.get("num_processed", 0),
                            num_mapped=o.get("num_mapped", 0),
                            num_decoy=o.get("num_decoy", 0),
                            num_filter=o.get("num_filter", 0)
                        )
                        
                        data.append({
                            "QC": qc["icon"],
                            "Sample": o["sample_id"],
                            "Mapped%": f"{qc['mapping_rate']:.2f}%",
                            "Decoy%": f"{qc['decoy_rate']:.1f}%",
                            "Filt%": f"{qc['filter_rate']:.1f}%",
                            "Mapped Reads": f"{o.get('num_mapped', 0):,}",
                            "Alerts": ", ".join(qc["alerts"]) if qc["alerts"] else "OK"
                        })
                    st.table(pd.DataFrame(data))
                    st.caption("🔴:<1%, 🟡:<10% mapping rate | Decoy: >20% | Filter: >30%")

        # 3. Detailed Sample Cards
        if "outputs" in run_summary:
            st.write("#### サンプル別詳細")
            all_outputs = run_summary["outputs"]
            
            # Show 3 cards per row
            for i in range(0, len(all_outputs), 3):
                cols = st.columns(3)
                for j in range(3):
                    if i + j < len(all_outputs):
                        o = all_outputs[i + j]
                        if not o.get("is_success"):
                            with cols[j]:
                                with st.container(border=True):
                                    st.markdown(f"#### ❌ {o['sample_id']}")
                                    st.error(f"Error: {o.get('error_reason', 'Unknown error')}")
                            continue

                        qc = evaluate_sample_qc(
                            sample_id=o["sample_id"],
                            num_processed=o.get("num_processed", 0),
                            num_mapped=o.get("num_mapped", 0),
                            num_decoy=o.get("num_decoy", 0),
                            num_filter=o.get("num_filter", 0)
                        )
                        
                        with cols[j]:
                            with st.container(border=True):
                                st.markdown(f"#### {qc['icon']} {o['sample_id']}")
                                
                                # Metrics column
                                m_c1, m_c2 = st.columns(2)
                                with m_c1:
                                    st.caption("**Mapped:**")
                                    st.caption("**Decoy:**")
                                    st.caption("**Filtered:**")
                                with m_c2:
                                    st.caption(f"`{qc['mapping_rate']:.2f}%`")
                                    st.caption(f"`{qc['decoy_rate']:.1f}%`")
                                    st.caption(f"`{qc['filter_rate']:.1f}%`")
                                
                                # Alerts
                                if qc["alerts"]:
                                    for alert in qc["alerts"]:
                                        st.caption(f":orange[• {alert}]")
                                
                                # Files
                                st.write("---")
                                qpath = o.get("quant_path", "")
                                if qpath:
                                    q_p = Path(qpath)
                                    log_path = q_p.parent / "logs" / "salmon_quant.log"
                                    st.markdown(f"📄 [quant.sf](file://{qpath})")
                                    st.markdown(f"📋 [salmon_log](file://{log_path})")

        st.write("---")

    if output_files:
        st.write("### 成果物ダウンロード (results/)")
        df_out = pd.DataFrame(output_files)
        df_display = df_out.rename(columns={"label": "名称", "path": "フルパス"})
        st.dataframe(df_display, use_container_width=True, hide_index=True)

    if log_summary:
        with st.expander("Master Execution Log (Summary)", expanded=True):
            st.text_area("Log summary", value=log_summary, height=400, label_visibility="collapsed", key="run_log_summary")
def format_elapsed_time(seconds: float | None, status: str) -> str:
    """Formats elapsed time based on status for UI display."""
    if status == "queued":
        return "waiting..."
    if seconds is None or seconds <= 0:
        return "--"
    
    m, s = divmod(int(seconds), 60)
    h, m = divmod(m, 60)
    if h > 0:
        return f"{h}h {m}m"
    if m > 0:
        return f"{m}m {s}s"
    return f"{s}s"


def render_job_history_panel(
    recent_runs: list[Any], 
    selected_job_id: str | None
) -> str | None:
    """
    Renders the job history panel (formerly sidebar) and returns the ID of the selected job.
    """
    st.markdown("### 🕒 Job History")
    
    if not recent_runs:
        st.caption("No valid runs found in this root.")
        return selected_job_id

    # List navigator
    new_selection = selected_job_id
    
    for run in recent_runs:
        job_id = run.run_dir.name
        status = (run.status or "unknown").lower()
        
        # Status Icon selection
        status_icon = "❔"
        if status == "running":
            status_icon = "🔵"
        elif status == "queued":
            status_icon = "⏳"
        elif status == "completed":
            status_icon = "✅"
        elif status == "failed":
            status_icon = "❌"
            
        elapsed_str = format_elapsed_time(run.elapsed_seconds, status)
        
        # Highlight if selected
        is_selected = (selected_job_id == job_id)
        
        # Label: [ICON] RUN_NAME
        label = f"{status_icon} {run.run_name}"
        
        col_name, col_time = st.columns([3, 1])
        with col_name:
            if st.button(
                label, 
                key=f"history_job_{job_id}", 
                help=f"ID: {job_id}",
                use_container_width=True,
                type="primary" if is_selected else "secondary"
            ):
                new_selection = job_id
            
            # Show STATUS explicitly as text
            status_display = f"STATUS: {status}"
            if is_selected:
                st.caption(f"**{status_display}**  •  selected")
            else:
                st.caption(status_display)
                
        with col_time:
            st.caption(elapsed_str)
            
    return new_selection
