from __future__ import annotations

import streamlit as st

from src.config import get_default_session_state
from src.fastq_discovery import collect_fastq_metadata, discover_fastq_files
from src.sample_parser import (
    apply_sample_table_edits,
    build_sample_table,
    detect_lane_groups,
    group_fastq_by_sample,
    infer_sample_layout,
)
from src.strandedness import infer_strandedness
from src.validators import (
    validate_input_directory,
    validate_output_directory,
    validate_run_conditions,
    validate_salmon_index,
    validate_tx2gene_file,
)
from src.salmon_runner import run_salmon_quant
from src.gene_aggregator import (
    load_tx2gene_map,
    aggregate_transcript_to_gene,
    build_transcript_quant_table,
    save_quant_tables,
)
from ui.sections import (
    render_analysis_section,
    render_app_header,
    render_fastq_section,
    render_reference_section,
    render_result_section,
    render_run_section,
    render_sample_section,
)


def init_session_state() -> None:
    defaults = get_default_session_state()
    for key, value in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = value


def run_app() -> None:
    st.set_page_config(page_title="iwa-rnaseq-counter", layout="wide")
    render_app_header()

    analysis_values = render_analysis_section()
    st.session_state.analysis_name = analysis_values["analysis_name"]
    st.session_state.input_dir = analysis_values["input_dir"]
    st.session_state.output_dir = analysis_values["output_dir"]

    input_validation = validate_input_directory(st.session_state.input_dir)
    output_validation = validate_output_directory(st.session_state.output_dir)

    if analysis_values["scan_requested"] and input_validation["is_valid"]:
        fastq_files = discover_fastq_files(st.session_state.input_dir)
        st.session_state.scanned_fastq_files = fastq_files
        st.session_state.fastq_df = collect_fastq_metadata(fastq_files)

        sample_groups = group_fastq_by_sample(fastq_files)
        layout_map = infer_sample_layout(sample_groups)
        lane_map = detect_lane_groups(sample_groups)
        st.session_state.sample_df = build_sample_table(sample_groups, layout_map, lane_map)

    render_fastq_section(st.session_state.fastq_df)

    sample_section_values = render_sample_section(st.session_state.sample_df)
    if st.session_state.sample_df is not None:
        st.session_state.sample_df = apply_sample_table_edits(
            st.session_state.sample_df,
            sample_section_values.get("edited_sample_df"),
        )

    reference_values = render_reference_section()
    st.session_state.salmon_index_path = reference_values["salmon_index_path"]
    st.session_state.tx2gene_path = reference_values["tx2gene_path"]

    index_validation = validate_salmon_index(st.session_state.salmon_index_path)
    tx2gene_validation = validate_tx2gene_file(st.session_state.tx2gene_path)

    strandedness_result = None
    if st.session_state.sample_df is not None and index_validation["is_valid"]:
        strandedness_result = infer_strandedness(
            sample_df=st.session_state.sample_df,
            input_dir=st.session_state.input_dir,
            salmon_index_path=st.session_state.salmon_index_path,
        )
        st.session_state.strandedness_prediction = strandedness_result

    run_validation = validate_run_conditions(
        input_dir=st.session_state.input_dir,
        output_dir=st.session_state.output_dir,
        sample_df=st.session_state.sample_df,
        salmon_index_path=st.session_state.salmon_index_path,
        tx2gene_path=st.session_state.tx2gene_path,
        strandedness_mode=st.session_state.strandedness_mode,
        strandedness_result=st.session_state.strandedness_prediction,
    )
    st.session_state.run_validation_result = run_validation

    run_values = render_run_section(
        strandedness_result=st.session_state.strandedness_prediction,
        run_validation=run_validation,
    )
    st.session_state.strandedness_mode = run_values["strandedness_mode"]

    if run_values["run_requested"]:
        with st.status("Running Salmon quantification...", expanded=True) as status:
            st.write("Step 1: Running Salmon...")
            run_result = run_salmon_quant(
                sample_df=st.session_state.sample_df,
                salmon_index_path=st.session_state.salmon_index_path,
                run_output_dir=st.session_state.output_dir,
                strandedness_mode=st.session_state.strandedness_mode,
            )
            st.session_state.run_status = "completed" if run_result["is_success"] else "failed"
            st.session_state.latest_log_summary = run_result["log_summary"]
            
            if run_result["is_success"]:
                st.write("Step 2: Aggregating results...")
                tx2gene_df = load_tx2gene_map(st.session_state.tx2gene_path)
                transcript_df = build_transcript_quant_table(run_result["outputs"])
                gene_df = aggregate_transcript_to_gene(run_result["outputs"], tx2gene_df)
                
                st.write("Step 3: Saving output files...")
                save_result = save_quant_tables(
                    transcript_df=transcript_df,
                    gene_df=gene_df,
                    run_output_dir=st.session_state.output_dir
                )
                
                st.session_state.output_files = [
                    {"file": "transcript_quant.csv", "path": save_result["transcript_quant_csv"]},
                    {"file": "gene_quant.csv", "path": save_result["gene_quant_csv"]},
                ]
                status.update(label="Analysis completed successfully!", state="complete", expanded=False)
            else:
                st.error("\n".join(run_result["errors"]))
                status.update(label="Analysis failed.", state="error", expanded=True)

    with st.expander("Validation details"):
        st.write({
            "input_validation": input_validation,
            "output_validation": output_validation,
            "index_validation": index_validation,
            "tx2gene_validation": tx2gene_validation,
            "run_validation": run_validation,
        })

    render_result_section(
        run_status=st.session_state.run_status,
        output_files=st.session_state.output_files,
        log_summary=st.session_state.latest_log_summary,
    )


def main() -> None:
    init_session_state()
    run_app()


if __name__ == "__main__":
    main()
