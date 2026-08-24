"""Streamlit dashboard surfacing Part 2, Part 3, and Part 4 results."""

import sys
from pathlib import Path

import pandas as pd
import streamlit as st

sys.path.insert(0, str(Path(__file__).parent.parent))

import load_data
from analysis import db
from analysis.part2_frequencies import compute_frequencies
from analysis.part3_stats import get_melanoma_pbmc_miraclib, run_stats, make_boxplots, summarize
from analysis.part4_subset import (
    get_baseline_subset,
    by_project,
    by_response,
    by_sex,
    average_bcell_male_responders,
)

st.set_page_config(page_title="Teiko Cell Count Analysis", layout="wide")
st.title("Immune Cell Count Analysis")

if not db.DB_PATH.exists():
    with st.spinner("First run: building cell_counts.db from cell-count.csv..."):
        load_data.main()

conn = db.get_connection()

tab2, tab3, tab4 = st.tabs(["Part 2: Frequencies", "Part 3: Responder Comparison", "Part 4: Baseline Subset"])

with tab2:
    st.header("Relative frequency per sample x population")
    freq_df = compute_frequencies(conn)

    samples = sorted(freq_df["sample"].unique())
    selected_samples = st.multiselect("Filter by sample (leave empty for all)", samples)
    populations = sorted(freq_df["population"].unique())
    selected_pops = st.multiselect("Filter by population (leave empty for all)", populations)

    display_df = freq_df
    if selected_samples:
        display_df = display_df[display_df["sample"].isin(selected_samples)]
    if selected_pops:
        display_df = display_df[display_df["population"].isin(selected_pops)]

    st.dataframe(display_df, width='stretch', height=400)
    st.caption(f"{len(display_df)} of {len(freq_df)} rows shown")

with tab3:
    st.header("Responders vs non-responders (melanoma, PBMC, miraclib)")
    merged = get_melanoma_pbmc_miraclib(conn)
    stats_df = run_stats(merged)

    boxplot_path = db.OUTPUTS_DIR / "part3_boxplots.png"
    if not boxplot_path.exists():
        db.ensure_outputs_dir()
        make_boxplots(merged, boxplot_path)
    st.image(str(boxplot_path), width='stretch')

    st.subheader("Mann-Whitney U test per population")
    st.caption(
        "Mann-Whitney U was chosen over Welch's t-test because relative "
        "frequencies are bounded proportions with no guarantee of normality; "
        "the rank-based test is robust to that without assuming a normal distribution."
    )
    st.dataframe(stats_df, width='stretch')

    st.subheader("Conclusion")
    st.markdown(summarize(stats_df).replace("\n", "\n\n"))

with tab4:
    st.header("Baseline (time=0) melanoma / PBMC / miraclib subset")
    subset = get_baseline_subset(conn)

    col1, col2, col3 = st.columns(3)
    with col1:
        st.subheader("By project")
        st.dataframe(by_project(subset), width='stretch')
    with col2:
        st.subheader("By response")
        st.dataframe(by_response(subset), width='stretch')
    with col3:
        st.subheader("By sex")
        st.dataframe(by_sex(subset), width='stretch')

    st.metric("Total samples in baseline subset", len(subset))

    st.subheader("Average B cell count: melanoma male responders, time=0")
    st.caption("Across all sample types and treatments (not restricted to PBMC/miraclib).")
    n_rows, avg_bcell = average_bcell_male_responders(conn)
    st.metric("Average B cell count", f"{avg_bcell:.2f}", help=f"n={n_rows} matching sample rows")

conn.close()
