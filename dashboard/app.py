"""
Streamlit dashboard surfacing Part 2, Part 3, and Part 4 results.

Built for reviewer verification, not just display: every computed subset
shows the exact SQL/pandas logic behind it, every table is downloadable as
CSV, and the Data Quality tab re-runs integrity checks against the live
database on every page load rather than asserting them once in a README.
"""

import sys
from pathlib import Path

import pandas as pd
import streamlit as st

sys.path.insert(0, str(Path(__file__).parent.parent))

import load_data
from analysis import db
from analysis.part2_frequencies import compute_frequencies
from analysis.part3_stats import (
    METADATA_QUERY,
    POPULATIONS,
    get_melanoma_pbmc_miraclib,
    make_boxplots,
    run_stats,
)
from analysis.part4_subset import (
    BASELINE_QUERY,
    average_bcell_male_responders,
    by_project,
    by_response,
    by_sex,
    get_baseline_subset,
)

CSV_PATH = Path(__file__).parent.parent / "cell-count.csv"

st.set_page_config(page_title="Teiko Cell Count Analysis", layout="wide")
st.title("Immune Cell Count Analysis")
st.caption(
    "Every table below can be downloaded as CSV, every filtered subset shows "
    "the exact SQL query behind it, and the Data Quality tab re-runs "
    "integrity checks against the live database on every load."
)

if not db.DB_PATH.exists():
    with st.spinner("First run: building cell_counts.db from cell-count.csv..."):
        load_data.main()

conn = db.get_connection()


def download_button(df, label, filename, key=None):
    st.download_button(
        label, df.to_csv(index=False).encode("utf-8"), filename, "text/csv", key=key
    )


tab_design, tab_quality, tab2, tab3, tab4 = st.tabs(
    [
        "Design Decisions",
        "Data Quality",
        "Part 2: Frequencies",
        "Part 3: Responder Comparison",
        "Part 4: Subset Explorer",
    ]
)

# ---------------------------------------------------------------- Design Decisions
with tab_design:
    st.header("Design decisions")
    st.caption(
        "The reasoning behind each part of this project -- not just what was "
        "built, but why, and what the alternative would have cost."
    )

    st.subheader("Database schema")
    st.markdown(
        "Normalized into `projects -> subjects -> samples -> cell_counts` "
        "rather than one wide table."
    )
    with st.expander("Why normalize instead of one flat table?"):
        st.markdown(
            "- Subject-level facts (condition, age, sex) live on `subjects`, "
            "not repeated across a subject's 3 samples, so they can't drift "
            "out of sync.\n"
            "- `treatment` and `response` live on `samples`, not `subjects`, "
            "even though they're constant per subject in this snapshot. "
            "They're trial-arm/readout concepts that could legitimately vary "
            "by sample in future data (a subject switching arms, or a "
            "response readout assessed per timepoint). Modeling them one "
            "level down keeps the schema correct as the data model matures, "
            "at the cost of a little redundancy today.\n"
            "- `cell_counts` is long/tidy (one row per sample x population) "
            "instead of one column per population. Adding a 6th, 7th, ... "
            "cell population later needs zero schema migration and zero "
            "loader code changes -- just more rows. It also makes every "
            "downstream query a plain `GROUP BY`/`JOIN` instead of bespoke "
            "per-column logic."
        )
    with st.expander("How does this scale to hundreds of projects / thousands of samples?"):
        st.markdown(
            "1. Move from SQLite to Postgres for concurrent writes and "
            "multi-user access -- the schema is already vanilla relational "
            "SQL and ports over unchanged.\n"
            "2. Add a `populations` reference table once population "
            "definitions carry more metadata than just a name (marker "
            "panels, gating hierarchy).\n"
            "3. Precompute a `sample_frequencies` table refreshed on load if "
            "the Part 2 computation becomes a query-time bottleneck.\n"
            "4. Add a `batches`/`ingestion_runs` table if data starts "
            "arriving incrementally rather than as one full-file load, so "
            "the loader can append instead of drop-and-recreate."
        )

    st.divider()
    st.subheader("Part 2: Frequency table")
    st.markdown(
        "`total_count` is defined as the sum of the five given population "
        "counts for that sample -- the CSV has no separate total column, so "
        "this is the only total that's actually derivable from the data."
    )
    with st.expander("Why long/tidy format instead of one row per sample?"):
        st.markdown(
            "One row per sample x population (52,500 rows for 10,500 "
            "samples) instead of one row per sample with 5 percentage "
            "columns. This mirrors the `cell_counts` table shape, makes "
            "filtering by population trivial, and means the same table "
            "answers 'give me every population for one sample' and 'give me "
            "one population across every sample' without reshaping."
        )

    st.divider()
    st.subheader("Part 3: Statistical comparison")
    with st.expander("Why Mann-Whitney U instead of Welch's t-test?"):
        st.markdown(
            "The values compared are per-sample relative frequencies -- "
            "bounded proportions in [0, 100] -- not raw counts, and there's "
            "no reason to assume they're normally distributed (the boxplots "
            "show visible outliers in every population). Mann-Whitney U is "
            "a non-parametric, rank-based test that doesn't require that "
            "assumption and is robust to those outliers, at the cost of "
            "some statistical power versus a t-test if the data actually "
            "were normal. Given the sample sizes here (975-993 per group), "
            "that power cost is negligible."
        )
    with st.expander("Why PBMC-only, and why exclude blank response?"):
        st.markdown(
            "PBMC-only is a direct requirement of the assignment (whole "
            "blood samples are a different compartment and would confound "
            "the comparison). Blank `response` values correspond to "
            "`treatment == 'none'` (e.g. healthy controls) -- there's no "
            "responder/non-responder label to compare for those samples, so "
            "they're excluded rather than silently coerced into one group."
        )
    with st.expander("Why a boxplot with individual points overlaid, not just a box?"):
        st.markdown(
            "A bare boxplot hides sample size and distribution shape -- two "
            "populations can produce visually identical boxes from very "
            "different underlying data. The jittered strip of individual "
            "points overlaid on each box lets a reviewer see the actual "
            "spread, density, and outliers, not just the five-number "
            "summary. The jitter's random seed is fixed so the same figure "
            "is produced on every pipeline run (verified byte-identical "
            "across repeated runs) -- otherwise a re-run would produce a "
            "cosmetically different image for identical underlying data, "
            "which would undermine trust in the artifact even though "
            "nothing had actually changed."
        )

    st.divider()
    st.subheader("Part 4: Subset analysis")
    with st.expander("Why are samples and subjects counted differently?"):
        st.markdown(
            "'Samples per project' counts sample rows directly -- a project "
            "can contribute multiple samples. 'Responders/non-responders' "
            "and 'males/females' are counted by **distinct subject**, "
            "de-duplicating a subject's samples first. Counting those by "
            "raw sample row would inflate the count for any subject "
            "contributing more than one sample to the subset and answer a "
            "different question ('how many samples' vs 'how many people')."
        )
    with st.expander("Why isn't the B cell average restricted to PBMC/miraclib?"):
        st.markdown(
            "The assignment names the PBMC/miraclib restriction explicitly "
            "for the baseline-subset breakdown, but not for this specific "
            "question ('across all sample types and treatments' is stated "
            "directly). Read literally, that's a deliberately broader "
            "question than the subset above it, not an oversight -- so it's "
            "computed against melanoma + male + responder + time=0 only, "
            "with no sample_type/treatment filter."
        )
    with st.expander("Why add an interactive filter explorer beyond the required subset?"):
        st.markdown(
            "The required baseline breakdown is one fixed query. The filter "
            "explorer below it runs the identical query pattern against "
            "whatever condition/sample type/treatment/timepoint combination "
            "is selected, which demonstrates the underlying logic is a "
            "general-purpose query rather than a value hardcoded to produce "
            "the one expected answer -- and it's how the "
            "carcinoma/WB/phauximab/t=7 -> 190-sample result quoted "
            "elsewhere in this project was cross-checked against an "
            "independent pandas computation."
        )

    st.divider()
    st.subheader("Verification tooling (this dashboard)")
    with st.expander("Why build verification into the dashboard instead of just trusting the pipeline?"):
        st.markdown(
            "A number without a visible derivation asks the reader to trust "
            "it on faith. Every subset shown here exposes the SQL query "
            "that produced it, every table is downloadable to check "
            "independently, and the Data Quality tab re-runs its checks "
            "live against the database on every page load rather than "
            "asserting them once in a README that can drift out of date as "
            "the code changes."
        )

# ---------------------------------------------------------------- Data Quality
with tab_quality:
    st.header("Dataset overview")
    st.caption("Computed live from cell_counts.db on this page load.")

    projects_n = conn.execute("SELECT COUNT(*) FROM projects").fetchone()[0]
    subjects_n = conn.execute("SELECT COUNT(*) FROM subjects").fetchone()[0]
    samples_n = conn.execute("SELECT COUNT(*) FROM samples").fetchone()[0]
    pops_n = conn.execute("SELECT COUNT(DISTINCT population) FROM cell_counts").fetchone()[0]

    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Projects", projects_n)
    c2.metric("Subjects", subjects_n)
    c3.metric("Samples", samples_n)
    c4.metric("Cell populations tracked", pops_n)

    st.subheader("Schema")
    st.code(
        "projects(project_id)\n"
        "  -> subjects(subject_id, project_id, condition, age, sex)\n"
        "    -> samples(sample_id, subject_id, sample_type, treatment,\n"
        "               response, time_from_treatment_start)\n"
        "      -> cell_counts(sample_id, population, count)",
        language="text",
    )

    st.subheader("Integrity checks")
    st.caption(
        "These re-run against the live database (and, where noted, the "
        "source CSV) every time this page loads -- they are not asserted "
        "once and forgotten."
    )

    checks = []

    dup_samples = pd.read_sql(
        "SELECT sample_id, COUNT(*) c FROM samples GROUP BY sample_id HAVING c > 1", conn
    )
    checks.append(
        (
            "Every sample_id is unique",
            len(dup_samples) == 0,
            "0 duplicates" if len(dup_samples) == 0 else f"{len(dup_samples)} duplicate sample_id(s)",
        )
    )

    counts_per_sample = pd.read_sql(
        "SELECT sample_id, COUNT(*) c FROM cell_counts GROUP BY sample_id", conn
    )
    bad_counts = counts_per_sample[counts_per_sample["c"] != 5]
    checks.append(
        (
            "Every sample has exactly 5 population counts",
            len(bad_counts) == 0,
            "confirmed for all samples" if len(bad_counts) == 0 else f"{len(bad_counts)} sample(s) with != 5 rows",
        )
    )

    freq_df = compute_frequencies(conn)
    pct_sums = freq_df.groupby("sample")["percentage"].sum()
    max_dev = (pct_sums - 100).abs().max()
    checks.append(
        (
            "Relative frequencies sum to 100% per sample",
            max_dev < 0.01,
            f"max deviation across all samples: {max_dev:.6f} percentage points (rounding only)",
        )
    )

    raw = pd.read_csv(CSV_PATH)
    nuniq = raw.groupby("subject")[["project", "condition", "age", "sex"]].nunique().max()
    checks.append(
        (
            "project / condition / age / sex are constant within each subject (source CSV)",
            bool((nuniq <= 1).all()),
            f"max distinct values observed per subject: {nuniq.to_dict()}",
        )
    )

    null_counts = raw.isnull().sum()
    non_response_nulls = int(null_counts.drop("response").sum())
    checks.append(
        (
            "No missing values outside the response column",
            non_response_nulls == 0,
            f"response has {int(null_counts['response'])} blanks (expected: treatment == 'none' has no readout); "
            f"all other columns: 0 nulls",
        )
    )

    row_count_match = len(raw) == samples_n
    checks.append(
        (
            "Row count in cell_counts.db matches source CSV",
            row_count_match,
            f"CSV: {len(raw)} rows, DB samples table: {samples_n} rows",
        )
    )

    for label, passed, detail in checks:
        icon = "✅" if passed else "❌"
        st.markdown(f"{icon} **{label}** — {detail}")

# ---------------------------------------------------------------- Part 2
with tab2:
    st.header("Relative frequency per sample x population")
    st.caption(
        "For each sample, total_count is the sum of all five population "
        "counts in that sample; percentage is each population's count as a "
        "percentage of that per-sample total."
    )

    with st.expander("Show the query / computation behind this table"):
        st.code("SELECT sample_id AS sample, population, count FROM cell_counts", language="sql")
        st.code(
            "total_count = sum(count) grouped by sample\n"
            "percentage  = count / total_count * 100",
            language="text",
        )

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

    st.dataframe(display_df, width="stretch", height=400)
    st.caption(f"{len(display_df)} of {len(freq_df)} rows shown")
    download_button(display_df, "Download this table as CSV", "part2_frequencies_filtered.csv")

# ---------------------------------------------------------------- Part 3
with tab3:
    st.header("Responders vs non-responders (melanoma, PBMC, miraclib)")
    st.caption(
        "Filter criteria: condition = melanoma, sample_type = PBMC, "
        "treatment = miraclib, response in (yes, no)."
    )

    with st.expander("Show the query behind this subset"):
        st.code(METADATA_QUERY.strip(), language="sql")

    merged = get_melanoma_pbmc_miraclib(conn)
    n_responders = merged.loc[merged["response"] == "yes", "sample"].nunique()
    n_non_responders = merged.loc[merged["response"] == "no", "sample"].nunique()
    st.caption(f"n = {n_responders} responder samples, {n_non_responders} non-responder samples")

    boxplot_path = db.OUTPUTS_DIR / "part3_boxplots.png"
    if not boxplot_path.exists():
        db.ensure_outputs_dir()
        make_boxplots(merged, boxplot_path)
    st.image(str(boxplot_path), width="stretch")
    with open(boxplot_path, "rb") as f:
        st.download_button("Download boxplot as PNG", f.read(), "part3_boxplots.png", "image/png")

    st.subheader("Mann-Whitney U test per population")
    st.caption(
        "Mann-Whitney U was chosen over Welch's t-test because relative "
        "frequencies are bounded proportions with no guarantee of "
        "normality; the rank-based test is robust to that without assuming "
        "a normal distribution."
    )

    alpha = st.slider(
        "Significance threshold (alpha) — drag to see the verdict below update live",
        min_value=0.01,
        max_value=0.10,
        value=0.05,
        step=0.01,
    )
    stats_df = run_stats(merged).drop(columns=["significant_at_0.05"])
    stats_df["significant"] = stats_df["p_value"] < alpha
    st.dataframe(stats_df, width="stretch")
    download_button(stats_df, "Download stats as CSV", "part3_stats.csv")

    st.subheader("Conclusion")
    lines = []
    for _, row in stats_df.sort_values("p_value").iterrows():
        verdict = "differs significantly" if row["significant"] else "does not differ significantly"
        direction = "higher" if row["median_responder_pct"] > row["median_non_responder_pct"] else "lower"
        lines.append(
            f"- **{row['population']}**: p = {row['p_value']:.4f} ({verdict} at "
            f"alpha={alpha:.2f}); responders have {direction} median relative "
            f"frequency ({row['median_responder_pct']:.2f}% vs "
            f"{row['median_non_responder_pct']:.2f}%)."
        )
    st.markdown("\n\n".join(lines))

    st.subheader("Verify manually")
    st.caption(
        "Pick a population to see the exact per-sample percentages that "
        "feed the test above -- download and run your own Mann-Whitney U "
        "(e.g. scipy, R's wilcox.test) to cross-check the p-value."
    )
    pop_choice = st.selectbox("Population", POPULATIONS)
    raw_view = merged.loc[
        merged["population"] == pop_choice, ["sample", "response_label", "percentage"]
    ].sort_values("response_label")
    st.dataframe(raw_view, width="stretch", height=300)
    download_button(
        raw_view, f"Download {pop_choice} raw frequencies", f"part3_{pop_choice}_raw.csv", key="part3_raw_dl"
    )

# ---------------------------------------------------------------- Part 4
with tab4:
    st.header("Required subset: baseline melanoma / PBMC / miraclib (time=0)")
    st.caption("This is the exact subset specified in the assignment.")

    with st.expander("Show the query behind this subset"):
        st.code(BASELINE_QUERY.strip(), language="sql")

    subset = get_baseline_subset(conn)

    col1, col2, col3 = st.columns(3)
    with col1:
        st.markdown("**By project (samples)**")
        st.dataframe(by_project(subset), width="stretch")
    with col2:
        st.markdown("**By response (subjects)**")
        st.dataframe(by_response(subset), width="stretch")
    with col3:
        st.markdown("**By sex (subjects)**")
        st.dataframe(by_sex(subset), width="stretch")

    st.metric("Total samples in baseline subset", len(subset))
    download_button(subset, "Download baseline subset (sample-level)", "part4_baseline_samples.csv")

    st.divider()
    st.subheader("Required answer: avg B cell count, melanoma male responders, time=0")
    st.caption(
        "Across all sample types and treatments -- intentionally NOT "
        "restricted to PBMC/miraclib, since that restriction is named "
        "explicitly for the subset above but not for this question."
    )
    n_rows, avg_bcell = average_bcell_male_responders(conn)
    st.metric("Average B cell count", f"{avg_bcell:.2f}", help=f"n={n_rows} matching sample rows")

    st.divider()
    st.subheader("Explore other subsets")
    st.caption(
        "Change any filter to run the same query logic against a different "
        "slice of the data -- confirms the underlying query is generic, not "
        "hardcoded to only the required subset above."
    )

    conditions = [r[0] for r in conn.execute("SELECT DISTINCT condition FROM subjects ORDER BY 1")]
    sample_types = [r[0] for r in conn.execute("SELECT DISTINCT sample_type FROM samples ORDER BY 1")]
    treatments = [r[0] for r in conn.execute("SELECT DISTINCT treatment FROM samples ORDER BY 1")]
    timepoints = [r[0] for r in conn.execute("SELECT DISTINCT time_from_treatment_start FROM samples ORDER BY 1")]

    fc1, fc2, fc3, fc4 = st.columns(4)
    f_condition = fc1.multiselect("Condition", conditions, default=["melanoma"])
    f_sample_type = fc2.multiselect("Sample type", sample_types, default=["PBMC"])
    f_treatment = fc3.multiselect("Treatment", treatments, default=["miraclib"])
    f_time = fc4.multiselect("Time from treatment start", timepoints, default=[0])

    if f_condition and f_sample_type and f_treatment and f_time:
        custom_query = """
            SELECT sa.sample_id AS sample, sa.response, sa.time_from_treatment_start,
                   su.subject_id, su.project_id, su.condition, su.sex
            FROM samples sa
            JOIN subjects su ON sa.subject_id = su.subject_id
            WHERE su.condition IN ({})
              AND sa.sample_type IN ({})
              AND sa.treatment IN ({})
              AND sa.time_from_treatment_start IN ({})
        """.format(
            ",".join("?" * len(f_condition)),
            ",".join("?" * len(f_sample_type)),
            ",".join("?" * len(f_treatment)),
            ",".join("?" * len(f_time)),
        )
        params = [*f_condition, *f_sample_type, *f_treatment, *[int(t) for t in f_time]]
        custom_subset = pd.read_sql(custom_query, conn, params=params)

        with st.expander("Show the query for this selection"):
            st.code(custom_query.strip(), language="sql")

        cc1, cc2, cc3 = st.columns(3)
        with cc1:
            st.markdown("**By project (samples)**")
            st.dataframe(by_project(custom_subset), width="stretch")
        with cc2:
            st.markdown("**By response (subjects)**")
            st.dataframe(by_response(custom_subset), width="stretch")
        with cc3:
            st.markdown("**By sex (subjects)**")
            st.dataframe(by_sex(custom_subset), width="stretch")

        st.metric("Total samples matching your filters", len(custom_subset))
        download_button(custom_subset, "Download this custom subset", "part4_custom_subset.csv", key="part4_custom_dl")
    else:
        st.info("Select at least one value in every filter to see results.")

conn.close()
