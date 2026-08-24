# Teiko Technical — Immune Cell Count Analysis

Analysis pipeline and dashboard for the `cell-count.csv` clinical trial
dataset: cell population frequencies, a responder vs non-responder
statistical comparison, and a baseline-subset breakdown.

## Setup & Run

Requires Python 3.11+.

```bash
make setup       # pip install -r requirements.txt
make pipeline    # load_data.py + Part 2/3/4 analysis scripts -> outputs/
make dashboard   # streamlit run dashboard/app.py
```

- `make setup` installs dependencies.
- `make pipeline` builds `cell_counts.db` from `cell-count.csv` (drops and
  recreates the DB each run, so it's safe to re-run) and writes every
  Part 2/3/4 result as a CSV/PNG under `outputs/`.
- `make dashboard` starts a local Streamlit server presenting all three
  parts interactively. Run `make pipeline` at least once first so the DB
  exists.

This works unattended in a fresh GitHub Codespace: clone, `make setup`,
`make pipeline`, `make dashboard`.

**Live dashboard:** https://teiko-technical-ktfspx47mungbqpdzj5cxl.streamlit.app/

## Database Schema

```
projects(project_id PK)
subjects(subject_id PK, project_id FK, condition, age, sex)
samples(sample_id PK, subject_id FK, sample_type, treatment, response, time_from_treatment_start)
cell_counts(sample_id FK, population, count)   -- PK (sample_id, population)
```

Indexes on every foreign key and every column commonly filtered on
(`condition`, `treatment`, `response`, `sample_type`, `time_from_treatment_start`,
`population`).

**Why normalized, not one wide table:**

- `projects -> subjects -> samples -> cell_counts` avoids repeating
  subject-level facts (condition, age, sex) across a subject's three
  samples, so they can't drift out of sync.
- `treatment` and `response` live on `samples`, not `subjects`, even
  though in this snapshot they're constant per subject. They're
  trial-arm/readout concepts that could legitimately vary by sample in
  future data (e.g. a subject switching treatment arms, or a response
  readout assessed per timepoint). Modeling them one level down keeps the
  schema correct as the data model matures, at the cost of a little
  redundancy today.
- `cell_counts` is long/tidy (one row per sample x population) instead of
  one column per population. Adding a 6th, 7th, ... cell population later
  needs zero schema migration and zero loader code changes — just more
  rows. It also makes every downstream query (frequencies, filtering by
  population, stats) a plain `GROUP BY`/`JOIN` instead of bespoke
  per-column logic.

**Scaling to hundreds of projects / thousands of samples / new analysis
types:**

1. Move from SQLite to Postgres for concurrent writes and multi-user
   access — the schema above is already vanilla relational SQL and
   ports over unchanged.
2. Add a `populations` reference table once population definitions carry
   more metadata than just a name (e.g. marker panels, parent/child gating
   hierarchy).
3. Precompute a `sample_frequencies` table refreshed on load if the Part 2
   computation becomes a query-time bottleneck at scale.
4. Add a `batches`/`ingestion_runs` table if data starts arriving
   incrementally rather than as one full-file load, so `load_data.py`
   can append instead of drop-and-recreate.

## Code Structure

```
load_data.py               # builds schema, loads CSV into cell_counts.db (idempotent)
analysis/
  db.py                     # shared SQLite connection + outputs/ dir helper
  part2_frequencies.py      # tidy sample x population frequency table
  part3_stats.py            # responder vs non-responder boxplots + Mann-Whitney U
  part4_subset.py           # baseline subset breakdown + B cell average
dashboard/
  app.py                    # Streamlit app, reuses the analysis/ functions directly
outputs/                    # CSV/PNG artifacts written by the pipeline
cell-count.csv              # source data
Makefile
requirements.txt
```

Each `analysis/part*.py` module is runnable standalone
(`python analysis/part2_frequencies.py`) or as part of `make pipeline`,
and writes its results to `outputs/` rather than only printing — so
outputs exist as on-disk artifacts, not just terminal output. The
dashboard imports the same compute functions (`compute_frequencies`,
`run_stats`, `get_baseline_subset`, etc.) instead of duplicating query
logic, so the CLI pipeline and the interactive dashboard can never
disagree.

## Results Summary

**Part 2** — `outputs/part2_frequencies.csv`: 52,500 rows (10,500 samples
x 5 populations), columns `sample, total_count, population, count,
percentage`.

**Part 3** — melanoma / PBMC / miraclib, responders vs non-responders
(Mann-Whitney U test; see justification below), `outputs/part3_boxplots.png`
and `outputs/part3_stats.csv`:

| population | p-value | significant (alpha=0.05) |
|---|---|---|
| cd4_t_cell | 0.0133 | yes |
| b_cell | 0.0557 | no |
| nk_cell | 0.1211 | no |
| monocyte | 0.1631 | no |
| cd8_t_cell | 0.6391 | no |

Only **cd4_t_cell** relative frequency differs significantly between
responders and non-responders (responders trend higher). The other four
populations show no statistically significant difference at this sample
size.

*Why Mann-Whitney U rather than Welch's t-test:* the values being compared
are per-sample relative frequencies (bounded proportions in `[0, 100]`),
not raw counts, and there's no reason to assume they're normally
distributed. Mann-Whitney U is a non-parametric rank-based test that
doesn't require that assumption and is robust to the outliers visible in
the boxplots, at the cost of some statistical power versus a t-test if the
data actually were normal.

**Part 4** — baseline (`time_from_treatment_start=0`) melanoma / PBMC /
miraclib subset:

- Total samples: **656** (prj1: 384, prj3: 272, prj2: 0)
- By response (distinct subjects): 331 responders / 325 non-responders
- By sex (distinct subjects): 344 male / 312 female
- **Average B cell count, melanoma male responders, time=0, across all
  sample types and treatments: 10206.15** (n=485 matching sample rows)

Full breakdowns: `outputs/part4_baseline_by_project.csv`,
`part4_baseline_by_response.csv`, `part4_baseline_by_sex.csv`,
`part4_avg_bcell_answer.csv`.

## Assumptions

- Part 4's final question ("average B cell count for melanoma male
  responders at time=0, across all sample types and treatments") is
  computed without restricting to PBMC/miraclib — that restriction is
  named explicitly in the earlier baseline-subset question but not in
  this one, so it's read as intentionally broader.
- "Responders vs non-responders" in Part 3 excludes samples with a blank
  `response` (i.e. `treatment == "none"`, which has no response readout by
  definition).
- The `.db` file is git-ignored since it's fully derived from
  `cell-count.csv` via `load_data.py`; `outputs/` is committed so the
  CSV/PNG artifacts are visible without running the pipeline.

## Dashboard

`make dashboard` runs `streamlit run dashboard/app.py`, presenting three
tabs: the Part 2 frequency table (filterable by sample/population), the
Part 3 boxplots and significance table, and Part 4's baseline breakdown
with the B cell average.

Live link: https://teiko-technical-ktfspx47mungbqpdzj5cxl.streamlit.app/
