"""
Load cell-count.csv into a normalized SQLite database (cell_counts.db).

Schema: projects -> subjects -> samples -> cell_counts (long/tidy format).
See README.md for the full rationale. Running this script twice is safe:
the existing .db file is dropped and rebuilt from scratch each time.
"""

import sqlite3
from pathlib import Path

import pandas as pd

CSV_PATH = Path(__file__).parent / "cell-count.csv"
DB_PATH = Path(__file__).parent / "cell_counts.db"

POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

SCHEMA = """
CREATE TABLE projects (
    project_id TEXT PRIMARY KEY
);

CREATE TABLE subjects (
    subject_id TEXT PRIMARY KEY,
    project_id TEXT NOT NULL REFERENCES projects(project_id),
    condition  TEXT NOT NULL,
    age        INTEGER,
    sex        TEXT
);

CREATE TABLE samples (
    sample_id                 TEXT PRIMARY KEY,
    subject_id                TEXT NOT NULL REFERENCES subjects(subject_id),
    sample_type               TEXT NOT NULL,
    treatment                 TEXT,
    response                  TEXT,
    time_from_treatment_start INTEGER
);

CREATE TABLE cell_counts (
    sample_id  TEXT NOT NULL REFERENCES samples(sample_id),
    population TEXT NOT NULL,
    count      INTEGER NOT NULL,
    PRIMARY KEY (sample_id, population)
);

CREATE INDEX idx_subjects_project      ON subjects(project_id);
CREATE INDEX idx_subjects_condition    ON subjects(condition);
CREATE INDEX idx_samples_subject       ON samples(subject_id);
CREATE INDEX idx_samples_treatment     ON samples(treatment);
CREATE INDEX idx_samples_response      ON samples(response);
CREATE INDEX idx_samples_type          ON samples(sample_type);
CREATE INDEX idx_samples_time          ON samples(time_from_treatment_start);
CREATE INDEX idx_cellcounts_sample     ON cell_counts(sample_id);
CREATE INDEX idx_cellcounts_population ON cell_counts(population);
"""


def main():
    if DB_PATH.exists():
        DB_PATH.unlink()

    df = pd.read_csv(CSV_PATH)

    conn = sqlite3.connect(DB_PATH)
    conn.executescript(SCHEMA)

    projects = df[["project"]].drop_duplicates().rename(columns={"project": "project_id"})
    projects.to_sql("projects", conn, if_exists="append", index=False)

    subjects = (
        df[["subject", "project", "condition", "age", "sex"]]
        .drop_duplicates(subset="subject")
        .rename(columns={"subject": "subject_id", "project": "project_id"})
    )
    subjects.to_sql("subjects", conn, if_exists="append", index=False)

    samples = df[
        ["sample", "subject", "sample_type", "treatment", "response", "time_from_treatment_start"]
    ].rename(columns={"sample": "sample_id", "subject": "subject_id"})
    samples.to_sql("samples", conn, if_exists="append", index=False)

    cell_counts = df.melt(
        id_vars=["sample"],
        value_vars=POPULATIONS,
        var_name="population",
        value_name="count",
    ).rename(columns={"sample": "sample_id"})
    cell_counts.to_sql("cell_counts", conn, if_exists="append", index=False)

    conn.commit()

    counts = {
        "projects": conn.execute("SELECT COUNT(*) FROM projects").fetchone()[0],
        "subjects": conn.execute("SELECT COUNT(*) FROM subjects").fetchone()[0],
        "samples": conn.execute("SELECT COUNT(*) FROM samples").fetchone()[0],
        "cell_counts": conn.execute("SELECT COUNT(*) FROM cell_counts").fetchone()[0],
    }
    conn.close()

    print(f"Loaded {DB_PATH.name}:")
    for table, n in counts.items():
        print(f"  {table}: {n} rows")


if __name__ == "__main__":
    main()
