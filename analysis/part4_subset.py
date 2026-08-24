"""
Part 4: baseline (time_from_treatment_start=0) melanoma/PBMC/miraclib subset,
broken down by project, response, and sex; plus the average B cell count for
melanoma male responders at time=0, across all sample types and treatments.

Writes:
  outputs/part4_baseline_by_project.csv
  outputs/part4_baseline_by_response.csv
  outputs/part4_baseline_by_sex.csv
  outputs/part4_avg_bcell_answer.csv
"""

import pandas as pd

try:
    from . import db
except ImportError:
    import db


BASELINE_QUERY = """
SELECT sa.sample_id AS sample, sa.response, sa.time_from_treatment_start,
       su.subject_id, su.project_id, su.condition, su.sex
FROM samples sa
JOIN subjects su ON sa.subject_id = su.subject_id
WHERE su.condition = 'melanoma'
  AND sa.sample_type = 'PBMC'
  AND sa.treatment = 'miraclib'
  AND sa.time_from_treatment_start = 0
"""

BCELL_QUERY = """
SELECT su.subject_id, su.sex, sa.response, cc.count
FROM samples sa
JOIN subjects su ON sa.subject_id = su.subject_id
JOIN cell_counts cc ON cc.sample_id = sa.sample_id
WHERE su.condition = 'melanoma'
  AND su.sex = 'M'
  AND sa.response = 'yes'
  AND sa.time_from_treatment_start = 0
  AND cc.population = 'b_cell'
"""


def get_baseline_subset(conn):
    return pd.read_sql(BASELINE_QUERY, conn)


def by_project(subset):
    return (
        subset.groupby("project_id")
        .size()
        .rename("sample_count")
        .reset_index()
        .rename(columns={"project_id": "project"})
    )


def by_response(subset):
    by_subject = subset.drop_duplicates("subject_id")
    return (
        by_subject.groupby("response", dropna=False)
        .size()
        .rename("subject_count")
        .reset_index()
    )


def by_sex(subset):
    by_subject = subset.drop_duplicates("subject_id")
    return by_subject.groupby("sex").size().rename("subject_count").reset_index()


def average_bcell_male_responders(conn):
    df = pd.read_sql(BCELL_QUERY, conn)
    return len(df), round(df["count"].mean(), 2)


def main():
    conn = db.get_connection()
    subset = get_baseline_subset(conn)
    n_rows, avg_bcell = average_bcell_male_responders(conn)
    conn.close()

    out_dir = db.ensure_outputs_dir()

    proj_df = by_project(subset)
    resp_df = by_response(subset)
    sex_df = by_sex(subset)

    proj_df.to_csv(out_dir / "part4_baseline_by_project.csv", index=False)
    resp_df.to_csv(out_dir / "part4_baseline_by_response.csv", index=False)
    sex_df.to_csv(out_dir / "part4_baseline_by_sex.csv", index=False)

    answer_df = pd.DataFrame(
        [{"matching_rows": n_rows, "avg_b_cell_count": avg_bcell}]
    )
    answer_df.to_csv(out_dir / "part4_avg_bcell_answer.csv", index=False)

    print(f"Baseline subset: {len(subset)} samples")
    print(proj_df.to_string(index=False))
    print(resp_df.to_string(index=False))
    print(sex_df.to_string(index=False))
    print(
        f"\nAverage B cell count, melanoma male responders, time=0, "
        f"all sample types/treatments: {avg_bcell:.2f} (n={n_rows} rows)"
    )


if __name__ == "__main__":
    main()
