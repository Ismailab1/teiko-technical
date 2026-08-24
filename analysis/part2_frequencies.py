"""
Part 2: relative frequency of each cell population within each sample.

For every sample, total_count is the sum of all five population counts in
that sample, and percentage is each population's count as a percentage of
that total. Writes outputs/part2_frequencies.csv with columns:
sample, total_count, population, count, percentage.
"""

import pandas as pd

try:
    from . import db
except ImportError:
    import db


def compute_frequencies(conn):
    cell_counts = pd.read_sql(
        "SELECT sample_id AS sample, population, count FROM cell_counts", conn
    )
    totals = cell_counts.groupby("sample")["count"].sum().rename("total_count")
    result = cell_counts.merge(totals, on="sample")
    result["percentage"] = (result["count"] / result["total_count"] * 100).round(4)
    result = result[["sample", "total_count", "population", "count", "percentage"]]
    return result.sort_values(["sample", "population"]).reset_index(drop=True)


def main():
    conn = db.get_connection()
    result = compute_frequencies(conn)
    conn.close()

    out_dir = db.ensure_outputs_dir()
    out_path = out_dir / "part2_frequencies.csv"
    result.to_csv(out_path, index=False)
    print(f"Wrote {len(result)} rows to {out_path}")
    print(result.head())


if __name__ == "__main__":
    main()
