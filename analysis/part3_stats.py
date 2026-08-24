"""
Part 3: compare relative frequencies between responders and non-responders,
restricted to melanoma PBMC samples treated with miraclib.

For each of the five populations, runs the Mann-Whitney U test (chosen over
Welch's t-test because sample-level relative frequencies are proportions
bounded in [0, 100] with no guarantee of normality, and Mann-Whitney is
robust to that without requiring a normality assumption). Writes:
  outputs/part3_boxplots.png
  outputs/part3_stats.csv
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu

try:
    from . import db
    from .part2_frequencies import compute_frequencies
except ImportError:
    import db
    from part2_frequencies import compute_frequencies

POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]
ALPHA = 0.05


def get_melanoma_pbmc_miraclib(conn):
    freqs = compute_frequencies(conn)

    meta = pd.read_sql(
        """
        SELECT sa.sample_id AS sample, sa.response, sa.sample_type, sa.treatment,
               su.condition
        FROM samples sa
        JOIN subjects su ON sa.subject_id = su.subject_id
        WHERE su.condition = 'melanoma'
          AND sa.sample_type = 'PBMC'
          AND sa.treatment = 'miraclib'
          AND sa.response IN ('yes', 'no')
        """,
        conn,
    )

    merged = freqs.merge(meta, on="sample")
    merged["response_label"] = merged["response"].map({"yes": "responder", "no": "non-responder"})
    return merged


def run_stats(merged):
    rows = []
    for pop in POPULATIONS:
        sub = merged[merged["population"] == pop]
        responders = sub[sub["response"] == "yes"]["percentage"]
        non_responders = sub[sub["response"] == "no"]["percentage"]
        stat, p = mannwhitneyu(responders, non_responders, alternative="two-sided")
        rows.append(
            {
                "population": pop,
                "n_responders": len(responders),
                "n_non_responders": len(non_responders),
                "median_responder_pct": round(responders.median(), 4),
                "median_non_responder_pct": round(non_responders.median(), 4),
                "u_statistic": stat,
                "p_value": round(p, 6),
                "significant_at_0.05": p < ALPHA,
            }
        )
    return pd.DataFrame(rows).sort_values("p_value").reset_index(drop=True)


def make_boxplots(merged, out_path):
    import numpy as np

    np.random.seed(0)  # deterministic jitter in the stripplot overlay, for reproducible output
    fig, axes = plt.subplots(1, len(POPULATIONS), figsize=(20, 5), sharey=False)
    for ax, pop in zip(axes, POPULATIONS):
        sub = merged[merged["population"] == pop]
        sns.boxplot(
            data=sub,
            x="response_label",
            y="percentage",
            hue="response_label",
            order=["non-responder", "responder"],
            legend=False,
            ax=ax,
        )
        sns.stripplot(
            data=sub,
            x="response_label",
            y="percentage",
            order=["non-responder", "responder"],
            color="black",
            alpha=0.4,
            size=3,
            ax=ax,
        )
        ax.set_title(pop)
        ax.set_xlabel("")
        ax.set_ylabel("relative frequency (%)" if pop == POPULATIONS[0] else "")
    fig.suptitle("Relative frequency by population: responders vs non-responders\n(melanoma, PBMC, miraclib)")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def summarize(stats_df):
    lines = []
    for _, row in stats_df.iterrows():
        verdict = "differs significantly" if row["significant_at_0.05"] else "does not differ significantly"
        direction = "higher" if row["median_responder_pct"] > row["median_non_responder_pct"] else "lower"
        lines.append(
            f"- {row['population']}: p = {row['p_value']:.4f} ({verdict} at alpha=0.05); "
            f"responders have {direction} median relative frequency "
            f"({row['median_responder_pct']:.2f}% vs {row['median_non_responder_pct']:.2f}%)."
        )
    return "\n".join(lines)


def main():
    conn = db.get_connection()
    merged = get_melanoma_pbmc_miraclib(conn)
    conn.close()

    stats_df = run_stats(merged)

    out_dir = db.ensure_outputs_dir()
    stats_path = out_dir / "part3_stats.csv"
    plot_path = out_dir / "part3_boxplots.png"

    stats_df.to_csv(stats_path, index=False)
    make_boxplots(merged, plot_path)

    summary = summarize(stats_df)
    summary_path = out_dir / "part3_summary.txt"
    summary_path.write_text(summary, encoding="utf-8")

    print(f"Wrote {stats_path}")
    print(f"Wrote {plot_path}")
    print(f"Wrote {summary_path}")
    print()
    print(summary)


if __name__ == "__main__":
    main()
