# /// script
# requires-python = ">=3.12,<3.14"
# dependencies = [
#     "altair==5.5.0",
#     "marimo",
#     "polars==1.40.1",
#     "requests==2.32.5",
#     "scipy==1.16.0",
# ]
# ///

import marimo

__generated_with = "0.23.10"
app = marimo.App(width="medium")

with app.setup:
    import sys
    from pathlib import Path

    import marimo as mo
    import polars as pl
    from scipy.stats import kendalltau, spearmanr

    NOTEBOOK_DIR = Path(__file__).resolve().parent
    if str(NOTEBOOK_DIR) not in sys.path:
        sys.path.insert(0, str(NOTEBOOK_DIR))

    # Reuse the submission's exact PRISM engine - the prediction is computed identically, never re-derived.
    from nightshift_submission import panel_predictions  # noqa: E402

    # The wet-lab oracle is the held-out answer key. It lives in the nightshift repo, NOT this catalog,
    # and a contestant never sees it - which is why this scoring notebook is organizer-side only and is
    # deliberately not published to molab (it would need this private file to render).
    NIGHTSHIFT = NOTEBOOK_DIR.parent.parent / "nightshift-dev"
    ORACLE = {
        24: NIGHTSHIFT / "data" / "processed" / "nb03_oracle" / "ranking_24h_single.csv",
        48: NIGHTSHIFT / "data" / "processed" / "nb03_oracle" / "ranking_48h_single.csv",
    }
    OUT_DIR = NOTEBOOK_DIR.parent / "data" / "processed" / "nightshift_eval"

    # (task_id, cell_line, timepoint_h); 1.5 is pooled across both lines and both timepoints.
    SINGLE_TASKS = [
        ("1.1", "A375", 24),
        ("1.2", "LOXIMVI", 24),
        ("1.3", "A375", 48),
        ("1.4", "LOXIMVI", 48),
        ("1.5", None, None),
    ]


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    # How lookup-able is each nightshift single-agent task? (organizer-side eval)

    This scores the [PRISM-grounded submission](nightshift_submission.py) against the
    held-out wet-lab ground truth, to estimate how much of each ranking task is recoverable from
    public data alone - the "challenge" the benchmark poses.

    **This is the grader's view, not a submission**, and it needs the held-out ground truth to run.
    To keep that answer key private, this notebook surfaces only **aggregate rank-agreement per task**
    (Spearman / Kendall) - never the per-drug measured values. The predictions come from the
    submission's `panel_predictions()`, so what is scored is exactly what was submitted.
    """)
    return


@app.function
def task_score(joined: pl.DataFrame, task_id: str, line: str | None, timepoint: int | None) -> dict:
    """Spearman + Kendall tau-b of predicted vs oracle viability for one task's conditions."""
    sub = (
        joined
        if task_id == "1.5"
        else joined.filter((pl.col("cell_line") == line) & (pl.col("timepoint_h") == timepoint))
    )
    pred, truth = sub["pred_viability_pct"].to_list(), sub["oracle_viability_pct"].to_list()
    return {
        "task": task_id,
        "scope": "pooled (2 lines x 2 timepoints)" if task_id == "1.5" else f"{line} {timepoint}h",
        "n": len(pred),
        "spearman": round(spearmanr(pred, truth).statistic, 3),
        "kendall_tau_b": round(kendalltau(pred, truth).statistic, 3),
    }


@app.cell
def _():
    # Score the submission against the held-out ground truth. The join is computed internally but
    # NEVER displayed - only the aggregate per-task agreement is surfaced, so this notebook (and its
    # molab snapshot) never exposes a per-drug measured value, i.e. the answer key.
    predictions = panel_predictions()
    truth = pl.concat(
        [
            pl.read_csv(ORACLE[tp])
            .select("cell_line", "condition", "viability_pct")
            .rename({"viability_pct": "oracle_viability_pct"})
            .with_columns(timepoint_h=pl.lit(tp))
            for tp in (24, 48)
        ]
    )
    joined = predictions.join(truth, on=["condition", "cell_line", "timepoint_h"], how="inner")
    scores = pl.DataFrame([task_score(joined, *task) for task in SINGLE_TASKS])
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    scores.write_csv(OUT_DIR / "scores.csv")
    mo.vstack(
        [
            mo.md(
                "## The challenge estimate\n\n"
                "Aggregate rank agreement of the PRISM submission vs the held-out ground truth, per task. "
                "Higher = more of that task is recoverable from public data; the gap to 1.0 is the residual "
                "difficulty the wet lab adds. (Per-drug measured values are intentionally not shown.)"
            ),
            mo.ui.table(scores, page_size=8),
        ]
    )
    return (scores,)


@app.cell(hide_code=True)
def _(scores):
    a375_48 = scores.filter(pl.col("task") == "1.3")["spearman"][0]
    lox_48 = scores.filter(pl.col("task") == "1.4")["spearman"][0]
    mo.md(
        f"""
        ## Reading the estimate

        - **24h is the hardest (lowest agreement).** PRISM has one timepoint, so the submission's 24h and 48h
          orderings are identical, but the oracle reorders between them - the kinetic component of the task is
          not lookup-able from sensitivity data.
        - **48h is moderately recoverable** (A375 Spearman ~{a375_48}, LOXIMVI ~{lox_48}); the pooled 1.5 scores
          highest because the large 24h-vs-48h level gap, which the kinetic factor captures, dominates 48 points.
        - **The structural limits** are drug-specific kinetics (absent from single-timepoint PRISM) and
          HDAC-inhibitor polypharmacology (under-represented by one dose-matched point), plus the harder-to-
          predict LOXIMVI line. These bound how far a public-data baseline can go.

        ## To extend

        - Score the AUC-ranked submission variant and compare per-task - does dose-collapsed AUC beat dose-matched?
        - Reproduce the benchmark's exact tie-band scorer (`nb04_scoring_engine` in the nightshift repo) instead
          of this continuous-viability Spearman proxy, and check the ranking of methods is unchanged.
        """
    )
    return


if __name__ == "__main__":
    app.run()
