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

    import altair as alt
    import marimo as mo
    import polars as pl
    from scipy.stats import kendalltau, spearmanr

    NOTEBOOK_DIR = Path(__file__).resolve().parent
    if str(NOTEBOOK_DIR) not in sys.path:
        sys.path.insert(0, str(NOTEBOOK_DIR))

    # Reuse the submission's exact PRISM engine - the prediction is computed identically, never re-derived.
    from nightshift_single_agent_submission import panel_predictions  # noqa: E402

    # The wet-lab oracle is the held-out answer key. It lives in the nightshift repo, NOT this catalog,
    # and a contestant never sees it - which is why this scoring notebook is organizer-side only and is
    # deliberately not published to molab (it would need this private file to render).
    NIGHTSHIFT = NOTEBOOK_DIR.parent.parent / "nightshift-dev"
    ORACLE = {
        24: NIGHTSHIFT / "data" / "processed" / "nb03_oracle" / "ranking_24h_single.csv",
        48: NIGHTSHIFT / "data" / "processed" / "nb03_oracle" / "ranking_48h_single.csv",
    }
    OUT_DIR = NOTEBOOK_DIR.parent / "data" / "processed" / "nightshift_single_agent_eval"

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

    This scores the [PRISM-grounded submission](nightshift_single_agent_submission.py) against the
    **private wet-lab oracle**, to estimate how much of each ranking task is recoverable from public
    data alone - the "challenge" the benchmark poses.

    **This is the grader's view, not a submission.** It reads the held-out oracle, which a contestant
    never has; that is exactly why the submission notebook does *not* do this. It is kept separate and
    is not a molab artifact (it needs the private oracle to run).

    The predictions come from the submission's `panel_predictions()` - the same PRISM engine, so the
    thing being scored here is identical to what was submitted.
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
    # The submission's predictions (fresh PRISM pull) joined to the oracle ground truth.
    predictions = panel_predictions()
    oracle = pl.concat(
        [
            pl.read_csv(ORACLE[tp])
            .select("cell_line", "condition", "viability_pct")
            .rename({"viability_pct": "oracle_viability_pct"})
            .with_columns(timepoint_h=pl.lit(tp))
            for tp in (24, 48)
        ]
    )
    joined = predictions.join(oracle, on=["condition", "cell_line", "timepoint_h"], how="inner")
    mo.vstack(
        [
            mo.md("## Predicted (PRISM) next to oracle (wet lab)\n\nOne row per condition x line x timepoint."),
            mo.ui.table(joined.sort(["cell_line", "timepoint_h", "oracle_viability_pct"]), page_size=12),
        ]
    )
    return (joined,)


@app.cell
def _(joined):
    scores = pl.DataFrame([task_score(joined, *task) for task in SINGLE_TASKS])
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    scores.write_csv(OUT_DIR / "scores.csv")
    joined.write_csv(OUT_DIR / "predicted_vs_oracle.csv")
    mo.vstack(
        [
            mo.md(
                "## The challenge estimate\n\n"
                "Rank agreement of the PRISM submission vs the oracle, per task. Higher = more of that task "
                "is recoverable from public data; the gap to 1.0 is the residual difficulty the wet lab adds."
            ),
            mo.ui.table(scores, page_size=8),
        ]
    )
    return (scores,)


@app.cell
def _(joined):
    chart = (
        alt.Chart(joined.with_columns(label=pl.format("{} {}h", pl.col("cell_line"), pl.col("timepoint_h"))))
        .mark_circle(size=90, opacity=0.8)
        .encode(
            x=alt.X("pred_viability_pct:Q", title="predicted % viability (PRISM)"),
            y=alt.Y("oracle_viability_pct:Q", title="oracle % viability (wet lab)"),
            color=alt.Color("label:N", title="line / timepoint"),
            tooltip=["condition", "cell_line", "timepoint_h", "pred_viability_pct", "oracle_viability_pct"],
        )
        .properties(width=440, height=380)
    )
    mo.vstack(
        [
            mo.md(
                "## Predicted vs measured\n\n"
                "A positive trend means PRISM and the wet lab agree on the ordering; the 24h points (predicted "
                ">55%) sit apart from the 48h points along x by the kinetic factor. Outliers are the genuine misses."
            ),
            mo.ui.altair_chart(chart),
        ]
    )
    return


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
        - **The standing misses** are Panobinostat (its dose-matched PRISM point understates pan-HDAC potency)
          and, generally, LOXIMVI (the worse-predicted line). These bound how far a public-data baseline can go.

        ## To extend

        - Score the AUC-ranked submission variant and compare per-task - does dose-collapsed AUC beat dose-matched?
        - Reproduce the benchmark's exact tie-band scorer (`nb04_scoring_engine` in the nightshift repo) instead
          of this continuous-viability Spearman proxy, and check the ranking of methods is unchanged.
        """
    )
    return


if __name__ == "__main__":
    app.run()
