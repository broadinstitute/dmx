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

    # Reuse the submission's exact engines - the prediction is computed identically, never re-derived.
    from nightshift_submission import combo_viability, panel_predictions  # noqa: E402

    # The wet-lab oracle is the held-out answer key. It lives in the nightshift repo, NOT this catalog,
    # and a contestant never sees it - which is why this scoring notebook is organizer-side only.
    NIGHTSHIFT = NOTEBOOK_DIR.parent.parent / "nightshift-dev"
    ORACLE = {
        24: NIGHTSHIFT / "data" / "processed" / "nb03_oracle" / "ranking_24h_single.csv",
        48: NIGHTSHIFT / "data" / "processed" / "nb03_oracle" / "ranking_48h_single.csv",
    }
    ORACLE_COMBO = NIGHTSHIFT / "data" / "processed" / "nb03_oracle" / "ranking_48h_combo.csv"
    OUT_DIR = NOTEBOOK_DIR.parent / "data" / "processed" / "nightshift_eval"

    # The rank-scorable (quantitative) tasks. Single agents (1.x) + combinations (2.x); 3.x is a single
    # nomination with no matching measured triple and 4.1 is free text, so neither is rank-scorable here.
    SINGLE_TASKS = [
        ("1.1", "A375", 24),
        ("1.2", "LOXIMVI", 24),
        ("1.3", "A375", 48),
        ("1.4", "LOXIMVI", 48),
        ("1.5", None, None),
    ]
    COMBO_TASKS = [("2.1", "A375"), ("2.2", "LOXIMVI"), ("2.3", None)]


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    # How lookup-able is each nightshift quantitative task? (organizer-side eval)

    This scores the [PRISM-grounded submission](nightshift_submission.py) against the held-out wet-lab
    ground truth, to estimate how much of each *rank-scorable* task is recoverable from public data
    alone - the "challenge" the benchmark poses. It covers the single-agent rankings (1.1-1.5) and the
    combination rankings (2.1-2.3). The 3-drug nominations (3.x) and the strategy essay (4.1) have no
    rank-scorable ground truth, so they are out of scope here.

    **This is the grader's view, not a submission**, and it needs the held-out ground truth to run. To
    keep that answer key private, this notebook surfaces only **aggregate rank-agreement per task**
    (Spearman / Kendall) - never the per-condition measured values. The predictions come from the
    submission's own `panel_predictions()` / `combo_viability()`, so what is scored is exactly what was
    submitted.
    """)
    return


@app.function
def rank_score(task_id: str, family: str, scope: str, pred: list, truth: list) -> dict:
    """Spearman + Kendall tau-b of predicted vs measured viability for one task (None if n < 3)."""
    ok = len(pred) >= 3
    return {
        "task": task_id,
        "family": family,
        "scope": scope,
        "n": len(pred),
        "spearman": round(spearmanr(pred, truth).statistic, 3) if ok else None,
        "kendall_tau_b": round(kendalltau(pred, truth).statistic, 3) if ok else None,
    }


@app.cell
def _():
    # Predictions come straight from the submission's engines (no re-derivation).
    predictions = panel_predictions()
    pred_via = {
        (r["condition"], r["cell_line"], r["timepoint_h"]): r["pred_viability_pct"]
        for r in predictions.iter_rows(named=True)
    }

    # Single agents: join predicted to the per-line/timepoint oracle.
    truth = pl.concat(
        [
            pl.read_csv(ORACLE[tp])
            .select("cell_line", "condition", "viability_pct")
            .rename({"viability_pct": "oracle"})
            .with_columns(timepoint_h=pl.lit(tp))
            for tp in (24, 48)
        ]
    )
    single = predictions.join(truth, on=["condition", "cell_line", "timepoint_h"], how="inner")
    single_rows = []
    for task_id, line, tp in SINGLE_TASKS:
        sub = (
            single if task_id == "1.5" else single.filter((pl.col("cell_line") == line) & (pl.col("timepoint_h") == tp))
        )
        scope = "pooled (2 lines x 2 timepoints)" if task_id == "1.5" else f"{line} {tp}h"
        single_rows.append(
            rank_score(task_id, "single-agent", scope, sub["pred_viability_pct"].to_list(), sub["oracle"].to_list())
        )

    # Combinations (48h): predict each pair's viability via the submission's Bliss-independence model.
    combo = (
        pl.read_csv(ORACLE_COMBO).select("cell_line", "condition", "viability_pct").rename({"viability_pct": "oracle"})
    )
    combo = pl.DataFrame(
        [
            {**r, "pred": combo_viability(r["condition"], r["cell_line"], 48, pred_via)}
            for r in combo.iter_rows(named=True)
        ]
    )
    combo_rows = []
    for task_id, line in COMBO_TASKS:
        sub = combo if task_id == "2.3" else combo.filter(pl.col("cell_line") == line)
        scope = "pooled (2 lines)" if task_id == "2.3" else f"{line} 48h"
        combo_rows.append(rank_score(task_id, "combination", scope, sub["pred"].to_list(), sub["oracle"].to_list()))

    scores = pl.DataFrame(single_rows + combo_rows)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    scores.write_csv(OUT_DIR / "scores.csv")
    mo.vstack(
        [
            mo.md(
                "## The challenge estimate (all quantitative tasks)\n\n"
                "Aggregate rank agreement of the submission vs the held-out ground truth, per task. Higher = more of "
                "that task is recoverable from public data; the gap to 1.0 (or a negative value) is the residual "
                "difficulty the wet lab adds. (Per-condition measured values are intentionally not shown.)"
            ),
            mo.ui.table(scores, page_size=10),
        ]
    )
    return (scores,)


@app.cell(hide_code=True)
def _(scores):
    s_48 = scores.filter(pl.col("task").is_in(["1.3", "1.4"]))["spearman"].mean()
    c_pooled = scores.filter(pl.col("task") == "2.3")["spearman"][0]
    mo.md(
        f"""
        ## Reading the estimate

        - **Single agents are moderately lookup-able** (48h Spearman ~{s_48:.2f}); 24h is hardest because PRISM has one
          timepoint, so the submission cannot reorder 24h vs 48h while the oracle does.
        - **Combinations are where the public-data baseline breaks** (pooled 2.3 Spearman ~{c_pooled}). The submission
          predicts each pair under Bliss INDEPENDENCE - the product of the single-agent viabilities - which by
          construction cannot see synergy. A pathway-synergistic pair is exactly what a single-agent model ranks
          *weakest*, so the combo ranking comes out near-zero or inverted. This is the honest headline: synergy is not
          derivable from single-agent public data, which is precisely what makes the combination tasks (and the live
          nominate-a-combination task) carry the benchmark.
        - Combos have only n = 4 per line (8 pooled), so the combination rho is noisy - read it as direction, not a
          precise estimate.

        ## To extend

        - Score an AUC-ranked single-agent variant and compare; reproduce nb04's exact tie-band scorer instead of this
          continuous-viability Spearman proxy.
        - For combinations, test whether a measured combo-screen prior (GDSC/NCI-ALMANAC-style) recovers the synergy
          that Bliss independence misses.
        """
    )
    return


if __name__ == "__main__":
    app.run()
