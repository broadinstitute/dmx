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
    import copy
    import json
    import math
    import re
    import sys
    from pathlib import Path

    import altair as alt
    import marimo as mo
    import polars as pl
    from scipy.stats import kendalltau, spearmanr

    NOTEBOOK_DIR = Path(__file__).resolve().parent
    if str(NOTEBOOK_DIR) not in sys.path:
        sys.path.insert(0, str(NOTEBOOK_DIR))

    from nb02_dataset_discovery import dataset_features  # noqa: E402
    from nb03_gene_dependency_profile import matrix_feature  # noqa: E402

    # The nightshift benchmark (task prompts + the local oracle) lives next to this catalog.
    NIGHTSHIFT = NOTEBOOK_DIR.parent.parent / "nightshift-dev"
    TASKS_DIR = NIGHTSHIFT / "data" / "raw" / "tasks"
    ORACLE = {
        24: NIGHTSHIFT / "data" / "processed" / "nb03_oracle" / "ranking_24h_single.csv",
        48: NIGHTSHIFT / "data" / "processed" / "nb03_oracle" / "ranking_48h_single.csv",
    }

    OUT_DIR = NOTEBOOK_DIR.parent / "data" / "processed" / "nightshift_single_agent_submission"
    SUB_DIR = OUT_DIR / "submission"

    CELL_LINES = {"A375": "ACH-000219", "LOXIMVI": "ACH-000750"}

    # PRISM Repurposing Secondary, dose-resolved: value = log2 fold-change vs DMSO (more negative
    # = more killing). The engine of this submission - measured drug response in the same two lines.
    PRISM_VIAB = "576e1cb6-ac8d-4e29-bf15-0552c8665d72"

    # The five single-agent tasks. (task_id, cell_line, timepoint_h); 1.5 is pooled across both.
    SINGLE_TASKS = [
        ("1.1", "A375", 24),
        ("1.2", "LOXIMVI", 24),
        ("1.3", "A375", 48),
        ("1.4", "LOXIMVI", 48),
        ("1.5", None, None),
    ]

    # Kinetic factor: PRISM is a ~5-day endpoint; the tasks read at 24h / 48h. We treat PRISM's
    # killing as the 48h magnitude (factor 1.0) and attenuate the 24h prediction toward baseline.
    # This is a stated prior, NOT fit to the oracle; it only shifts absolute level and how 24h/48h
    # interleave in the pooled 1.5 ranking - it cannot reorder drugs within one (line, timepoint).
    KINETIC_FACTOR = {24: 0.55, 48: 1.0}

    # The 12 single agents with their fixed nightshift dose and PRISM label (Sapanisertib is
    # screened under its code MLN0128). Doses taken verbatim from the task prompts.
    PANEL = [
        {"drug": "Panobinostat", "prism": "PANOBINOSTAT", "dose_uM": 0.05},
        {"drug": "Trametinib", "prism": "TRAMETINIB", "dose_uM": 0.01},
        {"drug": "Dabrafenib", "prism": "DABRAFENIB", "dose_uM": 0.1},
        {"drug": "Encorafenib", "prism": "ENCORAFENIB", "dose_uM": 0.1},
        {"drug": "Cobimetinib", "prism": "COBIMETINIB", "dose_uM": 0.1},
        {"drug": "Binimetinib", "prism": "BINIMETINIB", "dose_uM": 0.1},
        {"drug": "TAK-733", "prism": "TAK-733", "dose_uM": 0.03},
        {"drug": "Vemurafenib", "prism": "VEMURAFENIB", "dose_uM": 1.0},
        {"drug": "Regorafenib", "prism": "REGORAFENIB", "dose_uM": 5.0},
        {"drug": "Sapanisertib", "prism": "MLN0128", "dose_uM": 0.5},
        {"drug": "Capivasertib", "prism": "CAPIVASERTIB", "dose_uM": 5.0},
        {"drug": "Alpelisib", "prism": "ALPELISIB", "dose_uM": 5.0},
    ]


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    # A DepMap-grounded submission to nightshift single-agent tasks 1.1-1.5

    The nightshift tasks hand an agent a real melanoma viability experiment - cell line, the
    12 compounds, their exact doses, the timepoint - and ask it to **predict** the % viability
    and rank, from "known pharmacology, target, concentration, cell line, timepoint," with **no
    readout data and no requesting raw data**. So a public-data lookup is the sanctioned move.

    This notebook *is* that lookup, made concrete. Tasks 1.1-1.4 are A375/LOXIMVI at 24h/48h;
    1.5 pools all 48 conditions into one ranking. The engine is **PRISM Repurposing Secondary**
    (DepMap), read at the screened concentration nearest each drug's task dose:

    1. PRISM log2 fold-change at the matched dose -> predicted fractional viability `2**log2fc`.
    2. A stated 24h/48h kinetic factor turns that into a predicted % viability per timepoint.
    3. Rank each task's conditions by predicted viability (1 = strongest = lowest viability).

    It writes the five populated `output.json` submissions, then **self-scores** against the
    local oracle - the honest "how lookup-able is each task" readout. (A real contestant has no
    oracle; the self-score is ours, to see how far public data gets.)
    """)
    return


@app.function
def two_line_values(dataset_id: str, feature_label: str) -> dict[str, float | None]:
    """Pull one matrix feature and keep only the two nightshift cell lines, keyed by name."""
    frame = matrix_feature(dataset_id, feature_label, value_name="value")
    by_id = dict(zip(frame["depmap_id"], frame["value"], strict=False)) if frame.height else {}
    return {name: by_id.get(model_id) for name, model_id in CELL_LINES.items()}


@app.function
def nearest_dose_label(viab_features: pl.DataFrame, prism_name: str, target_uM: float) -> tuple[str, float] | None:
    """Pick the dose-resolved PRISM feature nearest (in log10 uM) to the nightshift dose.

    Labels are 'NAME <dose> uM' with per-compound rounding, so match the live feature list rather
    than reconstructing a label. Returns (label, picked_uM); a dose exactly between two grid points
    is a genuine tie (round before comparing) broken to the lower, more conservative dose.
    """
    pattern = re.compile(rf"^{re.escape(prism_name)} ([\d.]+) uM$")
    candidates = [(float(m.group(1)), label) for label in viab_features["label"] if (m := pattern.match(label))]
    if not candidates:
        return None
    candidates.sort(key=lambda d: (round(abs(math.log10(d[0]) - math.log10(target_uM)), 6), d[0]))
    dose_uM, label = candidates[0]
    return label, dose_uM


@app.function
def predicted_viability_pct(prism_log2fc: float, timepoint_h: int) -> float:
    """Predicted % viability at a timepoint from a PRISM log2 fold-change.

    PRISM fractional viability is 2**log2fc; the killed fraction (effect) is scaled by the
    timepoint's kinetic factor and subtracted from baseline. Clipped to [0, 110].
    """
    effect = max(0.0, 1.0 - 2.0**prism_log2fc)  # growth (>baseline) counts as zero killing
    viability = 100.0 * (1.0 - KINETIC_FACTOR[timepoint_h] * effect)
    return max(0.0, min(110.0, viability))


@app.function
def fill_submission(template: dict, pred_via: dict) -> dict:
    """Populate a task's output.json: rank + viability_pct from predicted viability.

    `pred_via` is keyed (condition, cell_line, timepoint_h) -> % viability. Ranks are computed
    among the template's own entries (1 = lowest viability = strongest), so per-line tasks rank
    within 12 and the pooled 1.5 ranks across all 48.
    """
    out = copy.deepcopy(template)
    task_line, task_tp = template.get("cell_line"), template.get("timepoint_h")
    entries = out["rankings"]
    values = [pred_via[(e["condition"], e.get("cell_line", task_line), e.get("timepoint_h", task_tp))] for e in entries]
    order = sorted(range(len(entries)), key=lambda i: values[i])
    for position, i in enumerate(order, start=1):
        entries[i]["rank"] = position
        entries[i]["viability_pct"] = round(values[i], 1)
    return out


@app.cell
def _():
    # Engine: PRISM viability at the dose nearest each drug's task dose, for both lines.
    viab_features = dataset_features(PRISM_VIAB)
    signal_rows = []
    for _c in PANEL:
        _pick = nearest_dose_label(viab_features, _c["prism"], _c["dose_uM"])
        _label, _picked_uM = _pick if _pick else (None, None)
        _vals = two_line_values(PRISM_VIAB, _label) if _label else {_k: None for _k in CELL_LINES}
        for _line, _log2fc in _vals.items():
            signal_rows.append(
                {
                    "drug": _c["drug"],
                    "cell_line": _line,
                    "dose_uM": _c["dose_uM"],
                    "prism_dose_uM": _picked_uM,
                    "prism_log2fc": _log2fc,
                    "prism_viab_pct": round(100.0 * 2.0**_log2fc, 1) if _log2fc is not None else None,
                }
            )
    prism_signal = pl.DataFrame(signal_rows)
    mo.md(
        "## The PRISM signal\n\n"
        "Dose-resolved PRISM log2 fold-change at the curve point nearest each task dose "
        "(`prism_dose_uM`), and the implied fractional viability. All 12/12 compounds resolve for "
        "both lines. This single number per (drug, line) is the only measured input to the prediction."
    )
    mo.ui.table(prism_signal.sort(["cell_line", "prism_log2fc"]), page_size=12)
    return (prism_signal,)


@app.cell
def _(prism_signal):
    # Predicted % viability for every (drug, line, timepoint): apply the kinetic factor.
    pred_rows = []
    for _row in prism_signal.iter_rows(named=True):
        if _row["prism_log2fc"] is None:
            continue
        for _tp in (24, 48):
            pred_rows.append(
                {
                    "condition": _row["drug"],
                    "cell_line": _row["cell_line"],
                    "timepoint_h": _tp,
                    "pred_viability_pct": round(predicted_viability_pct(_row["prism_log2fc"], _tp), 1),
                }
            )
    predictions = pl.DataFrame(pred_rows)
    pred_via = {(r["condition"], r["cell_line"], r["timepoint_h"]): r["pred_viability_pct"] for r in pred_rows}
    mo.md(
        "## Predicted viability (the submission's raw prediction)\n\n"
        "One predicted % viability per condition x cell line x timepoint. The kinetic factor makes "
        "24h less killed than 48h; it cannot reorder drugs within a (line, timepoint) - that order "
        "is fixed by the PRISM dose signal."
    )
    mo.ui.table(predictions.sort(["cell_line", "timepoint_h", "pred_viability_pct"]), page_size=12)
    return pred_via, predictions


@app.function
def build_reasoning(filled: dict) -> str:
    """Generate the reasoning_md trace for a filled task - method, grounding, ranking, caveats.

    This is the `reasoning_md` the Karman `submit_response` tool wants and the task prompts ask
    for; it is generated from the same prediction the output.json carries, so the two never drift.
    """
    entries = sorted(filled["rankings"], key=lambda e: e["rank"])
    pooled = "cell_line" in entries[0]
    header_cols = (
        "| rank | condition | cell line | timepoint | predicted viability % |"
        if pooled
        else "| rank | condition | conc | predicted viability % |"
    )
    sep = "|---|---|---|---|---|" if pooled else "|---|---|---|---|"
    rows = []
    for e in entries:
        if pooled:
            rows.append(
                f"| {e['rank']} | {e['condition']} | {e['cell_line']} | {e['timepoint_h']}h | {e['viability_pct']} |"
            )
        else:
            rows.append(f"| {e['rank']} | {e['condition']} | {e['concentration']} | {e['viability_pct']} |")
    table = "\n".join([header_cols, sep, *rows])
    return f"""# Reasoning - Night Shift task {filled["task"]}

{filled["description"]}

## Method (public-data lookup, no oracle access)

This ranking is grounded entirely in **DepMap PRISM Repurposing Secondary**, a measured
small-molecule viability screen, for the exact cell lines in this task (A375 = ACH-000219,
LOXIMVI = ACH-000750). For each compound I read PRISM's log2 fold-change vs DMSO at the
screened concentration nearest its stated dose, convert to fractional viability (`2**log2fc`),
and apply a fixed 24h/48h kinetic factor (PRISM is a ~5-day endpoint; 24h is attenuated toward
baseline). Conditions are ranked by predicted viability, rank 1 = strongest effect (lowest
viability). All 12 compounds resolve in PRISM for both lines (Sapanisertib via its code MLN0128).

## Predicted ranking

{table}

## Caveats (stated honestly)

- The ranking is driven by measured PRISM sensitivity; the absolute viability % is approximate
  (different assay, timepoint, and a coarse 4x dose grid).
- PRISM has a single timepoint, so this method cannot reorder drugs between 24h and 48h - the
  within-line ordering is identical at both, which understates fast-vs-slow kinetic differences.
- HDAC inhibitor potency (Panobinostat) is under-represented by its dose-matched PRISM point
  relative to its true effect, so its rank is the least certain.
"""


@app.cell
def _(pred_via):
    # Build + write the five populated output.json submissions (+ reasoning) from the templates.
    SUB_DIR.mkdir(parents=True, exist_ok=True)
    filled = {}
    for _task_id, _line, _tp in SINGLE_TASKS:
        _template = json.loads((TASKS_DIR / f"task_{_task_id}" / "output.json").read_text())
        _out = fill_submission(_template, pred_via)
        filled[_task_id] = _out
        _dest = SUB_DIR / f"task_{_task_id}"
        _dest.mkdir(parents=True, exist_ok=True)
        (_dest / "output.json").write_text(json.dumps(_out, indent=2))
        (_dest / "reasoning.md").write_text(build_reasoning(_out))
    # Show task 1.3 (A375 48h) as a worked example, ordered by predicted rank.
    example = (
        pl.DataFrame(filled["1.3"]["rankings"])
        .select("rank", "condition", "concentration", "viability_pct")
        .sort("rank")
    )
    mo.md(
        f"## Built {len(filled)} submissions -> `{SUB_DIR.relative_to(NOTEBOOK_DIR.parent)}/task_*/`\n\n"
        "Each task dir holds `output.json` (the filled template) and `reasoning.md` (the trace the "
        "Karman `submit_response` tool wants). Example: task 1.3 (A375, 48h), rank 1 = strongest effect:"
    )
    mo.ui.table(example, page_size=12)
    return (filled,)


@app.cell
def _(predictions):
    # Self-score: line each task's predicted viability up against the local oracle viability.
    oracle = pl.concat(
        [
            pl.read_csv(ORACLE[_tp])
            .select("cell_line", "condition", "viability_pct")
            .rename({"viability_pct": "oracle_viability_pct"})
            .with_columns(timepoint_h=pl.lit(_tp))
            for _tp in (24, 48)
        ]
    )
    joined = predictions.join(oracle, on=["condition", "cell_line", "timepoint_h"], how="inner")

    score_rows = []
    for _task_id, _line, _tp in SINGLE_TASKS:
        _sub = (
            joined
            if _task_id == "1.5"
            else joined.filter((pl.col("cell_line") == _line) & (pl.col("timepoint_h") == _tp))
        )
        _pred, _truth = _sub["pred_viability_pct"].to_list(), _sub["oracle_viability_pct"].to_list()
        score_rows.append(
            {
                "task": _task_id,
                "scope": "pooled (2 lines x 2 timepoints)" if _task_id == "1.5" else f"{_line} {_tp}h",
                "n": len(_pred),
                "spearman": round(spearmanr(_pred, _truth).statistic, 3),
                "kendall_tau_b": round(kendalltau(_pred, _truth).statistic, 3),
            }
        )
    scores = pl.DataFrame(score_rows)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    scores.write_csv(OUT_DIR / "scores.csv")
    joined.write_csv(OUT_DIR / "predicted_vs_oracle.csv")
    mo.md(
        "## Self-score against the local oracle\n\n"
        "Rank agreement of our predicted viability vs the measured oracle, per task (proxy for the "
        "benchmark's tie-band scorer). Higher = more of that task is recoverable from public data; "
        "the gap to 1.0 is the residual difficulty. **A real contestant cannot compute this** - the "
        "oracle is private; here it is our internal check."
    )
    mo.ui.table(scores, page_size=8)
    return joined, scores


@app.cell
def _(joined):
    chart = (
        alt.Chart(joined.with_columns(label=pl.format("{} {}h", pl.col("cell_line"), pl.col("timepoint_h"))))
        .mark_circle(size=90, opacity=0.8)
        .encode(
            x=alt.X("pred_viability_pct:Q", title="predicted % viability (PRISM-grounded)"),
            y=alt.Y("oracle_viability_pct:Q", title="oracle % viability (measured)"),
            color=alt.Color("label:N", title="line / timepoint"),
            tooltip=["condition", "cell_line", "timepoint_h", "pred_viability_pct", "oracle_viability_pct"],
        )
        .properties(width=420, height=360)
    )
    mo.md(
        "## Predicted vs measured viability\n\n"
        "Each point is one condition across all four per-task slices. A positive trend means the "
        "submission and the wet lab agree; the 24h points (predicted >55%) and 48h points separate "
        "along x by the kinetic factor."
    )
    mo.ui.altair_chart(chart)
    return


@app.cell
def _(filled, scores):
    a375_48 = scores.filter(pl.col("task") == "1.3")["spearman"][0]
    lox_48 = scores.filter(pl.col("task") == "1.4")["spearman"][0]
    pooled = scores.filter(pl.col("task") == "1.5")["spearman"][0]
    summary = {
        "description": (
            "A DepMap (PRISM Repurposing Secondary) grounded submission to nightshift single-agent "
            "tasks 1.1-1.5. Predicts per-condition % viability and rank from dose-matched PRISM "
            "fold-change plus a 24h/48h kinetic factor; self-scored against the local oracle."
        ),
        "numbers": {
            "tasks": [t[0] for t in SINGLE_TASKS],
            "spearman_A375_48h": a375_48,
            "spearman_LOXIMVI_48h": lox_48,
            "spearman_pooled_1_5": pooled,
        },
        "files": [f"submission/task_{t[0]}/{f}" for t in filled for f in ("output.json", "reasoning.md")]
        + ["scores.csv", "predicted_vs_oracle.csv"],
    }
    (OUT_DIR / "summary.json").write_text(json.dumps(summary, indent=2))
    mo.md(f"```json\n{json.dumps(summary, indent=2)}\n```")
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ## What this submission can and cannot do

    - **The ranking is the load-bearing prediction**, and it comes entirely from PRISM measured
      drug response in these two lines - the strongest public signal available. The grader uses
      `rank` when present, so the predicted % viability is a readable secondary.
    - **It cannot reorder drugs between 24h and 48h.** PRISM has one timepoint; our kinetic factor
      is a uniform scaling, so the within-line ordering is identical at both timepoints. The oracle
      *does* reorder (e.g. Panobinostat overtakes Sapanisertib by 48h) - that drug-specific kinetics
      is a real part of the challenge that public sensitivity data does not carry.
    - **Absolute viability is approximate.** PRISM is a ~5-day pooled screen, the oracle a 48h
      single-dose CellTiter-Glo, and the dose grid is coarse (4x steps), so the level can be off
      even where the ranking is right (Panobinostat is the clearest miss).
    - **LOXIMVI is the harder line** - its predicted ranking agrees with the oracle less than A375's.

    ## To extend

    - Swap the rank signal to dose-collapsed PRISM **AUC** and compare self-scores - AUC integrates
      the whole curve and ranked A375 better in earlier analysis, at the cost of a dose-specific level.
    - Encode a per-drug kinetic prior (BRAF/MEK inhibitors act slower than mTOR/HDAC) so the 24h and
      48h rankings can differ - the one thing this submission structurally cannot predict.
    - Ensemble PRISM with GDSC2 and CTD^2 (also in Breadbox for these lines) and submit the consensus.
    """)
    return


if __name__ == "__main__":
    app.run()
