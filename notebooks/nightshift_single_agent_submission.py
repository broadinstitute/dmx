# /// script
# requires-python = ">=3.12,<3.14"
# dependencies = [
#     "altair==5.5.0",
#     "marimo",
#     "polars==1.40.1",
#     "requests==2.32.5",
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

    NOTEBOOK_DIR = Path(__file__).resolve().parent
    if str(NOTEBOOK_DIR) not in sys.path:
        sys.path.insert(0, str(NOTEBOOK_DIR))

    # In the dmx repo this notebook reuses the catalog's Breadbox helpers from nb02/nb03. Shared as a
    # STANDALONE notebook (e.g. on molab) those siblings are absent, so fall back to inline definitions
    # of the same two helpers - a deliberate duplication that keeps the notebook self-contained.
    try:
        from nb02_dataset_discovery import dataset_features  # noqa: E402
        from nb03_gene_dependency_profile import matrix_feature  # noqa: E402
    except ImportError:
        import time  # noqa: E402

        import requests  # noqa: E402

        _BB_BASE = "https://depmap.org/portal/breadbox"
        _BB_HEADERS = {"Accept": "application/json", "User-Agent": "dmx/0.1 (+https://github.com/broadinstitute/dmx)"}

        def _bb(method, endpoint, *, body=None):
            # Identifying User-Agent (the portal 403s the default) + retry on transient 5xx / connection errors.
            for attempt in range(6):
                try:
                    resp = requests.request(
                        method,
                        f"{_BB_BASE}/{endpoint.lstrip('/')}",
                        json=body,
                        headers=_BB_HEADERS,
                        timeout=120 if method == "POST" else 60,
                    )
                    if resp.status_code < 500:
                        resp.raise_for_status()
                        return resp.json()
                except requests.RequestException:
                    pass
                time.sleep(3.0 * (attempt + 1))
            raise RuntimeError(f"Breadbox {method} {endpoint} failed after retries")

        def dataset_features(dataset_id):
            """Inline fallback for nb02: list the features (compound/gene labels) in a matrix dataset."""
            records = _bb("GET", f"datasets/features/{dataset_id}")
            return pl.DataFrame(records, infer_schema_length=10_000) if records else pl.DataFrame()

        def matrix_feature(dataset_id, feature_label, value_name="value"):
            """Inline fallback for nb03: fetch one feature column across models from a matrix dataset."""
            raw = _bb(
                "POST",
                f"datasets/matrix/{dataset_id}",
                body={"features": [feature_label], "feature_identifier": "label"},
            )
            values = raw.get(feature_label)
            if values is None and isinstance(raw, dict) and raw:
                values = next(iter(raw.values()))
            if not isinstance(values, dict):
                return pl.DataFrame()
            return pl.DataFrame({"depmap_id": list(values.keys()), value_name: list(values.values())})

    # The nightshift benchmark task prompts live next to this catalog. We read ONLY the answer-key-free
    # templates here - never the oracle. Scoring is the grader's job (see nightshift_single_agent_eval).
    NIGHTSHIFT = NOTEBOOK_DIR.parent.parent / "nightshift-dev"
    TASKS_DIR = NIGHTSHIFT / "data" / "raw" / "tasks"

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

    It writes the five populated `output.json` submissions plus a `reasoning.md` trace for each,
    using **only public data - no oracle access** - exactly as the task requires. Scoring against
    the wet-lab ground truth is a separate, organizer-side step (`nightshift_single_agent_eval`);
    a real contestant cannot do it, so it does not belong here.
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


@app.function
def prism_dose_signal() -> pl.DataFrame:
    """PRISM dose-resolved log2 fold-change at the dose nearest each panel drug's task dose, both lines.

    Importable so the evaluation notebook reuses the exact same engine - no re-implementation, no
    dependence on the (gitignored) written submission files.
    """
    viab_features = dataset_features(PRISM_VIAB)
    rows = []
    for c in PANEL:
        pick = nearest_dose_label(viab_features, c["prism"], c["dose_uM"])
        label, picked_uM = pick if pick else (None, None)
        vals = two_line_values(PRISM_VIAB, label) if label else {k: None for k in CELL_LINES}
        for line, log2fc in vals.items():
            rows.append(
                {
                    "drug": c["drug"],
                    "cell_line": line,
                    "dose_uM": c["dose_uM"],
                    "prism_dose_uM": picked_uM,
                    "prism_log2fc": log2fc,
                    "prism_viab_pct": round(100.0 * 2.0**log2fc, 1) if log2fc is not None else None,
                }
            )
    return pl.DataFrame(rows)


@app.function
def panel_predictions(signal: pl.DataFrame | None = None) -> pl.DataFrame:
    """Predicted % viability per (condition, cell_line, timepoint) from the PRISM signal + kinetic factor."""
    if signal is None:
        signal = prism_dose_signal()
    rows = []
    for row in signal.iter_rows(named=True):
        if row["prism_log2fc"] is None:
            continue
        for tp in (24, 48):
            rows.append(
                {
                    "condition": row["drug"],
                    "cell_line": row["cell_line"],
                    "timepoint_h": tp,
                    "pred_viability_pct": round(predicted_viability_pct(row["prism_log2fc"], tp), 1),
                }
            )
    return pl.DataFrame(rows)


@app.cell
def _():
    # Engine: PRISM viability at the dose nearest each drug's task dose, for both lines.
    prism_signal = prism_dose_signal()
    mo.vstack(
        [
            mo.md(
                "## The PRISM signal\n\n"
                "Dose-resolved PRISM log2 fold-change at the curve point nearest each task dose "
                "(`prism_dose_uM`), and the implied fractional viability. All 12/12 compounds resolve for "
                "both lines. This single number per (drug, line) is the only measured input to the prediction."
            ),
            mo.ui.table(prism_signal.sort(["cell_line", "prism_log2fc"]), page_size=12),
        ]
    )
    return (prism_signal,)


@app.cell
def _(prism_signal):
    # Predicted % viability for every (drug, line, timepoint) - reuse the already-pulled signal.
    predictions = panel_predictions(prism_signal)
    pred_via = {
        (r["condition"], r["cell_line"], r["timepoint_h"]): r["pred_viability_pct"]
        for r in predictions.iter_rows(named=True)
    }
    mo.vstack(
        [
            mo.md(
                "## Predicted viability (the submission's raw prediction)\n\n"
                "One predicted % viability per condition x cell line x timepoint. The kinetic factor makes "
                "24h less killed than 48h; it cannot reorder drugs within a (line, timepoint) - that order "
                "is fixed by the PRISM dose signal."
            ),
            mo.ui.table(predictions.sort(["cell_line", "timepoint_h", "pred_viability_pct"]), page_size=12),
        ]
    )
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
    mo.vstack(
        [
            mo.md(
                f"## Built {len(filled)} submissions -> `{SUB_DIR.relative_to(NOTEBOOK_DIR.parent)}/task_*/`\n\n"
                "Each task dir holds `output.json` (the filled template) and `reasoning.md` (the trace the "
                "Karman `submit_response` tool wants). Example: task 1.3 (A375, 48h), rank 1 = strongest effect:"
            ),
            mo.ui.table(example, page_size=12),
        ]
    )
    return (filled,)


@app.cell
def _(filled):
    # Render each task's reasoning trace INLINE, so the rationale reads directly in the notebook /
    # molab. This is the same text written to each task's reasoning.md and submitted with it.
    _traces = {f"Task {_tid}": mo.md(build_reasoning(_out)) for _tid, _out in filled.items()}
    mo.vstack(
        [
            mo.md(
                "## Reasoning traces (inline)\n\n"
                "The rationale submitted with each task, rendered here so it reads without opening the "
                "`reasoning.md` files. Expand a task for its method, predicted ranking, and caveats."
            ),
            mo.accordion(_traces),
        ]
    )
    return


@app.cell
def _(predictions):
    # The submission's own prediction - no ground truth: predicted viability per drug, by line/timepoint.
    chart = (
        alt.Chart(predictions)
        .mark_circle(size=90, opacity=0.85)
        .encode(
            x=alt.X("pred_viability_pct:Q", title="predicted % viability (PRISM-grounded; lower = stronger)"),
            y=alt.Y("condition:N", sort="x", title=None),
            color=alt.Color("timepoint_h:N", title="timepoint (h)"),
            tooltip=["condition", "cell_line", "timepoint_h", "pred_viability_pct"],
        )
        .properties(width=300, height=320)
        .facet(column=alt.Column("cell_line:N", title=None))
    )
    mo.vstack(
        [
            mo.md(
                "## Predicted viability across the panel\n\n"
                "The submission's prediction only (no ground truth). Lower = stronger predicted effect; "
                "24h sits higher than 48h by the kinetic factor. Drugs sharing a PRISM signal cluster together."
            ),
            mo.ui.altair_chart(chart),
        ]
    )
    return


@app.cell
def _(filled):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    summary = {
        "description": (
            "A DepMap (PRISM Repurposing Secondary) grounded submission to nightshift single-agent "
            "tasks 1.1-1.5. Predicts per-condition % viability and rank from dose-matched PRISM "
            "fold-change plus a 24h/48h kinetic factor. Public data only - no oracle access."
        ),
        "numbers": {
            "tasks": [t[0] for t in SINGLE_TASKS],
            "compounds": len(PANEL),
            "conditions_task_1_5": len(filled["1.5"]["rankings"]),
        },
        "files": [f"submission/task_{t[0]}/{f}" for t in filled for f in ("output.json", "reasoning.md")],
    }
    (OUT_DIR / "summary.json").write_text(json.dumps(summary, indent=2))
    mo.md(f"```json\n{json.dumps(summary, indent=2)}\n```")
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ## What this submission can and cannot do

    - **The ranking is the load-bearing prediction**, grounded entirely in PRISM measured drug
      response in these two lines - the strongest public signal available. The grader uses `rank`
      when present, so the predicted % viability is a readable secondary.
    - **It cannot reorder drugs between 24h and 48h.** PRISM has one timepoint; the kinetic factor is
      a uniform scaling, so the within-line ordering is identical at both timepoints. Drug-specific
      kinetics (some compounds act faster than others) is a known blind spot of this method.
    - **Absolute viability is approximate.** PRISM is a ~5-day pooled screen on a coarse 4x dose grid,
      not a 48h single-dose CellTiter-Glo, so the predicted level can be off even where the ranking holds.
    - **HDAC potency is under-represented.** Panobinostat's dose-matched PRISM point understates its
      true effect (pan-HDAC polypharmacology), so its rank is the least certain - flagged in every trace.

    ## To extend

    - Swap the rank signal to dose-collapsed PRISM **AUC** (integrates the whole curve) and submit that variant.
    - Encode a per-drug kinetic prior (BRAF/MEK inhibitors act slower than mTOR/HDAC) so the 24h and
      48h rankings can differ - the one thing this submission structurally cannot predict.
    - Ensemble PRISM with GDSC2 and CTD^2 (also in Breadbox for these lines) and submit the consensus.
    """)
    return


if __name__ == "__main__":
    app.run()
