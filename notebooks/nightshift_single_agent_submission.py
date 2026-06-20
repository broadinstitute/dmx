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

    # Run live on molab (WASM/Pyodide): patch the stdlib so `requests` works over browser fetch.
    try:  # no-op outside Pyodide
        import pyodide_http  # noqa: F401

        pyodide_http.patch_all()
    except Exception:
        pass

    NOTEBOOK_DIR = Path(__file__).resolve().parent
    if str(NOTEBOOK_DIR) not in sys.path:
        sys.path.insert(0, str(NOTEBOOK_DIR))

    # In the dmx repo this notebook reuses the catalog's hardened Breadbox primitives from nb02.
    # Shared as a STANDALONE notebook (e.g. on molab) that sibling is absent, so fall back to inline
    # definitions of bb_get / bb_post - a deliberate duplication that keeps the notebook self-contained.
    try:
        from nb02_dataset_discovery import bb_get, bb_post  # noqa: E402
    except ImportError:
        import time  # noqa: E402

        import requests  # noqa: E402

        _BB_BASE = "https://depmap.org/portal/breadbox"
        _BB_HEADERS = {"Accept": "application/json", "User-Agent": "dmx/0.1 (+https://github.com/broadinstitute/dmx)"}

        def _bb(method, endpoint, *, params=None, body=None):
            # Identifying User-Agent (the portal 403s the default) + retry on transient 5xx / connection errors.
            for attempt in range(6):
                try:
                    resp = requests.request(
                        method,
                        f"{_BB_BASE}/{endpoint.lstrip('/')}",
                        params=params,
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

        def bb_get(endpoint, params=None):
            return _bb("GET", endpoint, params=params)

        def bb_post(endpoint, body):
            return _bb("POST", endpoint, body=body)

    # The answer-key-free task templates are fetched from the live Karman benchmark server, so this
    # notebook is fully self-contained - no local files, no sibling checkout.
    TASKS_URL = "https://mcp.karmanai.org/tasks"

    OUT_DIR = NOTEBOOK_DIR.parent / "data" / "processed" / "nightshift_single_agent_submission"
    SUB_DIR = OUT_DIR / "submission"

    # The two melanoma lines and their DepMap model ids.
    CELL_LINES = {"A375": "ACH-000219", "LOXIMVI": "ACH-000750"}

    # DepMap datasets (ids are stable enough to hardcode for a report; re-discover via nb02 if a 404 appears).
    PRISM_VIAB = (
        "576e1cb6-ac8d-4e29-bf15-0552c8665d72"  # PRISM Repurposing Secondary (Viability), dose-resolved log2 FC
    )
    PRISM_AUC = "07b7bda9-ae00-43b3-bca1-336b9607f8f5"  # PRISM Repurposing Secondary (AUC), dose-collapsed
    HOTSPOT = "a952ab7b-56c8-4aeb-872e-8ee02eeae042"  # Hotspot Mutations (value > 0 = mutated)

    # Kinetic factor: PRISM is a ~5-day endpoint; the tasks read at 24h / 48h. We treat PRISM's
    # killing as the 48h magnitude (factor 1.0) and attenuate the 24h prediction toward baseline.
    # A stated prior, not tuned to any measured result; it only shifts level and how 24h/48h interleave.
    KINETIC_FACTOR = {24: 0.55, 48: 1.0}

    # The 12 single agents: fixed nightshift dose, target node, and PRISM label (Sapanisertib is
    # screened under its code MLN0128). Doses taken verbatim from the task prompts.
    PANEL = [
        {"drug": "Panobinostat", "target": "pan-HDAC", "prism": "PANOBINOSTAT", "dose_uM": 0.05},
        {"drug": "Trametinib", "target": "MEK", "prism": "TRAMETINIB", "dose_uM": 0.01},
        {"drug": "Dabrafenib", "target": "BRAF", "prism": "DABRAFENIB", "dose_uM": 0.1},
        {"drug": "Encorafenib", "target": "BRAF", "prism": "ENCORAFENIB", "dose_uM": 0.1},
        {"drug": "Cobimetinib", "target": "MEK", "prism": "COBIMETINIB", "dose_uM": 0.1},
        {"drug": "Binimetinib", "target": "MEK", "prism": "BINIMETINIB", "dose_uM": 0.1},
        {"drug": "TAK-733", "target": "MEK", "prism": "TAK-733", "dose_uM": 0.03},
        {"drug": "Vemurafenib", "target": "BRAF", "prism": "VEMURAFENIB", "dose_uM": 1.0},
        {"drug": "Regorafenib", "target": "multi-kinase", "prism": "REGORAFENIB", "dose_uM": 5.0},
        {"drug": "Sapanisertib", "target": "mTOR", "prism": "MLN0128", "dose_uM": 0.5},
        {"drug": "Capivasertib", "target": "AKT", "prism": "CAPIVASERTIB", "dose_uM": 5.0},
        {"drug": "Alpelisib", "target": "PI3K", "prism": "ALPELISIB", "dose_uM": 5.0},
    ]
    DRUGS = [c["drug"] for c in PANEL]

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
    # Predicting drug response in BRAF-mutant melanoma, from public data alone

    **A Night Shift submission, grounded end to end in [DepMap](https://depmap.org).**

    The [Night Shift / Karman](https://karmanai.org/) benchmark runs a real melanoma viability
    experiment - two cell lines, 12 compounds at fixed doses, read at 24h and 48h - and asks an agent
    to predict the ranking *before* the wet-lab readout exists, from public knowledge only. Tasks
    1.1-1.4 are A375 / LOXIMVI at 24h / 48h; 1.5 pools all 48 conditions.

    This notebook is that prediction, but it shows its work. Every number below is pulled live from
    DepMap's Breadbox API - the cell-line identities, the measured PRISM dose-response curves, and how
    these lines compare to the rest of the cancer cell-line panel - and only then turned into a ranked
    prediction. Nothing is hand-entered; the notebook is self-contained and re-runs from scratch.

    1. **The system** - confirm the two lines are BRAF-mutant melanoma, from DepMap.
    2. **The evidence** - the PRISM dose-response curves the prediction reads.
    3. **The context** - how sensitive these lines are versus ~700 others.
    4. **The prediction** - dose-matched read + a kinetic prior -> a ranked submission per task.
    5. **The reasoning** - the trace submitted with each task, inline.
    """)
    return


@app.function
def dataset_features(dataset_id: str) -> pl.DataFrame:
    """List the features (compound/gene labels) available in a Breadbox matrix dataset."""
    records = bb_get(f"datasets/features/{dataset_id}")
    return pl.DataFrame(records, infer_schema_length=10_000) if records else pl.DataFrame()


@app.function
def matrix_feature(dataset_id: str, feature_label: str, value_name: str = "value") -> pl.DataFrame:
    """Fetch one feature column (across all models) from a Breadbox matrix dataset."""
    raw = bb_post(f"datasets/matrix/{dataset_id}", {"features": [feature_label], "feature_identifier": "label"})
    values = raw.get(feature_label)
    if values is None and isinstance(raw, dict) and raw:
        values = next(iter(raw.values()))
    if not isinstance(values, dict):
        return pl.DataFrame()
    return pl.DataFrame({"depmap_id": list(values.keys()), value_name: list(values.values())})


@app.function
def two_line_values(dataset_id: str, feature_label: str) -> dict[str, float | None]:
    """Pull one matrix feature and keep only the two nightshift cell lines, keyed by name."""
    frame = matrix_feature(dataset_id, feature_label, value_name="value")
    by_id = dict(zip(frame["depmap_id"], frame["value"], strict=False)) if frame.height else {}
    return {name: by_id.get(model_id) for name, model_id in CELL_LINES.items()}


@app.function
def model_metadata(columns: list[str]) -> dict:
    """Fetch cell-line metadata columns from Breadbox (column-oriented {col: {model_id: value}})."""
    return bb_post("datasets/tabular/depmap_model_metadata", {"columns": columns})


@app.function
def nearest_dose(available: list[float], target_uM: float) -> float:
    """Pick the screened concentration nearest (in log10 uM) to a target; ties break to the lower dose."""
    return min(available, key=lambda d: (round(abs(math.log10(d) - math.log10(target_uM)), 6), d))


@app.function
def predicted_viability_pct(prism_log2fc: float, timepoint_h: int) -> float:
    """Predicted % viability at a timepoint from a PRISM log2 fold-change (2**log2fc, scaled by kinetics)."""
    effect = max(0.0, 1.0 - 2.0**prism_log2fc)  # growth (>baseline) counts as zero killing
    return max(0.0, min(110.0, 100.0 * (1.0 - KINETIC_FACTOR[timepoint_h] * effect)))


@app.function
def conc_label(dose_uM: float) -> str:
    """Human concentration string, e.g. 0.05 -> '50 nM', 5.0 -> '5 uM'."""
    return f"{dose_uM * 1000:g} nM" if dose_uM < 1 else f"{dose_uM:g} uM"


@app.function
def prism_dose_grid() -> pl.DataFrame:
    """Full PRISM dose-response: every screened concentration x compound x line, plus the matched read.

    Returns long rows (drug, target, cell_line, dose_uM, viability_pct, is_matched) where is_matched
    flags the screened point nearest each drug's nightshift dose - the single value the prediction reads.
    """
    features = dataset_features(PRISM_VIAB)
    labels = features["label"].to_list() if features.height else []
    rows = []
    for c in PANEL:
        pattern = re.compile(rf"^{re.escape(c['prism'])} ([\d.]+) uM$")
        doses = sorted(float(m.group(1)) for label in labels if (m := pattern.match(label)))
        if not doses:
            continue
        matched = nearest_dose(doses, c["dose_uM"])
        for dose in doses:
            label = next(la for la in labels if pattern.match(la) and float(pattern.match(la).group(1)) == dose)
            for line, log2fc in two_line_values(PRISM_VIAB, label).items():
                if log2fc is None:
                    continue
                rows.append(
                    {
                        "drug": c["drug"],
                        "target": c["target"],
                        "cell_line": line,
                        "dose_uM": dose,
                        "viability_pct": round(100.0 * 2.0**log2fc, 1),
                        "log2fc": log2fc,
                        "is_matched": dose == matched,
                    }
                )
    return pl.DataFrame(rows)


@app.function
def panel_predictions(grid: pl.DataFrame | None = None) -> pl.DataFrame:
    """Predicted % viability per (condition, cell_line, timepoint) from the matched PRISM read + kinetics."""
    if grid is None:
        grid = prism_dose_grid()
    matched = grid.filter(pl.col("is_matched"))
    rows = []
    for row in matched.iter_rows(named=True):
        for tp in (24, 48):
            rows.append(
                {
                    "condition": row["drug"],
                    "cell_line": row["cell_line"],
                    "timepoint_h": tp,
                    "pred_viability_pct": round(predicted_viability_pct(row["log2fc"], tp), 1),
                }
            )
    return pl.DataFrame(rows)


@app.function
def auc_percentiles() -> pl.DataFrame:
    """For each compound, where A375 / LOXIMVI sensitivity (PRISM AUC) falls among all screened lines."""
    rows = []
    for c in PANEL:
        frame = matrix_feature(PRISM_AUC, c["prism"], value_name="auc")
        by_id = dict(zip(frame["depmap_id"], frame["auc"], strict=False)) if frame.height else {}
        values = sorted(v for v in by_id.values() if v is not None)
        if not values:
            continue
        for line, model_id in CELL_LINES.items():
            x = by_id.get(model_id)
            if x is None:
                continue
            pctile = round(100.0 * sum(1 for a in values if a < x) / len(values), 1)
            rows.append(
                {
                    "drug": c["drug"],
                    "target": c["target"],
                    "cell_line": line,
                    "auc": round(x, 3),
                    "pctile": pctile,
                    "n_lines": len(values),
                }
            )
    return pl.DataFrame(rows)


@app.function
def load_task_template(task_id: str) -> dict:
    """Fetch a task's answer-key-free output.json template from the live Karman benchmark server."""
    return bb_get_json(f"{TASKS_URL}/{task_id}/output.json")


@app.function
def bb_get_json(url: str) -> dict:
    """Plain GET of a JSON document (used for the Karman task templates, not Breadbox)."""
    import requests  # local import: in-repo, requests lives only in the inline fallback branch

    resp = requests.get(url, headers={"Accept": "application/json", "User-Agent": "dmx/0.1"}, timeout=30)
    resp.raise_for_status()
    return resp.json()


@app.function
def fill_submission(template: dict, value_of) -> dict:
    """Populate a task's output.json: rank + viability_pct from value_of(condition, line, tp) (rank 1 = strongest)."""
    out = copy.deepcopy(template)
    task_line, task_tp = template.get("cell_line"), template.get("timepoint_h")
    entries = out["rankings"]
    values = [value_of(e["condition"], e.get("cell_line", task_line), e.get("timepoint_h", task_tp)) for e in entries]
    for position, i in enumerate(sorted(range(len(entries)), key=lambda i: values[i]), start=1):
        entries[i]["rank"] = position
        entries[i]["viability_pct"] = round(values[i], 1)
    return out


@app.function
def combo_viability(condition: str, cell_line: str, timepoint_h: int, pred_via: dict) -> float:
    """Bliss-independence prediction for a 2-drug combo: product of the single-agent fractional viabilities."""
    a, b = (d.strip() for d in condition.split("+"))
    return round(pred_via[(a, cell_line, timepoint_h)] / 100.0 * pred_via[(b, cell_line, timepoint_h)], 1)


@app.function
def build_reasoning(filled: dict) -> str:
    """Generate the reasoning_md trace for a filled task - method, grounding, ranking, caveats."""
    entries = sorted(filled["rankings"], key=lambda e: e["rank"])
    pooled = "cell_line" in entries[0]
    if pooled:
        head = "| rank | condition | cell line | timepoint | predicted viability % |\n|---|---|---|---|---|"
        body = "\n".join(
            f"| {e['rank']} | {e['condition']} | {e['cell_line']} | {e['timepoint_h']}h | {e['viability_pct']} |"
            for e in entries
        )
    else:
        head = "| rank | condition | conc | predicted viability % |\n|---|---|---|---|"
        body = "\n".join(
            f"| {e['rank']} | {e['condition']} | {e['concentration']} | {e['viability_pct']} |" for e in entries
        )
    return f"""# Reasoning - Night Shift task {filled["task"]}

{filled["description"]}

## Method (public-data lookup)

This ranking is grounded entirely in **DepMap PRISM Repurposing Secondary**, a measured small-molecule
viability screen, for the exact cell lines in this task (A375 = ACH-000219, LOXIMVI = ACH-000750). For
each compound I read PRISM's log2 fold-change vs DMSO at the screened concentration nearest its stated
dose, convert to fractional viability (`2**log2fc`), and apply a fixed 24h/48h kinetic factor (PRISM is
a ~5-day endpoint; 24h is attenuated toward baseline). Conditions are ranked by predicted viability,
rank 1 = strongest effect (lowest viability). All 12 compounds resolve in PRISM for both lines
(Sapanisertib via its code MLN0128).

## Predicted ranking

{head}
{body}

## Caveats

- The ranking is driven by measured PRISM sensitivity; the absolute viability % is approximate
  (different assay, timepoint, and a coarse 4x dose grid).
- PRISM has a single timepoint, so this method cannot reorder drugs between 24h and 48h.
- HDAC-inhibitor potency is under-represented by one dose-matched point (polypharmacology).
"""


@app.function
def build_combo_reasoning(filled: dict) -> str:
    """Reasoning trace for a combination task - the Bliss-independence model over the single-agent reads."""
    entries = sorted(filled["rankings"], key=lambda e: e["rank"])
    pooled = "cell_line" in entries[0]
    if pooled:
        head = "| rank | combination | cell line | predicted viability % |\n|---|---|---|---|"
        body = "\n".join(
            f"| {e['rank']} | {e['condition']} | {e['cell_line']} | {e['viability_pct']} |" for e in entries
        )
    else:
        head = "| rank | combination | conc | predicted viability % |\n|---|---|---|---|"
        body = "\n".join(
            f"| {e['rank']} | {e['condition']} | {e['concentration']} | {e['viability_pct']} |" for e in entries
        )
    return f"""# Reasoning - Night Shift task {filled["task"]}

{filled["description"]}

## Method (public-data lookup + Bliss independence)

Each combination is two of the 12 single agents at their stated doses. I predict each single agent's
48h viability from DepMap PRISM (exactly as in tasks 1.x), then combine under the **Bliss independence**
null model: the combined fractional viability is the product of the singles (v_AB = v_A x v_B). This is
the no-interaction expectation - it cannot predict synergy or antagonism, which are not derivable from
single-agent public data. Combinations are ranked by predicted viability, rank 1 = strongest effect.

## Predicted ranking

{head}
{body}

## Caveats

- Bliss independence is a null model: real synergy/antagonism will shift the true ranking.
- Inherits the single-agent caveats (assay / timepoint / dose-grid approximation).
"""


@app.cell
def _():
    # Section 1 - confirm the system from public data: lineage, disease, BRAF status.
    _meta = model_metadata(["CellLineName", "OncotreeLineage", "OncotreePrimaryDisease"])
    _braf = two_line_values(HOTSPOT, "BRAF")
    card = pl.DataFrame(
        [
            {
                "cell_line": _name,
                "DepMap id": _mid,
                "name": _meta["CellLineName"].get(_mid),
                "lineage": _meta["OncotreeLineage"].get(_mid),
                "disease": _meta["OncotreePrimaryDisease"].get(_mid),
                "BRAF": "hotspot-mutated" if (_braf.get(_name) or 0) > 0 else "wild-type",
            }
            for _name, _mid in CELL_LINES.items()
        ]
    )
    mo.vstack(
        [
            mo.md(
                "## 1. The system, confirmed from DepMap\n\n"
                "Both lines resolve in DepMap as **skin / melanoma** carrying a **BRAF hotspot mutation** - "
                "the BRAF-V600E melanoma setting the benchmark describes, verified rather than assumed."
            ),
            mo.ui.table(card, page_size=2),
        ]
    )
    return


@app.cell
def _():
    # Section 2 - the evidence: the full PRISM dose-response the prediction reads (the showpiece).
    dose_grid = prism_dose_grid()
    _curve = (
        alt.Chart(dose_grid)
        .mark_line(point=alt.OverlayMarkDef(size=22), strokeWidth=1.6)
        .encode(
            x=alt.X("dose_uM:Q", scale=alt.Scale(type="log"), title="dose (uM, log)"),
            y=alt.Y("viability_pct:Q", title="viability %", scale=alt.Scale(domain=[0, 120])),
            color=alt.Color("cell_line:N", title="line"),
            tooltip=["drug", "cell_line", "dose_uM", "viability_pct"],
        )
    )
    _matched = (
        alt.Chart(dose_grid)
        .transform_filter(alt.datum.is_matched)
        .mark_point(size=130, filled=True, opacity=0.9, shape="diamond")
        .encode(
            x=alt.X("dose_uM:Q", scale=alt.Scale(type="log")),
            y="viability_pct:Q",
            color=alt.Color("cell_line:N"),
            tooltip=["drug", "cell_line", "dose_uM", "viability_pct"],
        )
    )
    grid_chart = (
        alt.layer(_curve, _matched).properties(width=190, height=140).facet(facet="drug:N", columns=4, title=None)
    )
    mo.vstack(
        [
            mo.md(
                "## 2. The evidence: PRISM dose-response curves\n\n"
                "Each panel is one compound's measured 8-point viability curve in A375 and LOXIMVI (PRISM "
                "Repurposing Secondary). The **diamond** marks the screened concentration nearest the nightshift "
                "dose - the single value the prediction reads. Steep early-dropping curves (the MEK/HDAC/mTOR "
                "drugs) are the strong killers; flat curves (Alpelisib, Capivasertib) barely move."
            ),
            mo.ui.altair_chart(grid_chart),
        ]
    )
    return (dose_grid,)


@app.cell
def _():
    # An interactive focus: pick one compound and see its curve large, with the nightshift dose called out.
    drug_picker = mo.ui.dropdown(options=DRUGS, value="Panobinostat", label="Compound")
    drug_picker
    return (drug_picker,)


@app.cell
def _(dose_grid, drug_picker):
    _sel = drug_picker.value
    _one = dose_grid.filter(pl.col("drug") == _sel)
    _ns = next(c["dose_uM"] for c in PANEL if c["drug"] == _sel)
    _line = (
        alt.Chart(_one)
        .mark_line(point=True, strokeWidth=2)
        .encode(
            x=alt.X("dose_uM:Q", scale=alt.Scale(type="log"), title="dose (uM, log)"),
            y=alt.Y("viability_pct:Q", title="viability %", scale=alt.Scale(domain=[0, 120])),
            color=alt.Color("cell_line:N", title="line"),
            tooltip=["cell_line", "dose_uM", "viability_pct"],
        )
    )
    _rule = alt.Chart(pl.DataFrame({"d": [_ns]})).mark_rule(strokeDash=[5, 4], color="#888").encode(x="d:Q")
    chart_one = alt.layer(_line, _rule).properties(width=480, height=300)
    mo.vstack(
        [
            mo.md(f"### {_sel} - dose-response (dashed line = nightshift dose, {conc_label(_ns)})"),
            mo.ui.altair_chart(chart_one),
        ]
    )
    return


@app.cell
def _():
    # Section 3 - sensitivity in context: where these two lines fall among all screened lines.
    context = auc_percentiles()
    _ctx_chart = (
        alt.Chart(context)
        .mark_circle(size=140, opacity=0.85)
        .encode(
            x=alt.X(
                "pctile:Q",
                title="sensitivity percentile vs all DepMap lines (lower = more sensitive)",
                scale=alt.Scale(domain=[0, 100]),
            ),
            y=alt.Y("drug:N", sort=alt.EncodingSortField(field="pctile", op="min"), title=None),
            color=alt.Color("cell_line:N", title="line"),
            tooltip=["drug", "cell_line", "pctile", "auc", "n_lines"],
        )
        .properties(width=520, height=320)
    )
    _median = context["pctile"].median()
    mo.vstack(
        [
            mo.md(
                "## 3. Are these lines special? Sensitivity in context\n\n"
                f"For each compound, the percentile of A375 / LOXIMVI sensitivity (PRISM AUC) among all "
                f"~{int(context['n_lines'][0])} screened cell lines - **lower = more sensitive than most**. The "
                "BRAF/MEK/mTOR drugs land in the low single-to-double-digit percentiles: these melanoma lines are "
                "among the most sensitive in the whole panel, which is exactly the prior the prediction leans on. "
                f"(Median percentile across the panel: {_median:.0f}.)"
            ),
            mo.ui.altair_chart(_ctx_chart),
        ]
    )
    return


@app.cell
def _(dose_grid):
    # Section 4 - from the matched read to a ranked prediction.
    predictions = panel_predictions(dose_grid)
    pred_via = {
        (r["condition"], r["cell_line"], r["timepoint_h"]): r["pred_viability_pct"]
        for r in predictions.iter_rows(named=True)
    }
    _pchart = (
        alt.Chart(predictions)
        .mark_circle(size=80, opacity=0.85)
        .encode(
            x=alt.X("pred_viability_pct:Q", title="predicted % viability (lower = stronger)"),
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
                "## 4. The prediction\n\n"
                "Read PRISM at the matched dose (the diamond above), convert to % viability, apply the 24h/48h "
                "kinetic factor. Lower = stronger predicted effect; 24h sits above 48h by the kinetic factor."
            ),
            mo.ui.altair_chart(_pchart),
        ]
    )
    return (pred_via,)


@app.cell
def _(pred_via):
    # Build + write the five populated output.json submissions (+ reasoning) from the live templates.
    SUB_DIR.mkdir(parents=True, exist_ok=True)
    filled = {}
    for _task_id, _line, _tp in SINGLE_TASKS:
        _out = fill_submission(load_task_template(_task_id), lambda c, line, tp: pred_via[(c, line, tp)])
        filled[_task_id] = _out
        _dest = SUB_DIR / f"task_{_task_id}"
        _dest.mkdir(parents=True, exist_ok=True)
        (_dest / "output.json").write_text(json.dumps(_out, indent=2))
        (_dest / "reasoning.md").write_text(build_reasoning(_out))
    example = (
        pl.DataFrame(filled["1.3"]["rankings"])
        .select("rank", "condition", "concentration", "viability_pct")
        .sort("rank")
    )
    mo.vstack(
        [
            mo.md(
                f"### The five submissions -> `{SUB_DIR.relative_to(NOTEBOOK_DIR.parent)}/task_*/`\n\n"
                "Each task dir holds the filled `output.json` and a `reasoning.md` trace. Example: task 1.3 "
                "(A375, 48h), rank 1 = strongest effect."
            ),
            mo.ui.table(example, page_size=12),
        ]
    )
    return (filled,)


@app.cell
def _(filled):
    # Section 5 - the reasoning trace submitted with each task, rendered inline.
    _traces = {f"Task {_tid}": mo.md(build_reasoning(_out)) for _tid, _out in filled.items()}
    mo.vstack(
        [
            mo.md(
                "## 5. Reasoning traces (inline)\n\n"
                "The rationale submitted with each task - expand one for its method, ranking, and caveats."
            ),
            mo.accordion(_traces),
        ]
    )
    return


@app.cell
def _(pred_via):
    # Section 6 - combinations (tasks 2.1-2.3): Bliss independence over the single-agent predictions.
    combo_filled = {}
    for _tid in ("2.1", "2.2", "2.3"):
        _out = fill_submission(load_task_template(_tid), lambda c, line, tp: combo_viability(c, line, tp, pred_via))
        combo_filled[_tid] = _out
        _dest = SUB_DIR / f"task_{_tid}"
        _dest.mkdir(parents=True, exist_ok=True)
        (_dest / "output.json").write_text(json.dumps(_out, indent=2))
        (_dest / "reasoning.md").write_text(build_combo_reasoning(_out))
    combo_df = pl.DataFrame(
        [
            {
                "task": _tid,
                "combination": _e["condition"],
                "cell_line": _e.get("cell_line", _o.get("cell_line")),
                "rank": _e["rank"],
                "pred_viability_pct": _e["viability_pct"],
            }
            for _tid, _o in combo_filled.items()
            for _e in _o["rankings"]
        ]
    )
    _chart = (
        alt.Chart(combo_df.filter(pl.col("task") != "2.3"))
        .mark_bar(opacity=0.85)
        .encode(
            x=alt.X("pred_viability_pct:Q", title="predicted combo viability % (lower = stronger)"),
            y=alt.Y("combination:N", sort="x", title=None),
            color=alt.Color("cell_line:N", title="line"),
            tooltip=["combination", "cell_line", "rank", "pred_viability_pct"],
        )
        .properties(width=360, height=180)
        .facet(row=alt.Row("cell_line:N", title=None))
    )
    mo.vstack(
        [
            mo.md(
                "## 6. Combinations (tasks 2.1-2.3)\n\n"
                "Each combination is two of the 12 single agents at their stated doses. The predicted combo "
                "viability is the **Bliss-independence** product of the two single-agent predictions - the "
                "no-interaction expectation, so two individually-potent drugs rank strongest. True synergy "
                "(e.g. between complementary pathway arms) is not capturable from single-agent data; it is "
                "flagged in each combination's reasoning trace."
            ),
            mo.ui.altair_chart(_chart),
        ]
    )
    return (combo_filled,)


@app.cell
def _(combo_filled):
    # The combination reasoning traces, inline.
    _traces = {f"Task {_tid}": mo.md(build_combo_reasoning(_o)) for _tid, _o in combo_filled.items()}
    mo.vstack([mo.md("### Combination reasoning traces (tasks 2.1-2.3)"), mo.accordion(_traces)])
    return


@app.cell
def _(combo_filled, filled):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    _all_tasks = list(filled) + list(combo_filled)
    summary = {
        "description": (
            "A DepMap (PRISM Repurposing Secondary) grounded report + submission to nightshift tasks 1.1-1.5 "
            "(single agents) and 2.1-2.3 (combinations, via Bliss independence). Confirms the cell lines, shows "
            "the measured dose-response curves and sensitivity context, then predicts viability + rank. Public data only."
        ),
        "numbers": {"tasks": _all_tasks, "compounds": len(PANEL)},
        "files": [f"submission/task_{_t}/{_f}" for _t in _all_tasks for _f in ("output.json", "reasoning.md")],
    }
    (OUT_DIR / "summary.json").write_text(json.dumps(summary, indent=2))
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ## What this report shows, and its limits

    - **The prediction is evidence-backed, not a guess.** Sections 1-3 are all measured public data: the
      lines are BRAF-mutant melanoma, the dose-response curves are real PRISM measurements, and these lines
      are among the most drug-sensitive in DepMap. The ranking in section 4 is a direct read of those curves.
    - **It cannot reorder 24h vs 48h.** PRISM has one timepoint, so the kinetic factor is a uniform scaling;
      drug-specific kinetics (some compounds act faster) is a known blind spot.
    - **Absolute viability is approximate** - PRISM is a ~5-day pooled screen on a coarse 4x dose grid, not a
      48h single-dose CellTiter-Glo - and HDAC-inhibitor potency is under-represented by one matched point.

    ## To extend

    - Layer GDSC2 and CTD^2 dose-response (also in Breadbox for these lines) and submit a cross-screen consensus.
    - Add the CRISPR dependency of each drug's target gene (A375 is in the Chronos panel) as a mechanistic prior.
    - Encode a per-drug kinetic prior so the 24h and 48h rankings can differ - the one thing this cannot predict.
    """)
    return


if __name__ == "__main__":
    app.run()
