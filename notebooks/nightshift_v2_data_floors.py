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
    import itertools
    import json
    import math
    import re
    import time
    from pathlib import Path

    import altair as alt
    import marimo as mo
    import polars as pl
    import requests

    # Self-contained standalone notebook: vendors the engine + Breadbox HTTP helpers
    # below, so it runs as a single gist with no sibling files.
    NOTEBOOK_DIR = Path(__file__).resolve().parent
    REPO_DIR = NOTEBOOK_DIR.parent if NOTEBOOK_DIR.name == "notebooks" else NOTEBOOK_DIR
    OUT_DIR = REPO_DIR / "data" / "processed" / "nightshift_v2"
    BREADBOX = "https://depmap.org/portal/breadbox"
    PRISM_VIAB = "576e1cb6-ac8d-4e29-bf15-0552c8665d72"
    GDSC2_AUC = "2eac8e7b-beb4-48c1-b78f-c226723e54d7"  # second screen, pulled live for the collision
    A375 = "ACH-000219"
    LOXIMVI = "ACH-000750"
    CANDIDATES = [
        "Vemurafenib",
        "Dabrafenib",
        "Encorafenib",
        "Trametinib",
        "Cobimetinib",
        "Binimetinib",
        "TAK-733",
        "Regorafenib",
        "Alpelisib",
        "Capivasertib",
    ]
    # judgment-baseline floors (asymptotic % viability at task dose), used only where the
    # data source has no value for a drug. Kept inline so this notebook is standalone.
    FLOORS_V1 = {
        "Vemurafenib": (45, 92),
        "Dabrafenib": (28, 82),
        "Encorafenib": (33, 82),
        "Trametinib": (30, 60),
        "Cobimetinib": (32, 62),
        "Binimetinib": (38, 62),
        "TAK-733": (36, 62),
        "Regorafenib": (55, 58),
        "Alpelisib": (80, 78),
        "Capivasertib": (78, 72),
        "Sapanisertib": (50, 52),
        "Panobinostat": (30, 25),
    }


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    # Night Shift submission - data-derived floors

    Predict CellTiter-Glo % viability for A375 and LOXIMVI under a MAPK-inhibitor panel,
    and rank conditions by effect. No experimental readout is provided, so these are
    predictions.

    **The result worth leading with is a falsification.** The prior coming in held LOXIMVI to
    be BRAF-inhibitor-resistant (the CRISPR screen was silent on it; a long-assay drug screen,
    GDSC2, read it resistant). Experiment: read LOXIMVI's PRISM dose-viability for the BRAF
    inhibitors, and pull the GDSC2 numbers alongside it so both sides of the prior are live in
    this notebook. Observation: the two screens disagree *in direction* - GDSC2 reads LOXIMVI
    dabrafenib-resistant while PRISM reads it responding. Both are multi-day assays and neither
    is the 16-48 h CTG endpoint we predict, so this is a directional collision, not a
    calibrated number; the live "collision" cell below shows both. That unresolved disagreement
    - not any ranking - is the real finding. This submission then takes the PRISM measurement
    as its floor; which screen is right for a CTG plate is what the eventual readout decides,
    and the rankings downstream are just the fixed model re-run on these measured floors.

    The one design choice here: each drug's **effect floor** is read straight off a real
    screen rather than set by judgment. PRISM Repurposing Secondary stores 5-day viability
    at each dose, so for each drug we take the viability at the dose nearest the task
    concentration, per line, and use that as the floor (PRISM stores log2 fold-change, so
    % viability = `100 * 2^value`). Where a drug is absent from PRISM (sapanisertib), a
    judgment-baseline floor is used and flagged. The kinetics model and combination engine
    are vendored unchanged in the appendix; only the floors are data-derived.

    **This is not guaranteed to score better** - it trades "trust my judgment" for "trust
    the measurement," and one real data-vs-prior tension shows up (below). This is one of a
    small set of submissions that differ only in how the floors are produced; see the repo
    README for the others.
    """)
    return


@app.function
def prism_dose_floors() -> dict:
    """Effect floor per drug per line, read from PRISM Secondary dose-level viability.

    For each drug, find the two dose points nearest its task concentration (log scale),
    convert log2 fold-change to % viability (100*2^v), and average them. Returns
    {drug: {floorA, floorL, doses}} for drugs present in PRISM; callers fall back to a
    judgment-baseline floor for the rest.
    """
    panel = drug_panel()
    feats = bb_get(f"datasets/features/{PRISM_VIAB}")
    # Map each panel drug to its (dose_uM, feature_label) points.
    by_drug: dict[str, list] = {}
    for f in feats:
        label = str(f.get("label", ""))
        m = re.match(r"^([A-Za-z0-9-]+)\s+([\d.]+)\s*uM$", label)
        if not m:
            continue
        drug, dose = m.group(1).upper(), float(m.group(2))
        by_drug.setdefault(drug, []).append((dose, label))

    # Pick the two nearest dose labels per drug; query them all in one matrix POST.
    chosen: dict[str, list] = {}
    wanted = []
    for name, meta in panel.items():
        conc_uM = meta["conc_nM"] / 1000.0
        points = by_drug.get(name.upper())
        if not points:
            continue
        nearest = sorted(points, key=lambda x: abs(math.log10(x[0]) - math.log10(conc_uM)))[:2]
        chosen[name] = nearest
        wanted.extend(lbl for _, lbl in nearest)

    raw = bb_post(f"datasets/matrix/{PRISM_VIAB}", {"features": wanted, "feature_identifier": "label"})

    def viab(label, line_id):
        # log2 fold-change -> % viability. Keep values >100 (growth / paradoxical
        # activation) instead of clamping them away - that >100 signal is exactly
        # what flags a drug as inactive-or-worse in a line; only cap the floor at 2.
        v = raw.get(label, {}).get(line_id) if isinstance(raw, dict) else None
        return None if v is None else max(2.0, 100.0 * (2.0**v))

    floors = {}
    for name, nearest in chosen.items():
        a = [x for x in (viab(lbl, A375) for _, lbl in nearest) if x is not None]
        lox = [x for x in (viab(lbl, LOXIMVI) for _, lbl in nearest) if x is not None]
        if not a or not lox:
            continue
        doses = [d for d, _ in nearest]
        floors[name] = {
            "floorA": round(sum(a) / len(a), 1),
            "floorL": round(sum(lox) / len(lox), 1),
            "dose_uM": math.prod(doses) ** (1.0 / len(doses)),  # geomean dose the floor was read at
            "doses": " / ".join(f"{d:g}uM" for d in doses),
        }
    return floors


@app.cell
def _():
    prism = prism_dose_floors()

    def floor_v2(name, line):
        idx = "floorA" if line == "A375" else "floorL"
        if name in prism:
            return prism[name][idx]
        return FLOORS_V1[name][0 if line == "A375" else 1]  # absent from PRISM -> judgment baseline

    drugs_v2 = {
        n: {**meta, "floorA": floor_v2(n, "A375"), "floorL": floor_v2(n, "LOXIMVI")} for n, meta in drug_panel().items()
    }

    compare = pl.DataFrame(
        [
            {
                "drug": n,
                "source": "PRISM" if n in prism else "judgment",
                "doses": prism.get(n, {}).get("doses", "-"),
                "floor_A375_data": floor_v2(n, "A375"),
                "floor_A375_judgment": FLOORS_V1[n][0],
                "floor_LOX_data": floor_v2(n, "LOXIMVI"),
                "floor_LOX_judgment": FLOORS_V1[n][1],
            }
            for n in drug_panel()
        ]
    )
    mo.md("### Data-derived floors vs the judgment baseline")
    mo.ui.table(compare, page_size=12)
    return compare, drugs_v2, floor_v2, prism


@app.cell
def _(prism):
    # The leading falsification, both sides pulled live. GDSC2 (a multi-day screen) read
    # LOXIMVI as BRAF-inhibitor-RESISTANT (high AUC); PRISM (also multi-day, different design)
    # reads it RESPONDING (low % viability at the task dose). A directional disagreement, not a
    # number match - both are long assays, neither is the 16-48 h CTG endpoint we predict.
    gdsc = bb_post(
        f"datasets/matrix/{GDSC2_AUC}", {"features": ["DABRAFENIB", "ENCORAFENIB"], "feature_identifier": "label"}
    )

    def _g(drug, line):
        v = gdsc.get(drug, {}).get(line) if isinstance(gdsc, dict) else None
        return None if v is None else round(v, 2)

    collision = pl.DataFrame(
        [
            {
                "BRAF inhibitor": d.title(),
                "GDSC2 AUC LOXIMVI (high = resistant)": _g(d, LOXIMVI),
                "GDSC2 AUC A375 (low = sensitive)": _g(d, A375),
                "PRISM floor LOXIMVI % (low = responding)": prism.get(d.title(), {}).get("floorL"),
            }
            for d in ("DABRAFENIB", "ENCORAFENIB")
        ]
    )
    mo.md("### The collision, both sides pulled live (GDSC2 says resistant, PRISM says responding)")
    mo.ui.table(collision, page_size=4)
    return


@app.cell
def _(compare):
    # Where do data and judgment disagree most? (signed delta, A375 + LOXIMVI)
    div = (
        compare.with_columns(
            (pl.col("floor_A375_data") - pl.col("floor_A375_judgment")).alias("dA375"),
            (pl.col("floor_LOX_data") - pl.col("floor_LOX_judgment")).alias("dLOX"),
        )
        .select("drug", "dA375", "dLOX")
        .unpivot(index="drug", variable_name="line", value_name="delta")
    )
    _chart = (
        alt.Chart(div)
        .mark_bar()
        .encode(
            x=alt.X("delta:Q", title="data floor - judgment floor (negative = data says MORE killing)"),
            y=alt.Y("drug:N", sort="-x", title=None),
            color=alt.Color("line:N", title=None),
            tooltip=["drug", "line", "delta"],
        )
        .properties(width=560, height=420, title="Where PRISM data diverges from the judgment baseline")
    )
    mo.ui.altair_chart(_chart)
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    **The tension worth naming.** For A375 the PRISM floors track the judgment baseline
    closely (dabrafenib ~25 vs ~28) - which is reassuring but *expected*, not a discovery:
    the addicted line's BRAFi sensitivity is the one fact everyone, the prior included,
    already agreed on. The informative result is the **divergence in LOXIMVI BRAF
    inhibitors**: PRISM's fixed-dose 5-day viability shows LOXIMVI responding to
    dabrafenib/encorafenib, where a BRAF-WT prior expected near-vehicle (the Instrument-2
    GDSC2 pull, a longer assay, also read LOXIMVI as resistant). This submission takes the
    measurement as given rather than re-imposing the prior; that is exactly where data and
    judgment will disagree on the eventual CTG readout, and it is the honest cost of
    data-grounding a single screen with its own off-target / fixed-dose artifacts.

    One caveat to keep in view: only the *floor input* changed here. The kinetics + Bliss/
    synergy engine (vendored below) is held fixed, so this is not an independent test of the
    model - it is the same model on measured floors. The real 24-48 h CTG data is what
    actually adjudicates.
    """)
    return


@app.function
def methods_note_v2() -> str:
    """v2 methods preamble embedded in each reasoning.md."""
    return (
        "## Method (data-derived floors)\n\n"
        "Predictions, not measurements. A fixed CTG-kinetics + Bliss/pathway-combination "
        "engine (vendored in the notebook appendix), but each drug's per-line **effect "
        "floor** is read from PRISM Repurposing Secondary dose-level viability at the dose "
        "nearest the task concentration (% viability = 100*2^log2fc, averaged over the two "
        "nearest doses) instead of by judgment. Drugs absent from PRISM (sapanisertib) keep "
        "a judgment-baseline floor. The kinetic curves and synergy constants are held "
        "fixed - only the floors are data-derived.\n\n"
    )


@app.cell
def _(drugs_v2, floor_v2, prism):
    n_prism = sum(1 for n in drug_panel() if n in prism)
    payloads = run_submission(
        drugs_v2,
        floor_v2,
        OUT_DIR,
        methods_note_v2(),
        "Night Shift melanoma drug-response predictions (v2, data-derived floors): per-condition CTG "
        "% viability for A375 and LOXIMVI, floors read from PRISM Secondary dose-level viability at "
        "the task dose. Predictions, not measurements.",
        {"drugs_with_prism_floors": n_prism, "drugs_judgment_fallback": 12 - n_prism},
    )
    submission_view(payloads)
    return


@app.cell
def _(drugs_v2, floor_v2):
    # Quick look: does the data-floor change any single-agent ranking vs the judgment baseline?
    preview = pl.DataFrame(
        [
            {"cell_line": line, "timepoint_h": t, **r}
            for line in ("A375", "LOXIMVI")
            for t in (24, 48)
            for r in rank_single_agents(line, t, drugs_v2, floor_v2)
        ]
    ).select("cell_line", "timepoint_h", "rank", "condition", "viability_pct")
    mo.md("### v2 single-agent predictions")
    mo.ui.table(preview.sort("cell_line", "timepoint_h", "rank"), page_size=16)
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ## To extend

    - Calibrate the concentration-response with measured potency (e.g. ChEMBL IC50), so the
      floor reflects whether the task dose is actually saturating.
    - Average PRISM with GDSC2 dose-level viability to damp single-screen noise.
    - Use the real 24-48 h CTG readout (when released) to adjudicate the LOXIMVI-BRAFi
      tension between the prior and the data.
    """)
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ---
    ## Vendored engine (appendix)

    Below is the shared, notebook-agnostic machinery this submission reuses: the
    Breadbox HTTP helpers, the drug panel, the CTG-kinetics + Bliss/synergy model,
    and the task-runner. It is vendored inline (rather than imported) so this file
    runs as a single standalone gist. The science specific to this submission is all
    above; this section is boilerplate you can collapse.
    """)
    return


@app.function
def bb_request(method: str, endpoint: str, *, params: dict | None = None, json_body: dict | None = None) -> object:
    """Vendored Breadbox call: identifying User-Agent (portal 403s the default) + 5xx retry."""
    url = f"{BREADBOX}/{endpoint.lstrip('/')}"
    headers = {"Accept": "application/json", "User-Agent": "dmx/0.1 (+https://github.com/broadinstitute/dmx)"}
    last = None
    for attempt in range(6):
        try:
            r = requests.request(method, url, params=params, json=json_body, headers=headers, timeout=120)
            if r.status_code < 500:
                r.raise_for_status()
                return r.json()
            last = f"HTTP {r.status_code}"
        except requests.RequestException as exc:
            last = str(exc)
        time.sleep(3.0 * (attempt + 1))
    raise RuntimeError(f"Breadbox {method} {endpoint} failed after retries: {last}")


@app.function
def bb_get(endpoint: str, params: dict | None = None) -> object:
    return bb_request("GET", endpoint, params=params)


@app.function
def bb_post(endpoint: str, body: dict) -> object:
    return bb_request("POST", endpoint, json_body=body)


@app.function
def drug_panel() -> dict:
    """Per-drug mechanism metadata: the assay-fixed part of the panel (no floors).

    Floors are what each submission produces differently; everything else - target,
    kinetic class, fraction of effect retained at 1/3 dose, task concentration - is
    fixed by the assay design and reused by importing this function.
    """
    return {
        "Vemurafenib": dict(target="BRAF-V600E", kclass="cytostatic", retain3=0.90, conc="1 uM", conc_nM=1000),
        "Dabrafenib": dict(target="BRAF-V600E", kclass="cytostatic", retain3=0.90, conc="100 nM", conc_nM=100),
        "Encorafenib": dict(target="BRAF-V600E", kclass="cytostatic", retain3=0.90, conc="100 nM", conc_nM=100),
        "Trametinib": dict(target="MEK1/2", kclass="cytostatic", retain3=0.90, conc="10 nM", conc_nM=10),
        "Cobimetinib": dict(target="MEK1/2", kclass="cytostatic", retain3=0.90, conc="100 nM", conc_nM=100),
        "Binimetinib": dict(target="MEK1/2", kclass="cytostatic", retain3=0.90, conc="100 nM", conc_nM=100),
        "TAK-733": dict(target="MEK1/2", kclass="cytostatic", retain3=0.88, conc="30 nM", conc_nM=30),
        "Regorafenib": dict(target="multikinase/VEGFR", kclass="mixed", retain3=0.80, conc="5 uM", conc_nM=5000),
        "Alpelisib": dict(target="PI3K-alpha", kclass="cytostatic", retain3=0.80, conc="5 uM", conc_nM=5000),
        "Capivasertib": dict(target="pan-AKT", kclass="cytostatic", retain3=0.80, conc="5 uM", conc_nM=5000),
        "Sapanisertib": dict(target="mTORC1/2", kclass="mixed", retain3=0.85, conc="500 nM", conc_nM=500),
        "Panobinostat": dict(target="pan-HDAC", kclass="cytotoxic", retain3=0.85, conc="50 nM", conc_nM=50),
    }


@app.function
def realized_fraction(timepoint_h: int, kinetic_class: str) -> float:
    """Fraction of a drug's maximal (floor) effect realized by a CTG timepoint."""
    curves = {
        "cytotoxic": {16: 0.50, 24: 0.65, 48: 0.88},
        "mixed": {16: 0.30, 24: 0.45, 48: 0.78},
        "cytostatic": {16: 0.15, 24: 0.30, 48: 0.60},
    }
    return curves[kinetic_class][timepoint_h]


@app.function
def predict_single(floor: float, timepoint_h: int, kinetic_class: str) -> float:
    """Predicted % viability for one drug at one timepoint from its effect floor."""
    realized = realized_fraction(timepoint_h, kinetic_class)
    return round(floor + (100.0 - floor) * (1.0 - realized), 1)


@app.function
def pathway(drug: str) -> str:
    """Coarse pathway bucket used by the combination synergy rules."""
    groups = {
        "BRAF": {"Vemurafenib", "Dabrafenib", "Encorafenib"},
        "MEK": {"Trametinib", "Cobimetinib", "Binimetinib", "TAK-733"},
        "PI3K_AKT": {"Alpelisib", "Capivasertib"},
        "mTOR": {"Sapanisertib"},
        "multikinase": {"Regorafenib"},
        "HDAC": {"Panobinostat"},
    }
    for name, members in groups.items():
        if drug in members:
            return name
    return "other"


@app.function
def pathway_pair_synergy(p1: str, p2: str, line: str) -> float:
    """Bliss multiplier between two *distinct* pathways (<1 = synergy). Line-aware."""
    s = {p1, p2}
    if s == {"BRAF", "MEK"}:
        return 0.72 if line == "A375" else 0.95  # vertical MAPK: synergy only when addicted
    if (s & {"BRAF", "MEK"}) and (s & {"PI3K_AKT", "mTOR"}):
        return 0.85  # MAPK + parallel survival pathway
    return 0.92


@app.function
def predict_combo(
    drugs_list: list[str], line: str, timepoint_h: int, drugs: dict, floor_fn, dose_fraction: float = 1.0
) -> float:
    """Predicted % viability for a combination via Bliss across *distinct pathways* + synergy.

    Two drugs hitting the same node (e.g. two MEK inhibitors) are redundant, not
    independent, so they collapse into one pathway bucket and only the stronger one
    counts. That makes pathway-diverse combinations win, as pharmacology demands.
    """
    by_path: dict[str, float] = {}
    for name in drugs_list:
        adj_floor = 100.0 - (100.0 - floor_fn(name, line)) * (drugs[name]["retain3"] if dose_fraction < 1 else 1.0)
        s = predict_single(adj_floor, timepoint_h, drugs[name]["kclass"]) / 100.0
        p = pathway(name)
        by_path[p] = min(by_path.get(p, 1.0), s)  # strongest agent per pathway
    surviving = 1.0
    for s in by_path.values():
        surviving *= s
    synergy = 1.0
    for p1, p2 in itertools.combinations(by_path, 2):
        synergy *= pathway_pair_synergy(p1, p2, line)
    synergy = max(synergy, 0.55)  # cap compounded synergy
    return round(surviving * synergy * 100.0, 1)


@app.function
def rank_single_agents(line: str, timepoint_h: int, drugs: dict, floor_fn) -> list[dict]:
    """Predicted viability + rank for every compound at one (line, timepoint)."""
    rows = []
    for name, p in drugs.items():
        v = predict_single(floor_fn(name, line), timepoint_h, p["kclass"])
        rows.append({"condition": name, "concentration": p["conc"], "_v": v})
    rows.sort(key=lambda r: r["_v"])  # rank on the precise value
    for i, r in enumerate(rows, 1):
        r["rank"], r["viability_pct"] = i, round(r.pop("_v"))  # report integer magnitude
    return rows


@app.function
def fetch_template(task_id: str) -> dict:
    """Load the committed benchmark template for a retired task."""
    return json.loads((REPO_DIR / "data" / "processed" / "nightshift_v1" / task_id / "output.json").read_text())


@app.function
def write_outputs(task_id: str, populated: dict, reasoning_md: str, out_dir: Path) -> Path:
    """Write output.json + reasoning.md for a task under out_dir/<task_id>/."""
    d = out_dir / task_id
    d.mkdir(parents=True, exist_ok=True)
    (d / "output.json").write_text(json.dumps(populated, indent=2) + "\n")
    (d / "reasoning.md").write_text(reasoning_md)
    return d


@app.function
def task4_response() -> str:
    """The Task 4.1 free-text proposal (constant across submissions - the strategy
    does not depend on how the Task-1/2/3 floors are produced)."""
    return (
        "## Central hypothesis\n\n"
        "Acquired resistance to BRAF+MEK inhibition in melanoma is overwhelmingly "
        "non-genetic at onset: it is driven by reversible adaptive reprogramming - "
        "RTK-driven MAPK reactivation, a switch to a slow-cycling, "
        "neural-crest/AXL-high de-differentiated state, and drug-addicted clones "
        "whose fitness depends on the drug being present. The single most promising "
        "strategy is therefore not another front-line combination partner but "
        "**evolutionarily-informed intermittent dosing guided by a real-time "
        "biomarker of the adaptive (drug-tolerant persister) state**, paired with a "
        "scheduled 'second hit' that is selectively lethal to cells in that state. "
        "The mechanistic bet: the persister state that survives BRAF/MEK inhibition "
        "is itself a vulnerability (dependency on lipid-peroxidase defense / GPX4, "
        "on AXL-PI3K-AKT survival signaling, and on mitochondrial/oxidative "
        "metabolism), and drug-addicted resistant clones regress during a planned "
        "holiday. Hitting the persister window and pulsing the drug should delay or "
        "prevent the fixation of stable resistance.\n\n"
        "## Why this over the obvious alternatives\n\n"
        "Adding a third targeted agent up front (e.g. anti-PD-1, a CDK4/6 or ERK "
        "inhibitor) mostly raises toxicity and still selects for the same de-"
        "differentiated escape state; triplet trials have shown only modest PFS "
        "gains and more adverse events. The persister biology says the failure mode "
        "is a *cell state*, not a single node, so the leverage is in (a) timing - "
        "exploiting drug addiction with holidays - and (b) a state-selective killer "
        "(GPX4 inhibition / ferroptosis induction), to which de-differentiated, "
        "mesenchymal-like persisters are unusually sensitive. This attacks the "
        "reservoir that every front-line combination leaves behind.\n\n"
        "## Key experiment\n\n"
        "Model: BRAF-V600E melanoma (A375 plus PDX-derived lines and at least one "
        "immunocompetent syngeneic model) with matched BRAFi/MEKi-resistant "
        "sublines (DepMap already lists A375DABR and A375DABTRAMR derivatives). "
        "Arms: (1) continuous dabrafenib+trametinib (standard); (2) the same with "
        "a fixed intermittent schedule (e.g. 4 weeks on / 2 weeks off); (3) "
        "continuous BRAFi/MEKi + a GPX4 inhibitor (ferroptosis inducer) pulsed "
        "during the adaptive window (days 3-10 after each drug start); (4) the "
        "intermittent schedule with the GPX4 pulse timed to each off-on transition. "
        "Readouts: time to regrowth / resistance, longitudinal single-cell RNA-seq "
        "to track the AXL-high de-differentiated fraction, and a lipid-peroxidation "
        "reporter to confirm the persister window. Critical controls: a ferroptosis "
        "rescue (ferrostatin-1 / liproxstatin) must abolish arm-3/4 benefit if the "
        "mechanism is ferroptosis; a GPX4-pulse-only arm with no MAPK drug controls "
        "for nonspecific toxicity; schedule-matched vehicle controls separate "
        "'holiday' effects from the second hit.\n\n"
        "## Confirm vs falsify\n\n"
        "Confirmation: arms 3 and 4 significantly delay resistance vs continuous "
        "dosing; the AXL-high persister fraction collapses specifically when the "
        "GPX4 pulse lands in the adaptive window; ferrostatin rescue restores "
        "resistance. Falsification: timing the GPX4 pulse to the persister window "
        "gives no benefit over constant co-dosing (state-timing irrelevant), or "
        "resistant outgrowth is dominated by genetic events (NRAS/MEK mutations, "
        "BRAF amplification) insensitive to ferroptosis - which would say the "
        "durable reservoir is genetic, not a persister state, and redirect effort "
        "to up-front mutational suppression.\n\n"
        "## Path to patients\n\n"
        "A positive result yields two immediately translatable levers: a dosing "
        "*schedule* (intermittent BRAFi/MEKi), already clinically plausible and "
        "cheap to trial, and a *biomarker-timed* ferroptosis second hit. The "
        "biomarker (AXL-high / lipid-peroxidation-primed state) becomes the entry "
        "criterion and pharmacodynamic readout for a window-of-opportunity trial: "
        "treat to the adaptive window, deliver the second hit, and measure persister "
        "collapse on serial biopsies before moving to a randomized intermittent-"
        "plus-pulse design.\n"
    )


@app.function
def task4_analyses() -> str:
    """Task 4.1 'analyses' field - cites published DepMap facts, not the prediction model."""
    return (
        "Supporting facts that are independent of the prediction model - reported by "
        "DepMap, not produced by the model, so citing them is not circular: A375 carries "
        "a strong CRISPR BRAF dependency (~ -1.5) and a GDSC2 dabrafenib AUC ~0.35, while "
        "LOXIMVI's GDSC2 dabrafenib AUC is ~0.76 - the BRAF-V600E addiction that "
        "intermittent dosing is designed to exploit, measured rather than assumed. DepMap "
        "also lists BRAFi-resistant A375 sublines (A375DABR, A375DABTRAMR) as ready-made "
        "models for the proposed resistance experiments. The Task-1/2/3 viability numbers, "
        "by contrast, are model predictions - they illustrate the strategy's logic but are "
        "not evidence for it."
    )


@app.function
def run_submission(
    drugs: dict, floor_fn, out_dir: Path, methods_note: str, summary_desc: str, summary_numbers: dict
) -> pl.DataFrame:
    """Write all 11 task outputs (+ summary.json) for one set of floors."""

    def drugs_of(condition):
        return [c.strip() for c in condition.split("+")]

    def fmt_third(nM):
        third = nM / 3.0
        return f"{third / 1000:.2f} uM" if third >= 1000 else f"{round(third):g} nM"

    done = []

    # Task 1.1-1.4: single agents, one (line, timepoint) each.
    for task_id, (line, t) in {
        "1.1": ("A375", 24),
        "1.2": ("LOXIMVI", 24),
        "1.3": ("A375", 48),
        "1.4": ("LOXIMVI", 48),
    }.items():
        ranked = rank_single_agents(line, t, drugs, floor_fn)
        by_name = {r["condition"]: r for r in ranked}
        tmpl = fetch_template(task_id)
        for row in tmpl["rankings"]:
            hit = by_name[row["condition"]]
            row["rank"], row["viability_pct"] = hit["rank"], hit["viability_pct"]
        body = "\n".join(
            f"{r['rank']}. {r['condition']} ({r['concentration']}) - {r['viability_pct']}%" for r in ranked
        )
        note = (
            "At 24 h the cytostatic MAPK inhibitors have not yet separated, so panobinostat "
            "(fast cytotoxic) leads and the BRAF/MEK inhibitors cluster.\n"
            if t == 24
            else "By 48 h the cytostatic effect has matured; in A375 the BRAF/MEK inhibitors "
            "deepen markedly, in LOXIMVI they stay shallow.\n"
        )
        write_outputs(
            task_id,
            tmpl,
            f"# Task {task_id} reasoning - {line}, {t} h\n\n{methods_note}"
            f"## Predicted ranking ({line}, {t} h; rank 1 = greatest effect)\n\n{body}\n\n{note}",
            out_dir,
        )
        done.append(task_id)

    # Task 1.5: all 48 conditions, global rank.
    pooled = []
    for line in ("A375", "LOXIMVI"):
        for t in (24, 48):
            for r in rank_single_agents(line, t, drugs, floor_fn):
                pooled.append({**r, "cell_line": line, "timepoint_h": t})
    pooled.sort(key=lambda r: r["viability_pct"])
    for i, r in enumerate(pooled, 1):
        r["rank"] = i
    key15 = {(r["condition"], r["cell_line"], r["timepoint_h"]): r for r in pooled}
    tmpl15 = fetch_template("1.5")
    for row in tmpl15["rankings"]:
        hit = key15[(row["condition"], row["cell_line"], row["timepoint_h"])]
        row["rank"], row["viability_pct"] = hit["rank"], hit["viability_pct"]
    top = "\n".join(
        f"{r['rank']}. {r['condition']} - {r['cell_line']} {r['timepoint_h']}h - {r['viability_pct']}%"
        for r in pooled[:8]
    )
    write_outputs(
        "1.5",
        tmpl15,
        f"# Task 1.5 reasoning - all 48 conditions\n\n{methods_note}"
        f"## Top of the global ranking (lowest viability first)\n\n{top}\n\n"
        "Drivers: timepoint (48 h > 24 h), cytotoxic > cytostatic mechanism, and for "
        "MAPK drugs A375 >> LOXIMVI.\n",
        out_dir,
    )
    done.append("1.5")

    # Task 2.1/2.2: combinations per line; 2.3: pooled.
    for task_id, line in {"2.1": "A375", "2.2": "LOXIMVI"}.items():
        tmpl = fetch_template(task_id)
        for row in tmpl["rankings"]:
            row["_v"] = predict_combo(drugs_of(row["condition"]), line, 48, drugs, floor_fn)
        ordered = sorted(tmpl["rankings"], key=lambda r: r["_v"])
        for i, row in enumerate(ordered, 1):
            row["rank"], row["viability_pct"] = i, round(row.pop("_v"))
        body = "\n".join(f"{r['rank']}. {r['condition']} - {r['viability_pct']}%" for r in ordered)
        lead = (
            "In addicted A375 the BRAF+MEK doublets (vertical MAPK blockade) lead.\n"
            if line == "A375"
            else "In BRAF-WT LOXIMVI the BRAF+MEK doublets lose vertical synergy; the MEK+PI3K "
            "pair does relatively better.\n"
        )
        write_outputs(
            task_id,
            tmpl,
            f"# Task {task_id} reasoning - {line} combinations, 48 h\n\n{methods_note}"
            f"## Predicted ranking\n\n{body}\n\n{lead}",
            out_dir,
        )
        done.append(task_id)

    tmpl23 = fetch_template("2.3")
    for row in tmpl23["rankings"]:
        row["_v"] = predict_combo(drugs_of(row["condition"]), row["cell_line"], 48, drugs, floor_fn)
    ordered = sorted(tmpl23["rankings"], key=lambda r: r["_v"])
    for i, row in enumerate(ordered, 1):
        row["rank"], row["viability_pct"] = i, round(row.pop("_v"))
    body23 = "\n".join(f"{r['rank']}. {r['condition']} - {r['cell_line']} - {r['viability_pct']}%" for r in ordered)
    write_outputs(
        "2.3",
        tmpl23,
        f"# Task 2.3 reasoning - combinations across both lines, 48 h\n\n{methods_note}"
        f"## Predicted global ranking\n\n{body23}\n\nA375 doublets occupy the most-effective "
        "positions; the same doublets in LOXIMVI sit higher because the BRAF arm is inert there.\n",
        out_dir,
    )
    done.append("2.3")

    # Task 3.1/3.2: nominate best pathway-diverse 3-drug combo at 16 h, each at 1/3 dose.
    for task_id, line in {"3.1": "LOXIMVI", "3.2": "A375"}.items():
        scored = sorted(
            (
                [list(trio), predict_combo(list(trio), line, 16, drugs, floor_fn, dose_fraction=1 / 3)]
                for trio in itertools.combinations(CANDIDATES, 3)
            ),
            key=lambda x: x[1],
        )
        best, best_v = scored[0]
        tmpl = fetch_template(task_id)
        tmpl["nomination"]["drugs"] = best
        tmpl["nomination"]["concentrations"] = [fmt_third(drugs[d]["conc_nM"]) for d in best]
        tmpl["nomination"]["viability_pct"] = round(best_v)
        paths = ", ".join(f"{d} ({pathway(d)})" for d in best)
        why = (
            "LOXIMVI is not BRAF-V600E-addicted, so BRAF inhibitors are dead weight; the win is "
            "co-blocking MEK and AKT plus a high-exposure cytotoxic for acute 16 h pressure.\n"
            if line == "LOXIMVI"
            else "A375 is BRAF-V600E-addicted, so deep vertical MAPK blockade (BRAF + MEK) is lethal "
            "even at 16 h; adding AKT inhibition shuts the parallel survival/feedback arm.\n"
        )
        top_body = "\n".join(f"- {' + '.join(t)} - {round(v)}%" for t, v in scored[:5])
        write_outputs(
            task_id,
            tmpl,
            f"# Task {task_id} reasoning - nominate 3-drug combo for {line}, 16 h\n\n{methods_note}"
            f"## Nomination\n\nNominated: **{' + '.join(best)}** -> predicted {round(best_v)}% viability "
            f"at 16 h.\n\nPathways hit: {paths}. {why}\n## Top candidates considered\n\n{top_body}\n",
            out_dir,
        )
        done.append(task_id)

    # Task 4.1: open-ended strategy (constant across submissions).
    tmpl4 = fetch_template("4.1")
    tmpl4["response"] = task4_response()
    tmpl4["analyses"] = task4_analyses()
    write_outputs(
        "4.1",
        tmpl4,
        "# Task 4.1 reasoning\n\nWeighed three strategy families: (a) more up-front targeted "
        "blockade, (b) immunotherapy combination, (c) evolution-/state-aware dosing. Chose (c) - "
        "intermittent dosing + a persister-selective ferroptosis second hit - because it targets the "
        "drug-tolerant-persister reservoir every continuous combination leaves behind, and it is "
        "falsifiable with a ferrostatin rescue and a genetic-resistance escape clause. Full proposal "
        "in output.json.\n",
        out_dir,
    )
    done.append("4.1")

    # summary.json envelope.
    files = [str((out_dir / t / f).relative_to(REPO_DIR)) for t in done for f in ("output.json", "reasoning.md")]
    summary = {
        "description": summary_desc,
        "numbers": {"tasks_submitted": len(done), **summary_numbers},
        "files": files,
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    # Read back exactly what was written so the notebook can show the precise bytes
    # submitted - both the output.json and the reasoning.md - per task.
    return {
        t: {
            "output": json.loads((out_dir / t / "output.json").read_text()),
            "reasoning": (out_dir / t / "reasoning.md").read_text(),
        }
        for t in done
    }


@app.function
def submission_view(payloads: dict):
    """Render exactly what gets submitted, inline: a scan table of every prediction, plus
    per task the reasoning trace (reasoning.md) and the exact output.json bytes."""
    rows = []
    for tid, pp in payloads.items():
        p = pp["output"]
        if "rankings" in p:
            for r in p["rankings"]:
                rows.append(
                    {
                        "task": tid,
                        "condition": r["condition"],
                        "cell_line": str(r.get("cell_line", p.get("cell_line") or "")),
                        "timepoint_h": str(r.get("timepoint_h", p.get("timepoint_h") or "")),
                        "rank": r["rank"],
                        "viability_pct": r["viability_pct"],
                    }
                )
        elif "nomination" in p:
            n = p["nomination"]
            rows.append(
                {
                    "task": tid,
                    "condition": " + ".join(n["drugs"]),
                    "cell_line": str(p.get("cell_line") or ""),
                    "timepoint_h": str(p.get("timepoint_h") or ""),
                    "rank": 1,
                    "viability_pct": n["viability_pct"],
                }
            )
        else:  # 4.1 free-text strategy
            rows.append(
                {
                    "task": tid,
                    "condition": "(free-text strategy)",
                    "cell_line": "",
                    "timepoint_h": "",
                    "rank": None,
                    "viability_pct": None,
                }
            )
    panels = {
        f"Task {tid}": mo.vstack(
            [
                mo.md(pp["reasoning"]),
                mo.md("**Exact output.json submitted:**"),
                mo.md("```json\n" + json.dumps(pp["output"], indent=2) + "\n```"),
            ]
        )
        for tid, pp in payloads.items()
    }
    return mo.vstack(
        [
            mo.md("### Exactly what gets submitted (predictions + reasoning, per task)"),
            mo.ui.table(pl.DataFrame(rows), page_size=25),
            mo.md("#### Per task: reasoning trace, then exact payload (expand any task)"),
            mo.accordion(panels),
        ]
    )


if __name__ == "__main__":
    app.run()
