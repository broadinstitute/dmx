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
    import time
    from pathlib import Path

    import altair as alt
    import marimo as mo
    import polars as pl
    import requests

    # Self-contained: this notebook vendors its own Breadbox HTTP helpers (below) so it
    # runs as a standalone gist with no sibling files. REPO_DIR is just where outputs land.
    NOTEBOOK_DIR = Path(__file__).resolve().parent
    REPO_DIR = NOTEBOOK_DIR.parent if NOTEBOOK_DIR.name == "notebooks" else NOTEBOOK_DIR
    OUT_DIR = REPO_DIR / "data" / "processed" / "nightshift_v1"
    BREADBOX = "https://depmap.org/portal/breadbox"

    A375 = "ACH-000219"
    LOXIMVI = "ACH-000750"

    KARMAN = "https://mcp.karmanai.org"
    PRISM_SEC_AUC = "07b7bda9-ae00-43b3-bca1-336b9607f8f5"
    GDSC2_AUC = "2eac8e7b-beb4-48c1-b78f-c226723e54d7"

    # Task-3 nomination candidate list (panobinostat + sapanisertib excluded by the task).
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


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    # Night Shift submission - reasoning baseline

    Predict CellTiter-Glo (CTG) % viability for two melanoma lines (A375, LOXIMVI)
    under a MAPK-inhibitor panel at 16 / 24 / 48 h, and rank conditions by effect.
    No experimental readout is provided, so these are predictions.

    The predictions turn on two hypotheses, and the cells below are the loop that tests
    them - state a hypothesis, run a pull, look at what came back - not a pile of tables:

    - **H1 (cell-line identity).** A375 is BRAF-V600E-addicted; LOXIMVI is not
      BRAF-V600E-driven. Tested against DepMap CRISPR dependency and drug sensitivity in
      Instruments 1-2. The genetics pull settles A375 cleanly but leaves LOXIMVI
      *unsettled* - which is exactly the seam a later data source can pull apart.
    - **H2 (assay kinetics).** At a 16-48 h CTG endpoint, cytostatic MAPK inhibitors barely
      move viability while broadly cytotoxic agents (HDAC, mTOR) drop it fast, so **the
      lowest-viability drug at 24 h is usually not the most potent one.** This is a
      mechanistic prior the model encodes; with no readout it is *not* tested here, and its
      falsifier is named at the end.

    The terminal hypothesis - what gets committed - is the ranked, predicted viability. Each
    drug gets a per-line **effect floor** set by judgment (anchored to the DepMap data, but a
    judgment call) and a kinetic class, run through
    `floor + (100-floor) * (1 - realized_fraction(t, class))`. These rankings are model
    outputs, not independent findings; the "evidence vs assumption" cell at the end marks
    which parts are load-bearing data and which are my knobs. The mechanics are vendored in
    the appendix so this file runs as a standalone gist; it is one of a small set of
    submissions that swap how the floors are produced (see the repo README).
    """)
    return


@app.function
def methods_note_v1() -> str:
    """v1 methods preamble embedded in each reasoning.md."""
    return (
        "## Method (v1 - reasoning baseline)\n\n"
        "No experimental readout was provided; these are predictions. Each drug has a per-line "
        "**effect floor** (asymptotic % viability at its task dose, set by *judgment* anchored to "
        "PRISM/GDSC AUC + dose-vs-IC50 margin + pathway dependence) and a **kinetic class**. "
        "Prediction = `floor + (100-floor) * (1 - realized_fraction(t, class))`. Combinations use "
        "Bliss across distinct pathways with a line-aware synergy factor.\n\n"
    )


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ## Instrument 1 - what are these two lines addicted to?

    CRISPR (Chronos) dependency and hotspot calls from DepMap. Strongly negative
    dependency means the line dies without that gene.
    """)
    return


@app.cell
def _():
    def _line_values(gene, dataset):
        df = matrix_feature(dataset, gene, value_name="v")
        m = dict(zip(df["depmap_id"].to_list(), df["v"].to_list())) if not df.is_empty() else {}
        return {"A375": m.get(A375), "LOXIMVI": m.get(LOXIMVI)}

    genetics = pl.DataFrame(
        [
            {"feature": "BRAF dependency (Chronos)", **_line_values("BRAF", "Chronos_Combined")},
            {"feature": "MAP2K1 (MEK1) dependency", **_line_values("MAP2K1", "Chronos_Combined")},
            {"feature": "BRAF hotspot mutation", **_line_values("BRAF", "mutations_hotspot")},
        ]
    )
    mo.md("### Dependency & mutation: A375 vs LOXIMVI")
    mo.ui.table(genetics, page_size=10)
    return (genetics,)


@app.cell
def _(genetics):
    def _val(feature, line):
        return genetics.filter(pl.col("feature") == feature)[line][0]

    mo.md(
        "**H1, tested against CRISPR.** A375 BRAF dependency = "
        f"`{_val('BRAF dependency (Chronos)', 'A375')}`, MAP2K1 = "
        f"`{_val('MAP2K1 (MEK1) dependency', 'A375')}` - strongly negative, so H1's A375 half is "
        "confirmed: it dies without BRAF/MEK, a BRAF-V600E-addicted line. H1's LOXIMVI half is "
        f"*not* settled here: LOXIMVI BRAF dependency = `{_val('BRAF dependency (Chronos)', 'LOXIMVI')}` "
        "- it is absent from the CRISPR screen, so genetics is silent on it (silence is not "
        "confirmation). Whether LOXIMVI is BRAF-inhibitor-resistant rests entirely on the drug "
        "data below - and that is the exact seam where a different screen can later disagree."
    )
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ## Instrument 2 - empirical drug sensitivity for the actual panel

    PRISM (5-day) and GDSC2 are longer assays than our 16-48 h endpoint and at their
    own doses, so they anchor *relative* efficacy and the line-to-line contrast, not
    the % we predict. **AUC: lower = more sensitive** (1.0 = no effect).
    """)
    return


@app.cell
def _():
    panel_labels = [
        "VEMURAFENIB",
        "DABRAFENIB",
        "ENCORAFENIB",
        "TRAMETINIB",
        "COBIMETINIB",
        "BINIMETINIB",
        "REGORAFENIB",
        "ALPELISIB",
        "CAPIVASERTIB",
        "PANOBINOSTAT",
    ]

    def _auc_frame(dataset_id, label):
        raw = bb_post(f"datasets/matrix/{dataset_id}", {"features": panel_labels, "feature_identifier": "label"})
        rows = []
        for drug in panel_labels:
            v = raw.get(drug, {}) if isinstance(raw, dict) else {}
            rows.append({"drug": drug.title(), f"{label}_A375": v.get(A375), f"{label}_LOXIMVI": v.get(LOXIMVI)})
        return pl.DataFrame(rows)

    sensitivity = _auc_frame(PRISM_SEC_AUC, "PRISM").join(_auc_frame(GDSC2_AUC, "GDSC2"), on="drug", how="left")
    mo.md("### Drug-sensitivity AUC anchors (lower = more sensitive)")
    mo.ui.table(sensitivity.with_columns(pl.col(pl.Float64).round(3)), page_size=12)
    return (sensitivity,)


@app.cell
def _(sensitivity):
    def _auc(drug, col):
        v = sensitivity.filter(pl.col("drug") == drug)[col][0]
        return "n/a" if v is None else round(v, 2)

    mo.md(
        "**H1's LOXIMVI half, tested against drug response.** GDSC2 dabrafenib: A375 "
        f"`{_auc('Dabrafenib', 'GDSC2_A375')}` (sensitive) vs LOXIMVI "
        f"`{_auc('Dabrafenib', 'GDSC2_LOXIMVI')}` (resistant) - this longer assay settles LOXIMVI as "
        "BRAF-inhibitor-resistant, so H1 stands on this evidence. (Hold that lightly: a different "
        "screen disagrees, and a companion submission leads with that collision.) Consistent with H2, "
        f"PRISM panobinostat is the single most potent agent in both lines (A375 `{_auc('Panobinostat', 'PRISM_A375')}` "
        f"/ LOXIMVI `{_auc('Panobinostat', 'PRISM_LOXIMVI')}`) - a broadly cytotoxic HDAC inhibitor, not a "
        "MAPK drug. These live numbers anchor the floors below."
    )
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ## Effect floors (judgment, data-anchored)

    The floors below are curated estimates from the Instrument-2 AUCs + dose-vs-IC50
    margin + pathway dependence. BRAF inhibitors: deep in A375, shallow in LOXIMVI.
    PI3K/AKT: weak everywhere. Panobinostat: deep in both (cytotoxic). These judgment
    floors are this submission's one soft spot - the obvious next step is to read them
    straight from data instead.
    """)
    return


@app.cell
def _():
    # v1 judgment floors (asymptotic % viability at task dose), per line.
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
    drugs_v1 = {n: {**meta, "floorA": FLOORS_V1[n][0], "floorL": FLOORS_V1[n][1]} for n, meta in drug_panel().items()}

    def floor_v1(name, line):
        return drugs_v1[name]["floorA" if line == "A375" else "floorL"]

    drug_table = pl.DataFrame(
        [
            {
                "drug": n,
                "target": p["target"],
                "class": p["kclass"],
                "floor_A375": p["floorA"],
                "floor_LOXIMVI": p["floorL"],
                "conc": p["conc"],
            }
            for n, p in drugs_v1.items()
        ]
    )
    mo.md("### v1 drug parameter table")
    mo.ui.table(drug_table, page_size=12)
    return drugs_v1, floor_v1


@app.cell
def _(drugs_v1, floor_v1):
    single_preview = pl.DataFrame(
        [
            {"cell_line": line, "timepoint_h": t, **r}
            for line in ("A375", "LOXIMVI")
            for t in (24, 48)
            for r in rank_single_agents(line, t, drugs_v1, floor_v1)
        ]
    ).select("cell_line", "timepoint_h", "rank", "condition", "concentration", "viability_pct")
    mo.md("### All single-agent predictions (A375/LOXIMVI x 24/48 h)")
    mo.ui.table(single_preview.sort("cell_line", "timepoint_h", "rank"), page_size=16)
    return (single_preview,)


@app.cell
def _(single_preview):
    _chart = (
        alt.Chart(
            single_preview.with_columns(
                (pl.col("cell_line") + " " + pl.col("timepoint_h").cast(pl.Utf8) + "h").alias("cond")
            )
        )
        .mark_circle(size=70, opacity=0.85)
        .encode(
            x=alt.X("viability_pct:Q", title="predicted % viability", scale=alt.Scale(domain=[0, 105])),
            y=alt.Y("condition:N", sort=alt.EncodingSortField("viability_pct", op="min"), title=None),
            color=alt.Color("cond:N", title=None),
            tooltip=["condition", "cell_line", "timepoint_h", "viability_pct"],
        )
        .properties(width=560, height=420, title="Lower = more effect. MAPK drugs bite A375, spare LOXIMVI.")
    )
    mo.ui.altair_chart(_chart)
    return


@app.cell
def _(drugs_v1, floor_v1, genetics, sensitivity):
    def _g(feature, line):
        return genetics.filter(pl.col("feature") == feature)[line][0]

    def _s(drug, col):
        return sensitivity.filter(pl.col("drug") == drug)[col][0]

    status = run_submission(
        drugs_v1,
        floor_v1,
        OUT_DIR,
        methods_note_v1(),
        "Night Shift melanoma drug-response predictions (v1, reasoning baseline): per-condition CTG "
        "% viability for A375 (BRAF-V600E) and LOXIMVI from judgment effect floors + a CTG-kinetics "
        "model. Predictions, not measurements.",
        {
            "A375_BRAF_dependency_chronos": _g("BRAF dependency (Chronos)", "A375"),
            "GDSC2_dabrafenib_AUC_A375": _s("Dabrafenib", "GDSC2_A375"),
            "GDSC2_dabrafenib_AUC_LOXIMVI": _s("Dabrafenib", "GDSC2_LOXIMVI"),
        },
    )
    mo.md("### v1 submission complete - 11 tasks + summary.json under data/processed/nightshift_v1/")
    mo.ui.table(status, page_size=12)
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ## What is evidence, what is assumption

    - **Independent evidence (DepMap, live above):** A375 BRAF/MAP2K1 dependency; PRISM/GDSC AUCs;
      the GDSC2 dabrafenib A375-vs-LOXIMVI contrast. These could have contradicted the prior and did not.
    - **Assumptions (mine):** the per-line floors (judgment), the three kinetic curves, the synergy
      constants. "BRAF+MEK doublets win in A375" is the `0.72 vs 0.95` synergy input read back out.

    **Calibration.** High confidence in the *direction* of the big contrasts (MAPK inhibitors hit
    A375 far harder than LOXIMVI; panobinostat tops both; clinical doublets are the best A375 combos).
    Low confidence in absolute magnitudes and the fine ordering within the clustered 24 h MAPK
    inhibitors. **Falsifier:** a real 24 h A375 readout where a BRAF/MEK inhibitor drops viability
    below panobinostat would break the kinetics premise.

    ## To extend

    - Read the effect floors straight from dose-level viability data instead of judgment.
    - Calibrate the drug concentration-response with measured potency (e.g. ChEMBL IC50s).
    - Calibrate the kinetic curves against any one real CTG time-course.
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
def matrix_feature(dataset_id: str, feature_label: str, value_name: str = "value") -> pl.DataFrame:
    """One feature column from a Breadbox matrix dataset, as a depmap_id -> value frame."""
    raw = bb_post(f"datasets/matrix/{dataset_id}", {"features": [feature_label], "feature_identifier": "label"})
    values = raw.get(feature_label) if isinstance(raw, dict) else None
    if values is None and isinstance(raw, dict) and raw:
        values = next(iter(raw.values()))
    if not isinstance(values, dict):
        return pl.DataFrame()
    return pl.DataFrame({"depmap_id": list(values.keys()), value_name: list(values.values())})


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
    """GET the benchmark's output template so we populate, never restructure."""
    r = requests.get(f"{KARMAN}/tasks/{task_id}/output.json", headers={"Accept": "application/json"}, timeout=30)
    r.raise_for_status()
    return r.json()


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

    return pl.DataFrame({"task": done}).with_columns(pl.lit("written").alias("status"))


if __name__ == "__main__":
    app.run()
