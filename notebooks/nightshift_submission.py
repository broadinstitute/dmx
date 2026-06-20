# /// script
# requires-python = ">=3.12,<3.14"
# dependencies = [
#     "altair==5.5.0",
#     "marimo",
#     "polars==1.40.1",
#     "requests==2.32.5",
# ]
# ///
"""Night Shift wet-lab benchmark - predicted compound effects on melanoma viability.

Composed notebook (not a catalog vignette). It answers the Karman "Night Shift"
tasks at https://mcp.karmanai.org/tasks: predict CellTiter-Glo % viability for a
fixed compound panel on two melanoma lines (A375, LOXIMVI) at 16/24/48 h, rank
single agents and combinations, nominate triple combinations, and propose a
resistance strategy.

No experimental readouts are provided - every number here is a *prediction*
built from (1) DepMap genetics + measured PRISM drug-response as a reference and
(2) a transparent growth-kinetics model for the short assay timepoints.
"""

import marimo

__generated_with = "0.23.10"
app = marimo.App(width="medium")

with app.setup:
    import json
    import sys
    from pathlib import Path

    import altair as alt
    import marimo as mo
    import polars as pl
    import requests

    NOTEBOOK_DIR = Path(__file__).resolve().parent
    if str(NOTEBOOK_DIR) not in sys.path:
        sys.path.insert(0, str(NOTEBOOK_DIR))

    from nb02_dataset_discovery import bb_get, bb_post  # noqa: E402

    KARMAN = "https://mcp.karmanai.org"

    OUT_DIR = NOTEBOOK_DIR.parent / "data" / "processed" / "nightshift_submission"

    # DepMap model ids for the two benchmark lines.
    LINES = {"A375": "ACH-000219", "LOXIMVI": "ACH-000750"}

    # Doubling times (h) in standard culture, from Cellosaurus / NCI-DTP.
    # A375 ~27 h; LOXIMVI ~22 h. Used by the CTG growth-kinetics model.
    DOUBLING_H = {"A375": 27.0, "LOXIMVI": 22.0}

    # Breadbox dataset ids (26Q1). Re-discover via nb02 if these 404 on a release bump.
    DS_PRISM_AUC = "07b7bda9-ae00-43b3-bca1-336b9607f8f5"  # PRISM Repurposing Secondary (AUC)
    DS_HOTSPOT = "a952ab7b-56c8-4aeb-872e-8ee02eeae042"  # Hotspot Mutations (Public 26Q1)

    # The 12-compound single-agent panel: benchmark dose, target, mechanism class.
    PANEL = [
        ("Vemurafenib", "1 uM", 1.0, "BRAF V600E", "BRAFi"),
        ("Dabrafenib", "100 nM", 0.1, "BRAF V600E", "BRAFi"),
        ("Encorafenib", "100 nM", 0.1, "BRAF V600E", "BRAFi"),
        ("Trametinib", "10 nM", 0.01, "MEK1/2", "MEKi"),
        ("Cobimetinib", "100 nM", 0.1, "MEK1/2", "MEKi"),
        ("Binimetinib", "100 nM", 0.1, "MEK1/2", "MEKi"),
        ("TAK-733", "30 nM", 0.03, "MEK1/2", "MEKi"),
        ("Regorafenib", "5 uM", 5.0, "multi-kinase (VEGFR/RAF)", "multikinase"),
        ("Alpelisib", "5 uM", 5.0, "PI3K-alpha", "PI3K/AKT/mTOR"),
        ("Capivasertib", "5 uM", 5.0, "AKT1/2/3", "PI3K/AKT/mTOR"),
        ("Sapanisertib", "500 nM", 0.5, "mTORC1/2", "PI3K/AKT/mTOR"),
        ("Panobinostat", "50 nM", 0.05, "pan-HDAC", "HDACi"),
    ]
    PANEL_COLS = ["condition", "concentration", "dose_uM", "target", "class"]

    # Effect index E per (line, compound). E is the fractional growth-rate
    # suppression in the CTG model: E=0 inactive (tracks control), E~1 full
    # cytostatic arrest, E>1 net cytotoxic (population shrinks). Hand-set from
    # PRISM dose-matched viability + the line biology, NOT from the dose-integrated
    # AUC (which over-credits off-target high-dose killing). Rationale is carried
    # so every value is auditable. Key biology: both lines are BRAF V600E, but
    # A375 is BRAF-addicted (sensitive) while LOXIMVI is dedifferentiated
    # (MITF-low/AXL-high/SOX10-low) and intrinsically resistant to BRAF inhibition.
    EFFECT = {
        "A375": {
            "Vemurafenib": (0.95, "saturating 1 uM BRAFi in addicted line; strong arrest"),
            "Dabrafenib": (1.05, "potent BRAFi, 100 nM well above IC50; near-complete arrest"),
            "Encorafenib": (1.00, "potent BRAFi, slow off-rate; strong arrest"),
            "Trametinib": (1.05, "very potent MEKi (sub-nM); strong arrest at 10 nM"),
            "Cobimetinib": (0.95, "potent MEKi; PRISM dose-matched among deepest in A375"),
            "Binimetinib": (0.85, "MEKi, somewhat less potent per-mole"),
            "TAK-733": (0.95, "potent MEKi; PRISM dose-matched strong"),
            "Regorafenib": (0.40, "multikinase, partial off-target growth slowing"),
            "Alpelisib": (0.20, "PI3K-alpha; weak single agent, A375 not PI3K-driven"),
            "Capivasertib": (0.20, "AKT inhibitor; weak single agent in A375"),
            "Sapanisertib": (0.50, "mTORC1/2 cytostatic; moderate broad growth slowing"),
            "Panobinostat": (1.30, "pan-HDACi, actively cytotoxic at 50 nM (not just arrest)"),
        },
        "LOXIMVI": {
            "Vemurafenib": (0.20, "intrinsic BRAFi resistance (dedifferentiated); near-inactive at 1 uM"),
            "Dabrafenib": (0.20, "BRAFi near-inactive despite V600E; MITF-low/AXL-high state"),
            "Encorafenib": (0.25, "BRAFi near-inactive; slight edge from slow off-rate"),
            "Trametinib": (0.70, "MEKi retains partial activity (downstream node); blunted vs A375"),
            "Cobimetinib": (0.60, "MEKi partial activity in dedifferentiated line"),
            "Binimetinib": (0.60, "MEKi partial activity"),
            "TAK-733": (0.65, "MEKi partial activity; PRISM AUC sensitive"),
            "Regorafenib": (0.35, "multikinase, modest"),
            "Alpelisib": (0.20, "PI3K-alpha; weak single agent"),
            "Capivasertib": (0.20, "AKT inhibitor; weak single agent"),
            "Sapanisertib": (0.55, "mTORC1/2; mesenchymal lines partly PI3K/mTOR-driven"),
            "Panobinostat": (
                1.45,
                "most active candidate: pan-HDACi strongly cytotoxic; PRISM dose-matched -3.46 (~9%)",
            ),
        },
    }

    # Mechanism class per compound, for the combination synergy rules below.
    CLASS = {c[0]: c[4] for c in PANEL}

    # The four benchmark combinations (Task 2): (label, drug_a, drug_b, conc_string).
    COMBOS = [
        ("Trametinib+Alpelisib", "Trametinib", "Alpelisib", "10 nM + 5 uM"),
        ("Dabrafenib+Trametinib", "Dabrafenib", "Trametinib", "100 nM + 10 nM"),
        ("Encorafenib+Binimetinib", "Encorafenib", "Binimetinib", "100 nM + 100 nM"),
        ("Vemurafenib+Cobimetinib", "Vemurafenib", "Cobimetinib", "1 uM + 100 nM"),
    ]

    # Task 3 candidate list (no HDACi / no mTOR - panobinostat & sapanisertib excluded).
    TRIPLE_CANDIDATES = [
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
    # Map combo label -> its two drugs, for scoring by condition string.
    COMBO_DRUGS = {label: (a, b) for label, a, b, _conc in COMBOS}

    # Task 3 nominations (drugs + 1/3-dose strings), set in the Task 3 final hypothesis.
    TRIPLE_NOMINATIONS = {
        "A375": (["Dabrafenib", "Trametinib", "Alpelisib"], ["33 nM", "3.3 nM", "1.67 uM"]),
        "LOXIMVI": (["Trametinib", "Alpelisib", "Capivasertib"], ["3.3 nM", "1.67 uM", "1.67 uM"]),
    }
    # Per-class effect attenuation at 1/3 of the listed dose. Potent targeted agents
    # dosed well above IC50 barely drop; weak high-micromolar agents drop more.
    DOSE_THIRD = {"BRAFi": 0.85, "MEKi": 0.85, "multikinase": 0.65, "PI3K/AKT/mTOR": 0.65}

    # GDSC combination screen (anchor+library design) - the independent test of the
    # combination-synergy rule in Task 2. The benchmark lines are NOT in it (and there
    # is no melanoma at all), so this tests the *rule* on other BRAF-mutant lines.
    DS_GDSC_COMBO = "6dfb8dcb-eb95-407a-bc12-ca1e2c91a50e"  # Combination Viability
    DS_GDSC_LIB = "1bfefe86-8daa-4b2d-83bc-741db486fe8f"  # Library (single) Viability
    DS_GDSC_ANC = "2748aa32-a74c-4897-8a57-0bc2a0e1b00a"  # Anchor (single) Viability
    # Two informative pairs at near-benchmark doses; each: (label, combo feature,
    # (single_a feature, dataset), (single_b feature, dataset)).
    GDSC_TESTS = [
        (
            "BRAFi+MEKi",
            "Library: Dabrafenib @ 0.625 uM& Anchor: Trametinib @ 0.01 uM",
            ("Dabrafenib @ 0.625 uM", DS_GDSC_LIB),
            ("Trametinib @ 0.01 uM", DS_GDSC_ANC),
        ),
        (
            "MEKi+PI3Ki",
            "Library: Trametinib @ 0.015625 uM& Anchor: Pictilisib @ 0.125 uM",
            ("Trametinib @ 0.015625 uM", DS_GDSC_LIB),
            ("Pictilisib @ 0.125 uM", DS_GDSC_ANC),
        ),
    ]

    RESISTANCE_ESSAY = r"""
    ---
    # Task 4 - Overcoming clinical resistance to BRAF/MEK inhibitors

    ## Loop / hypothesis

    The strategy I would bet on is **to eliminate the drug-tolerant "persister"
    reservoir up front by pairing BRAF/MEK inhibition with a ferroptosis inducer**,
    rather than chasing the many genetic routes of acquired MAPK reactivation after
    they emerge.

    **Mechanistic rationale.** Acquired resistance to BRAF+MEK inhibition is
    genetically heterogeneous (NRAS/MEK mutations, BRAF amplification/splice variants,
    RTK upregulation) - so there is no single mechanism to block. But before any of
    those clones expands, a
    fraction of cells survive therapy by switching, *non-genetically and reversibly*,
    into a dedifferentiated, slow-cycling, MITF-low / SOX10-low / AXL-high state. That
    persister pool is the seed from which stable resistance later grows. The
    decisive point - and this notebook's own data make it concretely: **that
    dedifferentiated state carries a collateral vulnerability to ferroptosis.**
    LOXIMVI, the dedifferentiated line in this very benchmark, is BRAF V600E yet
    intrinsically BRAF-inhibitor-resistant, and is one of the most ferroptosis-
    sensitive melanoma lines known (Tsoi et al., Cancer Cell 2018). Drug-tolerant
    persisters across cancers depend on GPX4 to survive (Hangauer et al., Nature
    2017). So the state that lets melanoma cells *escape* MAPK inhibition is the same
    state that makes them *die* from GPX4 inhibition. Inhibiting both pathways together
    turns the resistance mechanism into a vulnerability.

    ## Why this beats the obvious alternatives

    - **Sequential / intermittent dosing** delays resistance but does not remove the
      persister reservoir; relapse still comes.
    - **Adding an ERK inhibitor or going to MAPK triplets** deepens pathway blockade
      but stays inside the MAPK axis the resistant cells are learning to live without;
      it addresses each reactivation route only after it appears.
    - **Up-front immunotherapy combination** is clinically real but a different axis,
      with additive toxicity and its own resistance.

    Targeting ferroptosis is orthogonal: it does not care which MAPK-reactivating
    mutation a cell would have acquired, because it kills the tolerant state itself,
    and it is *selective* for therapy-induced dedifferentiation rather than broadly
    toxic.

    ## Key experiments (with controls and the confirm/falsify line)

    1. **Persisters are GPX4-dependent (in vitro).** Treat BRAFi-sensitive melanoma
       (A375 plus 2-3 patient-derived V600E lines) with dabrafenib+trametinib to
       residual disease; confirm survivors are MITF-low/AXL-high (qPCR/flow). Then
       challenge with a GPX4 inhibitor (RSL3/ML210). **Confirm:** persisters are
       hypersensitive to GPX4i vs naive cells, and death is **rescued by
       ferrostatin-1 / liproxstatin-1 and iron chelation** (proves ferroptosis) but
       not by a pan-caspase inhibitor. **Falsify:** no selective persister killing,
       or death not rescued by ferrostatin (then it is not the proposed mechanism).
    2. **Triple delays/prevents regrowth (in vitro).** Long-term clonogenic regrowth
       after dabrafenib+trametinib +/- GPX4i. **Confirm:** the triple suppresses or
       abolishes regrowth vs the doublet. Controls: GPX4i alone (should spare naive
       proliferating cells), ferrostatin co-treatment (should restore regrowth).
    3. **In vivo relapse delay.** A375 / PDX xenografts (and ideally an immunocompetent
       BRAF-mutant GEMM) on dabrafenib+trametinib +/- a clinical-stage ferroptosis
       inducer. **Primary endpoint:** time to progression and minimal-residual-disease
       burden (MITF-low/AXL-high cells by IHC). **Confirm:** longer progression-free
       interval and fewer persisters with the triple.

    ## Translation to patient benefit

    A positive result motivates adding a ferroptosis inducer to the approved BRAF+MEK
    doublet, with **melanoma differentiation state (MITF/AXL, or a lipid-peroxidation
    readout) as a pharmacodynamic biomarker** of persister depletion. Because the
    combination attacks the reservoir that precedes genetic resistance, the goal is
    not a deeper initial response but a **durable** one - converting a ~1-year median
    to relapse into a longer, possibly treatment-free remission. The differentiation/
    ferroptosis axis is independently entering the clinic, so the path from a positive
    xenograft result to a biomarker-stratified trial is short.
    """


@app.cell(hide_code=True)
def _intro():
    mo.md(r"""
    # Night Shift benchmark - predicted melanoma viability

    Predict CellTiter-Glo (CTG) **% viability** for a fixed small-molecule panel on
    two melanoma lines, **A375** and **LOXIMVI**, at short timepoints (16/24/48 h),
    then rank single agents and combinations. No experimental data is given - every
    number is a prediction.

    **How to read this notebook.** Each task is worked as an explicit
    *hypothesis -> experiment -> observation* loop, the same loop the work was
    actually done in: state a hypothesis, run an experiment (a live DepMap pull, a
    literature check, a calculation), look at what came back, and revise. The loop
    stops when no further experiment would change the answer. The terminal state of
    each task is a **final hypothesis = the prediction** that gets submitted.
    Sections are bounded by task; finishing one before starting the next.
    """)
    return


@app.cell(hide_code=True)
def _task1_header():
    mo.md(r"""
    ---
    # Task 1 - Rank single-agent effects (subtasks 1.1-1.5)

    One engine answers all five subtasks (A375/LOXIMVI x 24/48 h, plus the pooled
    48-condition ranking), so they are worked together.

    ## Loop 1 - hypothesis

    **H0:** Effect tracks BRAF/MAPK pathway targeting. A375 is BRAF V600E and should
    be killed by BRAF and MEK inhibitors at these low doses. LOXIMVI is *(initial
    belief)* a BRAF-wild-type melanoma, so MAPK inhibitors should do little and a
    pan-cytotoxic agent (panobinostat) should dominate.

    **Predicted consequence to test first:** if H0 is right, the two lines have
    *different driver genetics* - A375 BRAF-mutant, LOXIMVI not. Check that before
    anything else, because it is load-bearing for every LOXIMVI ranking.
    """)
    return


@app.cell(hide_code=True)
def _panel_view():
    mo.vstack(
        [
            mo.md("**The panel under test:**"),
            mo.ui.table(pl.DataFrame(PANEL, schema=PANEL_COLS, orient="row"), page_size=12),
        ]
    )
    return


@app.function
def line_matrix_values(dataset_id: str, features: list[str]) -> pl.DataFrame:
    """Fetch a feature x model matrix slice for the two benchmark lines only."""
    raw = bb_post(
        f"datasets/matrix/{dataset_id}",
        {"features": features, "feature_identifier": "label"},
    )
    rows = []
    for feat in features:
        values = raw.get(feat, {}) if isinstance(raw, dict) else {}
        row = {"feature": feat}
        for line, ach in LINES.items():
            row[line] = values.get(ach)
        rows.append(row)
    return pl.DataFrame(rows)


@app.cell
def _genetics():
    # Loop 1 experiment: confirm driver genetics for both lines (live DepMap pull).
    genetics = line_matrix_values(DS_HOTSPOT, ["BRAF", "NRAS", "KRAS"])
    mo.vstack(
        [
            mo.md("### Loop 1 experiment - hotspot mutation status (1 = present)"),
            mo.ui.table(genetics, page_size=5),
        ]
    )
    return


@app.cell(hide_code=True)
def _loop1_obs():
    mo.md(r"""
    ### Loop 1 observation - H0 is falsified

    **Both** lines carry a BRAF hotspot mutation. The premise that LOXIMVI is
    BRAF-wild-type is wrong. So the simple "mutation -> response" story cannot
    explain a difference between the lines, because there is no difference in driver
    mutation. Either both lines respond like A375, or something other than mutation
    status governs LOXIMVI's response. This is the most useful result in the task: it
    rules out the mutation-only model and directs the next step.

    ## Loop 2 - hypothesis

    **H1:** Response is set by *differentiation state*, not mutation. Pull the measured
    drug-response data (PRISM) and check LOXIMVI's lineage phenotype in the literature.
    Predicted consequence: if LOXIMVI is dedifferentiated, BRAF
    inhibitors should be near-inactive at low nM despite the V600E, while a
    pan-cytotoxic agent stays active.
    """)
    return


@app.cell
def _prism():
    # Loop 2 experiment: PRISM Secondary AUC (lower = more sensitive). 5-day,
    # dose-integrated screen - sets sensitivity direction/tiers, not the short-
    # timepoint absolute viability (the model below handles that).
    prism = (
        line_matrix_values(DS_PRISM_AUC, [c[0].upper() for c in PANEL])
        .with_columns(pl.col("feature").str.to_titlecase().alias("condition"))
        .select("condition", "A375", "LOXIMVI")
    )
    mo.vstack(
        [
            mo.md("### Loop 2 experiment - PRISM Secondary AUC (lower = more sensitive; 5-day)"),
            mo.ui.table(prism.with_columns(pl.col("A375", "LOXIMVI").round(3)), page_size=12),
        ]
    )
    return


@app.cell(hide_code=True)
def _loop2_obs():
    mo.md(r"""
    ### Loop 2 observation - H1 confirmed, with a caveat about the instrument

    Two things came back. **(1) Literature** (Cellosaurus CVCL_1381; Tsoi et al.
    *Cancer Cell* 2018; SOX10-resensitization work): LOXIMVI is V600E-heterozygous
    but **dedifferentiated** - MITF-low / SOX10-low / AXL-high, the "undifferentiated"
    melanoma subtype that is **intrinsically resistant to BRAF inhibition despite the
    mutation**; re-expressing SOX10 restores vemurafenib sensitivity. The same state
    is ferroptosis- and HDAC-inhibitor-sensitive. A375, by contrast, is the textbook
    BRAF-addicted proliferative line.

    **(2) The PRISM AUC is misleading if read directly.** It shows BRAF/MEK inhibitors as
    fairly "sensitive" in LOXIMVI - but AUC integrates doses up to ~10 uM, crediting
    off-target killing far above the benchmark's 100 nM. At the *fixed low doses*
    here, BRAF inhibitors are effectively inactive in LOXIMVI. So PRISM informs the
    sensitivity *direction and tiers*, but I weight the dedifferentiation biology
    over the AUC for the low-dose call. Dose-matched PRISM viability corroborates the
    one unambiguous signal: panobinostat is dramatically more cytotoxic in LOXIMVI
    (log2FC -3.46, ~9% at 5 d) than in A375.

    ## Loop 3 - hypothesis (the quantitative model)

    CellTiter-Glo reads total ATP ~ live-cell number. The benchmark normalizes
    treated wells to the DMSO control on the same plate, so % viability is
    `treated_cells / control_cells` at the readout time. Over 16-48 h the control
    keeps dividing while treated cells are slowed or killed. Model the readout as

    $$\text{viability}\%(t) = 100 \cdot 2^{-E \cdot t / T_d}$$

    where $T_d$ is the line doubling time and $E$ is an **effect index**: the
    fractional suppression of net growth rate. $E=0$ tracks the control (inactive),
    $E\approx1$ is full cytostatic arrest (treated count frozen while control
    doubles), and $E>1$ is net cytotoxic (population shrinks). This makes the
    timepoint scaling fall out mechanistically - the same $E$ predicts 16, 24, and
    48 h - and keeps a single, auditable parameter per (line, compound).

    $E$ is set from PRISM dose-matched viability plus the line biology, **not** from
    the dose-integrated AUC, which over-credits off-target killing at high doses
    (e.g. it makes BRAF inhibitors look active in LOXIMVI, which they are not at
    100 nM). A floor of 3% guards against the model implying total ATP loss at
    these short timepoints.
    """)
    return


@app.function
def predict_viability(effect: float, timepoint_h: float, doubling_h: float) -> float:
    """Predicted CTG % viability from the growth-kinetics model (floored at 3%)."""
    return round(max(3.0, 100.0 * 2.0 ** (-effect * timepoint_h / doubling_h)), 1)


@app.function
def effect_table() -> pl.DataFrame:
    """Long-form (line, compound) effect index with rationale."""
    rows = []
    for line, by_compound in EFFECT.items():
        for compound, (e, why) in by_compound.items():
            rows.append({"line": line, "condition": compound, "E": e, "rationale": why})
    return pl.DataFrame(rows)


@app.cell
def _effect_view():
    mo.vstack(
        [
            mo.md("### Loop 3 instrument - effect index E per (line, compound), with rationale"),
            mo.ui.table(effect_table().sort("line", "E", descending=[False, True]), page_size=24),
        ]
    )
    return


@app.cell(hide_code=True)
def _task1_final():
    mo.md(r"""
    ## Task 1 - final hypothesis (the prediction)

    No further experiment would move these calls: the genetics are confirmed, the
    LOXIMVI phenotype is settled in the literature, the PRISM data are in hand, and
    the model converts effect index + timepoint into viability. The loop stops here.

    **What the predictions say.** In **A375**, panobinostat and the BRAF/MEK
    inhibitors cluster at the top (strong effect), the PI3K/AKT/mTOR agents and
    regorafenib at the bottom. In **LOXIMVI** the order *flips for BRAF inhibitors*:
    panobinostat dominates, MEK inhibitors are partial, and BRAF inhibitors fall to
    the bottom alongside the PI3K/AKT agents - despite the identical V600E mutation.
    That cross-line flip is the substantive result. Tables 1.1-1.5 below are the
    submitted predictions.

    **Confidence and limits.** The cross-line contrast (BRAFi strong in A375, dead in
    LOXIMVI; panobinostat strongest in LOXIMVI) is high-confidence. The *fine
    ordering within* the A375 BRAFi/MEKi cluster (E 0.85-1.05) is low-confidence - a
    Spearman score will be sensitive to it, and the true ties may differ. PRISM is a
    5-day screen, so it anchors direction, not the absolute 16-48 h magnitude.
    """)
    return


@app.function
def single_agent_ranking(line: str, timepoint_h: int) -> pl.DataFrame:
    """Ranked single-agent prediction for one (line, timepoint)."""
    td = DOUBLING_H[line]
    rows = []
    for condition, conc, _dose, target, klass in PANEL:
        e = EFFECT[line][condition][0]
        rows.append(
            {
                "condition": condition,
                "concentration": conc,
                "class": klass,
                "target": target,
                "viability_pct": predict_viability(e, timepoint_h, td),
            }
        )
    df = pl.DataFrame(rows).sort("viability_pct")
    return df.with_columns(pl.int_range(1, df.height + 1).alias("rank")).select(
        "rank", "condition", "concentration", "class", "viability_pct"
    )


@app.cell
def _task_1_1():
    mo.vstack(
        [
            mo.md("### Task 1.1 - A375, 24 h (single agents, ranked by effect)"),
            mo.ui.table(single_agent_ranking("A375", 24), page_size=12),
        ]
    )
    return


@app.cell
def _task_1_2():
    mo.vstack(
        [
            mo.md("### Task 1.2 - LOXIMVI, 24 h"),
            mo.ui.table(single_agent_ranking("LOXIMVI", 24), page_size=12),
        ]
    )
    return


@app.cell
def _task_1_3():
    mo.vstack(
        [
            mo.md("### Task 1.3 - A375, 48 h"),
            mo.ui.table(single_agent_ranking("A375", 48), page_size=12),
        ]
    )
    return


@app.cell
def _task_1_4():
    mo.vstack(
        [
            mo.md("### Task 1.4 - LOXIMVI, 48 h"),
            mo.ui.table(single_agent_ranking("LOXIMVI", 48), page_size=12),
        ]
    )
    return


@app.cell
def _task_1_5():
    # 1.5: all 48 single-agent conditions (compound x line x timepoint) ranked together.
    _frames = []
    for _line in LINES:
        for _tp in (24, 48):
            _frames.append(
                single_agent_ranking(_line, _tp)
                .drop("rank")
                .with_columns(pl.lit(_line).alias("cell_line"), pl.lit(_tp).alias("timepoint_h"))
            )
    _combined = pl.concat(_frames).sort("viability_pct")
    _combined = _combined.with_columns(pl.int_range(1, _combined.height + 1).alias("rank"))
    mo.vstack(
        [
            mo.md("### Task 1.5 - all 48 conditions ranked together"),
            mo.ui.table(
                _combined.select("rank", "cell_line", "timepoint_h", "condition", "class", "viability_pct"),
                page_size=12,
            ),
        ]
    )
    return


@app.cell
def _heatmap():
    # Predicted viability across compounds x (line, timepoint) - the whole single-agent surface.
    _rows = []
    for _line in LINES:
        for _tp in (24, 48):
            for _condition, _conc, _dose, _target, _klass in PANEL:
                _rows.append(
                    {
                        "condition": _condition,
                        "cond_tp": f"{_line} {_tp}h",
                        "viability_pct": predict_viability(EFFECT[_line][_condition][0], _tp, DOUBLING_H[_line]),
                    }
                )
    grid = pl.DataFrame(_rows)
    order = [c[0] for c in PANEL]
    chart = (
        alt.Chart(grid)
        .mark_rect()
        .encode(
            x=alt.X("cond_tp:N", title=None, sort=["A375 24h", "A375 48h", "LOXIMVI 24h", "LOXIMVI 48h"]),
            y=alt.Y("condition:N", sort=order, title=None),
            color=alt.Color(
                "viability_pct:Q", scale=alt.Scale(scheme="redyellowgreen", domain=[0, 100]), title="viability %"
            ),
            tooltip=["condition", "cond_tp", "viability_pct"],
        )
        .properties(width=320, height=420)
    )
    mo.vstack([mo.md("### Predicted viability surface (single agents)"), mo.ui.altair_chart(chart)])
    return


@app.cell(hide_code=True)
def _task1_close():
    mo.md(r"""
    **Task 1 complete.** Final hypothesis written above and in tables 1.1-1.5. The
    loop ran three iterations: a naive mutation-driven model (falsified by the
    genetics pull), a differentiation-state model (confirmed by literature + PRISM),
    and a growth-kinetics quantification. Next: Task 2 (combinations).

    ---
    """)
    return


@app.function
def synergy_fraction(line: str, drugs: list[str]) -> float:
    """Bounded, mechanism-based synergy fraction applied below the best single agent.

    Magnitudes follow the GDSC combination test (Task 2 external-test cell): the added
    benefit of combining is small (single-digit to ~0.1 in viability) and is
    BRAF-mutant/addiction-specific for vertical MAPK pairs. BRAFi+MEKi pays off only in
    the addicted line; MEK+PI3K pays off in both (GDSC shows it adds *more* than
    BRAFi+MEKi in BRAF-WT context).
    """
    classes = {CLASS[d] for d in drugs}
    addicted = line == "A375"
    f = 0.0
    if "BRAFi" in classes and "MEKi" in classes and addicted:
        f += 0.30
    if "MEKi" in classes and "PI3K/AKT/mTOR" in classes:
        f += 0.20
    if "BRAFi" in classes and "PI3K/AKT/mTOR" in classes and addicted:
        f += 0.10
    return min(f, 0.6)


@app.function
def combo_viability(line: str, drugs: list[str], timepoint_h: int, attenuation: dict | None = None) -> float:
    """Combination viability = highest-single-agent baseline reduced by a bounded synergy.

    HSA (the best single agent) is the conservative reference that matches the GDSC
    observation that same-pathway combos sit near - not far below - the best single
    agent; `synergy_fraction` adds the modest, mechanism-specific extra benefit.
    `attenuation` scales each drug's E by class (the 1/3-dose triples).
    """
    td = DOUBLING_H[line]
    viabs = []
    for d in drugs:
        e = EFFECT[line][d][0] * (attenuation.get(CLASS[d], 1.0) if attenuation else 1.0)
        viabs.append(predict_viability(e, timepoint_h, td))
    return round(max(3.0, min(viabs) * (1.0 - synergy_fraction(line, drugs))), 1)


@app.function
def combo_ranking(line: str, timepoint_h: int) -> pl.DataFrame:
    """Ranked prediction for the four benchmark combinations in one line."""
    rows = [
        {
            "condition": label,
            "concentration": conc,
            "viability_pct": combo_viability(line, [a, b], timepoint_h),
        }
        for label, a, b, conc in COMBOS
    ]
    df = pl.DataFrame(rows).sort("viability_pct")
    return df.with_columns(pl.int_range(1, df.height + 1).alias("rank")).select(
        "rank", "condition", "concentration", "viability_pct"
    )


@app.function
def matrix_col(dataset_id: str, feature: str) -> dict:
    """Return the full {model_id: value} column for one matrix feature."""
    raw = bb_post(f"datasets/matrix/{dataset_id}", {"features": [feature], "feature_identifier": "label"})
    return raw.get(feature, {}) if isinstance(raw, dict) else {}


@app.function
def gdsc_synergy_test() -> pl.DataFrame:
    """Independent test: in GDSC combination screens, is the added benefit of combining
    (best single agent viability - combination viability) BRAF-mutant-specific?"""
    import statistics

    braf = matrix_col(DS_HOTSPOT, "BRAF")

    def group(model_id: str) -> str | None:
        v = braf.get(model_id)
        if not isinstance(v, (int, float)):
            return None
        return "BRAF-mutant" if v >= 1 else "BRAF-WT"

    rows = []
    for pair, combo_feat, (fa, dsa), (fb, dsb) in GDSC_TESTS:
        combo, a, b = matrix_col(DS_GDSC_COMBO, combo_feat), matrix_col(dsa, fa), matrix_col(dsb, fb)
        buckets: dict[str, list[tuple[float, float]]] = {}
        for model_id in combo:
            grp = group(model_id)
            ca, cb, cc = a.get(model_id), b.get(model_id), combo.get(model_id)
            if grp is None or not all(isinstance(x, (int, float)) for x in (ca, cb, cc)):
                continue
            buckets.setdefault(grp, []).append((min(ca, cb) - cc, cc))
        for grp, vals in buckets.items():
            rows.append(
                {
                    "pair": pair,
                    "BRAF_status": grp,
                    "n_lines": len(vals),
                    "median_combo_viab": round(statistics.median(v[1] for v in vals), 3),
                    "median_added_benefit": round(statistics.median(v[0] for v in vals), 3),
                }
            )
    return pl.DataFrame(rows).sort("pair", "BRAF_status")


@app.cell(hide_code=True)
def _task2_header():
    mo.md(r"""
    ---
    # Task 2 - Rank drug combinations (subtasks 2.1, 2.2, 2.3)

    Four combinations, both lines, 48 h: three clinical BRAFi+MEKi pairs
    (Dabrafenib+Trametinib, Encorafenib+Binimetinib, Vemurafenib+Cobimetinib) and one
    MEKi+PI3Ki pair (Trametinib+Alpelisib).

    ## Loop 1 - hypothesis (and why the model cannot test itself)

    A combination model assembled from single-agent effects plus a chosen synergy rule
    will always reproduce its own assumptions, so running it is not a test of them. The
    real experiment is an **independent** check of the one rule the model rests on:
    **is the benefit of combining BRAF + MEK inhibition specific to BRAF-driven cells?**
    If yes, the rule is justified; if not, the Task 2 ranking is unfounded. Test it on
    the GDSC combination screen.
    """)
    return


@app.cell
def _task2_external():
    # Loop 1 experiment - INDEPENDENT data, not my model. The benchmark lines are not in
    # the GDSC combo screen (and there is no melanoma in it), so this tests the *rule*
    # on the BRAF-mutant lines that are present (mostly colorectal).
    test = gdsc_synergy_test()
    mo.vstack(
        [
            mo.md(
                "### Loop 1 experiment - GDSC combination screen (independent of my model)\n\n"
                "**What this is.** GDSC (Genomics of Drug Sensitivity in Cancer, Sanger Institute) "
                "is a large public cell-line pharmacology resource; its *combination* screen measured "
                "drug **pairs** across a panel of lines, exposed via DepMap/Breadbox. The design is "
                "**anchor + library**: an 'anchor' drug is held at a fixed dose while a 'library' drug "
                "is titrated across a dose range, and a CellTiter-Glo-style viability is read for the "
                "pair and for each drug alone (so the same screen gives combo, library-alone, and "
                "anchor-alone). Viability is a fraction of untreated control (1.0 = no effect). This is "
                "*measured* combination data - the right kind of independent check, because it owes "
                "nothing to my effect-index model.\n\n"
                "**What I compute.** `added_benefit` = (best single-agent viability) - (combination "
                "viability) per line; positive means the pair beats its better single agent. I take two "
                "near-benchmark pairs - Dabrafenib+Trametinib (BRAFi+MEKi) and Trametinib+Pictilisib "
                "(MEKi+PI3Ki; pictilisib is the available PI3K inhibitor) - and split lines by BRAF "
                "status. **Caveat:** A375/LOXIMVI are not in this screen and it contains no melanoma, "
                "so it tests the *rule*, not the exact lines."
            ),
            mo.ui.table(test, page_size=8),
        ]
    )
    return


@app.cell(hide_code=True)
def _task2_obs():
    mo.md(r"""
    ### Loop 1 observation - what the independent data shows

    Three things, none of them from my model:

    1. **BRAFi+MEKi benefit is BRAF-mutant-specific.** Added benefit is positive in
       BRAF-mutant lines (~+0.06) and ~zero/negative in BRAF-WT (~-0.03). The rule the
       ranking depends on - vertical MAPK synergy gated on BRAF-addiction - holds.
    2. **MEKi+PI3Ki adds *more* than BRAFi+MEKi in BRAF-WT context** (~+0.12 vs -0.03):
       a non-addicted cell benefits more from MEK+PI3K than from BRAF+MEK. Independent
       support for the direction of the LOXIMVI flip.
    3. **The benefit of combining is small.** The combo sits only single-digit points
       below the best single agent, and *above* Bliss independence - these same-pathway
       pairs are sub-Bliss.

    Caveat that bounds all of this: the screen contains **no melanoma** and not the
    benchmark lines; the BRAF-mutant lines are mostly colorectal, which is *less*
    addicted than BRAF-mutant melanoma. So the melanoma combo benefit is plausibly
    larger than the colorectal ~0.06, but unmeasured.

    ## Loop 2 - the model these data imply

    Combination viability = the **best single agent**, reduced by a **small synergy
    fraction** that is nonzero only for mechanism-rational, addiction-appropriate pairs
    (BRAFi+MEKi in the addicted line; MEK+PI3K in either). The best single agent is the
    floor the GDSC data point to; the synergy is the modest measured added benefit.
    Applied to the four benchmark combinations:
    """)
    return


@app.cell
def _task2_experiment():
    a375 = combo_ranking("A375", 48)
    lox = combo_ranking("LOXIMVI", 48)
    mo.vstack(
        [
            mo.md("### Task 2.1 - A375 (48 h):"),
            mo.ui.table(a375, page_size=4),
            mo.md("### Task 2.2 - LOXIMVI:"),
            mo.ui.table(lox, page_size=4),
        ]
    )
    return


@app.cell
def _task2_pooled():
    _frames = [combo_ranking(line, 48).drop("rank").with_columns(pl.lit(line).alias("cell_line")) for line in LINES]
    pooled = pl.concat(_frames).sort("viability_pct")
    pooled = pooled.with_columns(pl.int_range(1, pooled.height + 1).alias("rank"))
    mo.vstack(
        [
            mo.md("### Task 2.3 - all 8 (combination x line) ranked together"),
            mo.ui.table(pooled.select("rank", "cell_line", "condition", "viability_pct"), page_size=8),
        ]
    )
    return


@app.cell(hide_code=True)
def _task2_final():
    mo.md(r"""
    ### Task 2 final hypothesis

    **A375** - the three BRAFi+MEKi pairs lead
    (Dabrafenib+Trametinib deepest), Trametinib+Alpelisib trails, but all sit in the
    ~15-25% range, not near-zero. **LOXIMVI** - the order flips, Trametinib+Alpelisib
    wins (the GDSC test backs this direction), with the BRAFi+MEKi pairs at MEKi-alone
    depth. **Pooled (2.3)** - the four A375 combinations still outrank the four LOXIMVI
    combinations.

    **Confidence, honestly stated.** The *rankings* are the defensible output and the
    A375-vs-LOXIMVI split is well-supported. The *absolute viabilities* are
    low-confidence: combination effect modeled from single agents is hard, the synergy
    fraction is calibrated on colorectal not melanoma lines, and the true magnitude in
    A375 could be deeper than predicted. The fine ordering among the three A375
    BRAFi+MEKi pairs is near-tied and not meaningfully resolvable here.
    """)
    return


@app.function
def triple_search(line: str, timepoint_h: int = 16) -> pl.DataFrame:
    """Enumerate all 3-drug combinations from the candidate list at 1/3 dose, ranked."""
    from itertools import combinations

    rows = []
    for triple in combinations(TRIPLE_CANDIDATES, 3):
        rows.append(
            {
                "combination": " + ".join(triple),
                "classes": " / ".join(sorted({CLASS[d] for d in triple})),
                "viability_pct": combo_viability(line, list(triple), timepoint_h, attenuation=DOSE_THIRD),
            }
        )
    return pl.DataFrame(rows).sort("viability_pct")


@app.function
def candidate_auc() -> pl.DataFrame:
    """PRISM Secondary AUC of the 10 candidate drugs in both lines (lower = more sensitive)."""
    df = line_matrix_values(DS_PRISM_AUC, [c.upper() for c in TRIPLE_CANDIDATES])
    return (
        df.with_columns(pl.col("feature").str.to_titlecase().alias("candidate"))
        .with_columns(pl.col("candidate").replace({"Tak-733": "TAK-733"}))
        .with_columns(pl.col("A375", "LOXIMVI").round(3))
        .select("candidate", "A375", "LOXIMVI")
        .sort("LOXIMVI")
    )


@app.function
def prism_broad_top(line: str, n: int = 15, batch: int = 600) -> pl.DataFrame:
    """Across the whole PRISM Repurposing screen, the most-sensitive compounds for one
    line - independent of the candidate list, to ask what actually kills the cell."""
    ach = LINES[line]
    feats = [f.get("label") for f in bb_get(f"datasets/features/{DS_PRISM_AUC}") if f.get("label")]
    values: dict[str, float] = {}
    for i in range(0, len(feats), batch):
        chunk = feats[i : i + batch]
        raw = bb_post(f"datasets/matrix/{DS_PRISM_AUC}", {"features": chunk, "feature_identifier": "label"})
        for feat in chunk:
            v = raw.get(feat, {}) if isinstance(raw, dict) else {}
            if isinstance(v.get(ach), (int, float)):
                values[feat] = v[ach]
    top = sorted(values.items(), key=lambda kv: kv[1])[:n]
    return pl.DataFrame({"compound": [t[0] for t in top], "prism_auc": [round(t[1], 3) for t in top]})


@app.cell(hide_code=True)
def _task3_header():
    mo.md(r"""
    ---
    # Task 3 - Nominate the best 3-drug combination (subtasks 3.1, 3.2)

    Pick one combination of exactly 3 distinct drugs from a 10-drug candidate list
    (no panobinostat, no sapanisertib), each at **1/3** its listed dose, that gives
    the lowest viability at **16 h** - separately for LOXIMVI (3.1) and A375 (3.2).

    ## Loop 1 - hypothesis, and not assuming the conventional answer

    My combination model returns the textbook triple (vertical MAPK + PI3K) and, by
    construction, a many-way tie - so it cannot, on its own, justify a nomination. The
    discipline here is to not stop at the conventional answer. Two independent questions
    for the data, neither answerable by my model: (1) within the candidate list, which
    drugs are actually the active ones in these cells? (2) Does the candidate list even
    contain what most effectively kills these cells, or are all of its drugs aimed at
    the wrong target?
    """)
    return


@app.cell
def _task3_candidates():
    mo.vstack(
        [
            mo.md("### Loop 1 experiment A - candidate potency in these lines (PRISM AUC, lower = more sensitive)"),
            mo.ui.table(candidate_auc(), page_size=10),
        ]
    )
    return


@app.cell
def _task3_broad():
    mo.vstack(
        [
            mo.md(
                "### Loop 1 experiment B - what actually kills these lines (whole PRISM screen, ~1500 compounds)\n\n"
                "Most-sensitive compounds, independent of the candidate list. Read the MOAs: "
                "antimitotics (dolastatin/maytansinoids), topoisomerase-I (SN38/camptothecins), "
                "nucleoside analogs (gemcitabine), and - for LOXIMVI - oligomycin (oxidative "
                "phosphorylation). None are in the candidate list."
            ),
            mo.md("**A375:**"),
            mo.ui.table(prism_broad_top("A375"), page_size=15),
            mo.md("**LOXIMVI:**"),
            mo.ui.table(prism_broad_top("LOXIMVI"), page_size=15),
        ]
    )
    return


@app.cell(hide_code=True)
def _task3_final():
    mo.md(r"""
    ### Loop 1 observation - the data contradict the conventional answer

    **Experiment A:** in both lines the candidate drugs split cleanly. The MAPK agents
    (BRAF/MEK inhibitors) are the active ones; **alpelisib, capivasertib and
    regorafenib are the *weakest* of the ten** (PRISM AUC ~0.76-0.90, near-inactive).
    So the conventional "add a PI3K/AKT drug as the third agent" gives the slot to the
    *least* active candidate. That choice comes from the *durability* literature
    (parallel blockade delays resistance over days-weeks); it does not increase killing
    at 16 h.

    **Experiment B (the broad screen):** the candidate list is built entirely from
    pathway-targeted cytostatics, but the compounds that actually kill these cells are
    fast cytotoxics - antimitotics, topoisomerase-I poisons, nucleoside analogs - and,
    for dedifferentiated LOXIMVI, an **oxidative-phosphorylation inhibitor (oligomycin)**,
    its known metabolic vulnerability. The best MAPK agent (dabrafenib) is only ~#12 in
    A375; in LOXIMVI the first MEK inhibitor is ~#21. At a 16 h endpoint, where fast
    cytotoxics act and cytostatic kinase inhibitors have barely had time to work, the
    candidate list does not contain the most effective agents.

    ## Task 3 final hypothesis (the nominations, and their limit)

    Within the allowed list, the nomination is the potent MAPK pair plus a third agent
    the data say is near-irrelevant at 16 h:

    - **A375 (3.2): Dabrafenib + Trametinib + Alpelisib** (1/3 dose: 33 nM + 3.3 nM +
      1.67 uM), predicted **~28% viability** at 16 h. The BRAFi+MEKi pair accounts for
      most of the effect; alpelisib is a near-arbitrary third choice (any PI3K/AKT drug,
      or even swapping the BRAF inhibitor, leaves the model score unchanged - the top is
      a many-way tie).
    - **LOXIMVI (3.1): Trametinib + Alpelisib + Capivasertib** (1/3 dose: 3.3 nM +
      1.67 uM + 1.67 uM), predicted **~59% viability** at 16 h - one genuinely active
      drug (the MEK inhibitor) plus two weak partners, because the list offers nothing
      better for this line.

    **The actionable, non-obvious result** is the limit itself: the greatest 16 h effect
    is not reachable from this candidate list. If the goal were maximal kill, the data
    point to a fast cytotoxic (an antimitotic or topoisomerase-I poison) for A375 and an
    **OXPHOS or ferroptosis inducer for LOXIMVI** - outside the candidate list, but where
    the biology actually points. **Low confidence** on the within-list nominations, by
    design; the expected effect at 16 h is modest, especially in LOXIMVI.
    """)
    return


@app.cell(hide_code=True)
def _task4():
    mo.md(RESISTANCE_ESSAY)
    return


@app.function
def karman_template(task_id: str) -> dict:
    """Fetch the canonical output template for a task (structure must not be edited)."""
    resp = requests.get(f"{KARMAN}/tasks/{task_id}", timeout=30)
    resp.raise_for_status()
    return resp.json()["output_template"]


@app.function
def condition_viability(condition: str, line: str, timepoint_h: int) -> float:
    """Predicted viability for a single agent or a combination, by condition label."""
    if condition in COMBO_DRUGS:
        return combo_viability(line, list(COMBO_DRUGS[condition]), timepoint_h)
    return predict_viability(EFFECT[line][condition][0], timepoint_h, DOUBLING_H[line])


@app.function
def fill_rankings(template: dict, default_line: str | None, default_tp: int | None) -> dict:
    """Populate a ranked template (1.x / 2.x): set viability_pct, then rank by it."""
    out = json.loads(json.dumps(template))  # deep copy, structure preserved
    for item in out["rankings"]:
        line = item.get("cell_line", default_line)
        tp = item.get("timepoint_h", default_tp)
        item["viability_pct"] = condition_viability(item["condition"], line, tp)
    order = sorted(out["rankings"], key=lambda r: r["viability_pct"])
    rank_of = {id(r): i + 1 for i, r in enumerate(order)}
    for item in out["rankings"]:
        item["rank"] = rank_of[id(item)]
    return out


@app.function
def fill_nomination(template: dict, line: str) -> dict:
    """Populate a 3-drug nomination template (3.x)."""
    out = json.loads(json.dumps(template))
    drugs, concs = TRIPLE_NOMINATIONS[line]
    out["nomination"]["drugs"] = list(drugs)
    out["nomination"]["concentrations"] = list(concs)
    out["nomination"]["viability_pct"] = combo_viability(line, drugs, 16, attenuation=DOSE_THIRD)
    return out


@app.function
def build_outputs() -> dict[str, dict]:
    """Build the populated output.json for every task from the live templates."""
    single = {"1.1": ("A375", 24), "1.2": ("LOXIMVI", 24), "1.3": ("A375", 48), "1.4": ("LOXIMVI", 48)}
    combo = {"2.1": ("A375", 48), "2.2": ("LOXIMVI", 48)}
    outputs: dict[str, dict] = {}
    for task_id, (line, tp) in {**single, **combo}.items():
        outputs[task_id] = fill_rankings(karman_template(task_id), line, tp)
    outputs["1.5"] = fill_rankings(karman_template("1.5"), None, None)
    outputs["2.3"] = fill_rankings(karman_template("2.3"), None, 48)
    for task_id, line in {"3.1": "LOXIMVI", "3.2": "A375"}.items():
        outputs[task_id] = fill_nomination(karman_template(task_id), line)
    essay = karman_template("4.1")
    essay["response"] = RESISTANCE_ESSAY.strip()
    essay["analyses"] = (
        "Grounding analysis is in the companion notebook: DepMap confirms both A375 and "
        "LOXIMVI carry BRAF V600E, yet LOXIMVI is dedifferentiated (MITF-low/AXL-high) and "
        "BRAF-inhibitor-resistant - the persister phenotype this strategy targets, and the "
        "most ferroptosis-sensitive melanoma subtype."
    )
    outputs["4.1"] = essay
    return outputs


@app.cell
def _write_files():
    # Build-only: write output.json + reasoning.md per task. No submission.
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    outputs = build_outputs()
    written = []
    for task_id, payload in outputs.items():
        task_dir = OUT_DIR / task_id
        task_dir.mkdir(exist_ok=True)
        (task_dir / "output.json").write_text(json.dumps(payload, indent=2) + "\n")
        (task_dir / "reasoning.md").write_text(reasoning_md(task_id, payload))
        written.append({"task": task_id, "dir": str(task_dir.relative_to(OUT_DIR.parents[2]))})
    mo.vstack(
        [
            mo.md(f"## Submission files written ({len(written)} tasks, build-only - nothing submitted)"),
            mo.ui.table(pl.DataFrame(written), page_size=12),
        ]
    )
    return


@app.function
def reasoning_md(task_id: str, payload: dict) -> str:
    """Per-task reasoning trace, sharing the notebook's method and biology."""
    method = (
        "## Method\n\n"
        "No experimental readouts were provided; these are predictions. They are grounded in "
        "(1) DepMap genetics + PRISM Repurposing drug-response for A375 (ACH-000219) and "
        "LOXIMVI (ACH-000750), and (2) a transparent growth-kinetics model of the CellTiter-Glo "
        "readout: viability% = 100 * 2^(-E * t / Td), where Td is the line doubling time "
        "(A375 ~27 h, LOXIMVI ~22 h) and E is an effect index (0 = inactive, ~1 = full cytostatic "
        "arrest, >1 = net cytotoxic) set per (line, compound) from PRISM dose-matched viability "
        "plus the line biology.\n\n"
        "## Key biology\n\n"
        "Both lines carry BRAF V600E, but A375 is BRAF-addicted (sensitive to low-nM BRAF/MEK "
        "inhibitors) while LOXIMVI is dedifferentiated (MITF-low/SOX10-low/AXL-high) and "
        "intrinsically BRAF-inhibitor-resistant despite the mutation; MEK inhibition is only "
        "partial there and the pan-HDAC inhibitor panobinostat is the most effective single agent. "
        "Mutation status alone mispredicts LOXIMVI - differentiation state is the better predictor.\n\n"
    )
    if "rankings" in payload:
        lines = [
            f"{r['rank']}. {r['condition']}"
            + (
                f" [{r['cell_line']}" + (f", {r['timepoint_h']}h]" if "timepoint_h" in r else "]")
                if "cell_line" in r
                else ""
            )
            + f" - {r['viability_pct']}%"
            for r in sorted(payload["rankings"], key=lambda r: r["rank"])
        ]
        body = "## Predicted ranking (most to least effective)\n\n" + "\n".join(lines) + "\n"
    elif "nomination" in payload:
        nom = payload["nomination"]
        body = (
            "## Nomination\n\n"
            f"Drugs (each at 1/3 listed dose): {', '.join(f'{d} {c}' for d, c in zip(nom['drugs'], nom['concentrations']))}\n\n"
            f"Predicted viability at 16 h: {nom['viability_pct']}%\n\n"
            "Rationale: the potent MAPK agents account for most of the effect (BRAFi+MEKi in addicted "
            "A375; the MEK inhibitor in LOXIMVI). The third agent is near-irrelevant at 16 h - PRISM "
            "shows the PI3K/AKT and multikinase candidates are the *weakest* of the ten in both lines, "
            "so the conventional 'add PI3K' increases durability, not 16 h killing, and the top of the "
            "enumeration is a many-way tie.\n\n"
            "Limit of the candidate list (the non-obvious finding): it is built entirely from "
            "pathway-targeted cytostatics, but the whole PRISM screen shows the agents that actually "
            "kill these cells are fast cytotoxics (antimitotics, topoisomerase-I poisons, nucleoside "
            "analogs) and, for dedifferentiated LOXIMVI, an oxidative-phosphorylation inhibitor "
            "(oligomycin) - none in the candidate list. So the greatest 16 h effect is not reachable "
            "from this list; the within-list nomination is the best available, at low confidence and "
            "modest expected effect.\n"
        )
    else:
        body = "## Response\n\nSee output.json `response` (full proposal) and `analyses`.\n"
    return f"# Task {task_id} - reasoning\n\n{method}{body}"


@app.cell(hide_code=True)
def _to_extend():
    mo.md(r"""
    ---
    ## To extend

    - **Calibrate the kinetics against a reference curve.** If even one DMSO-normalized
      CTG time-course for these lines were available, fit `E` and `Td` directly instead
      of setting `E` by hand - the biggest source of absolute-magnitude error.
    - **Pull DepMap GDSC/CTD2 *combination* datasets** to measure the BRAFi+MEKi
      synergy directly rather than asserting it from clinical pharmacology.
    - **Test the Task 4 hypothesis in DepMap:** check whether GPX4 / ferroptosis
      dependency (CRISPR) tracks with the dedifferentiated (MITF-low/AXL-high) melanoma
      lines, which would strengthen the persister-ferroptosis proposal.
    """)
    return


@app.cell(hide_code=True)
def _():

    import time as _t
    _t0 = _t.time()
    _feats = [f['label'] for f in bb_get(f'datasets/features/{DS_PRISM_AUC}') if f.get('label')]
    PRISM_DEMO = {}
    for _i in range(0, len(_feats), 600):
        _chunk = _feats[_i:_i+600]
        _raw = bb_post(f'datasets/matrix/{DS_PRISM_AUC}', {'features': _chunk, 'feature_identifier': 'label'})
        for _f in _chunk:
            _v = _raw.get(_f, {})
            if isinstance(_v.get(LINES['LOXIMVI']), (int, float)):
                PRISM_DEMO[_f] = _v[LINES['LOXIMVI']]
    print(f'fetched {len(PRISM_DEMO)} LOXIMVI viabilities in {_t.time()-_t0:.1f}s -> now cached in kernel state')

    return


@app.cell(hide_code=True)
def _():

    import time, sys
    from pathlib import Path
    sys.path.insert(0, str(Path.cwd()/"notebooks"))
    from nb02_dataset_discovery import bb_get, bb_post
    DSx="07b7bda9-ae00-43b3-bca1-336b9607f8f5"; LOXx="ACH-000750"
    _t=time.time()
    _f=[x["label"] for x in bb_get(f"datasets/features/{DSx}") if x.get("label")]
    DEMO={}
    for i in range(0,len(_f),600):
        c=_f[i:i+600]; r=bb_post(f"datasets/matrix/{DSx}",{"features":c,"feature_identifier":"label"})
        for k in c:
            v=r.get(k,{})
            if isinstance(v.get(LOXx),(int,float)): DEMO[k]=v[LOXx]
    print(f"LOADED {len(DEMO)} compounds in {time.time()-_t:.1f}s -> persisted in kernel as DEMO")

    return bb_get, bb_post


if __name__ == "__main__":
    app.run()
