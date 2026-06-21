# dmx - DepMap eXplore

An experiment in agent-driven scientific data exploration, built around the public [DepMap](https://depmap.org/portal/) Breadbox API - genome-scale perturbation screens, omics, drug sensitivity, and cell-line metadata for cancer-dependency analysis.

dmx is a curated catalog of [marimo](https://marimo.io) notebooks for cancer-dependency analysis, plus a thin skill that lets an agent compose new analyses from them.
Each notebook is both a runnable demonstration and a source of pure functions other notebooks can [import and reuse](https://docs.marimo.io/guides/reusing_functions/) directly.
Given a new cancer-dependency question, the agent picks relevant notebooks, composes their functions into a new notebook, executes it in a live kernel, and hands back a self-contained, re-runnable result.

Breadbox is the REST API underlying the DepMap Portal: find datasets, search genes and compounds, retrieve matrix values, define model contexts, run two-class comparisons, and inspect precomputed associations.
By default dmx targets `https://depmap.org/portal/breadbox`, which does not require an API key for these read-only examples.
A non-default `User-Agent` header is required - the portal nginx rejects the bare `python-requests` agent with a 403.

## The catalog

Each notebook ships a committed session snapshot under [`notebooks/__marimo__/session/`](notebooks/__marimo__/session/) so the molab preview renders cell outputs without re-executing.

| Notebook | Role | Preview |
|---|---|---|
| `nb01_orientation.py` | Landing page: what Breadbox exposes, a live dimension-types check, and what questions this catalog should answer first | [![Open in molab](https://marimo.io/molab-shield.svg)](https://molab.marimo.io/github/broadinstitute/dmx/blob/main/notebooks/nb01_orientation.py) |
| `nb02_dataset_discovery.py` | Breadbox primitives: `bb_get`, `bb_post`, dataset discovery, dimension search, feature listing | [![Open in molab](https://marimo.io/molab-shield.svg)](https://molab.marimo.io/github/broadinstitute/dmx/blob/main/notebooks/nb02_dataset_discovery.py) |
| `nb03_gene_dependency_profile.py` | Gene -> Chronos dependency profile across cell lines, joined to model metadata | [![Open in molab](https://marimo.io/molab-shield.svg)](https://molab.marimo.io/github/broadinstitute/dmx/blob/main/notebooks/nb03_gene_dependency_profile.py) |
| `nb04_context_comparison.py` | Lineage / mutation-defined context -> compare dependency distributions inside vs outside | [![Open in molab](https://marimo.io/molab-shield.svg)](https://molab.marimo.io/github/broadinstitute/dmx/blob/main/notebooks/nb04_context_comparison.py) |
| `nb05_association_query.py` | Gene dependency or drug sensitivity slice -> precomputed Breadbox association hits across modalities | [![Open in molab](https://marimo.io/molab-shield.svg)](https://molab.marimo.io/github/broadinstitute/dmx/blob/main/notebooks/nb05_association_query.py) |

`catalog.toml`'s `[[vignette]]` table is the machine-readable version the compose skill reads: each notebook, its reusable helpers, and what it does.

## Submissions

Composed notebooks that answer the [Karman](https://mcp.karmanai.org/tasks) melanoma drug-response benchmark - predict CellTiter-Glo % viability for A375 (BRAF-V600E) and LOXIMVI under a MAPK-inhibitor panel.
Three submissions, each deleting one layer of hand-tuning and replacing it with measured data; v2 and v3 reuse the v1 engine unchanged and only change how the effect floors are produced.
Each notebook is self-contained (it vendors the shared engine rather than importing siblings), so any row opens and re-runs as a standalone molab gist; v3 reads ChEMBL from a local copy when present and falls back to the ChEMBL REST API otherwise.

| Submission | Approach | Preview |
|---|---|---|
| `nightshift_v1_reasoning.py` | Reasoning baseline: pharmacology + DepMap-anchored *judgment* floors + a CTG-kinetics model | [![Open in molab](https://marimo.io/molab-shield.svg)](https://molab.marimo.io/github/broadinstitute/dmx/blob/main/notebooks/nightshift_v1_reasoning.py) |
| `nightshift_v2_data_floors.py` | Floors read straight off PRISM Secondary dose-level viability at the task dose, not judgment | [![Open in molab](https://marimo.io/molab-shield.svg)](https://molab.marimo.io/github/broadinstitute/dmx/blob/main/notebooks/nightshift_v2_data_floors.py) |
| `nightshift_v3_chembl_calibrated.py` | v2 floors slid to the exact task dose along a ChEMBL-IC50 Hill curve, with the 1/3-dose retention from the same curve | [![Open in molab](https://marimo.io/molab-shield.svg)](https://molab.marimo.io/github/broadinstitute/dmx/blob/main/notebooks/nightshift_v3_chembl_calibrated.py) |

Each writes its populated `output.json` + `reasoning.md` per task under `data/processed/nightshift_v*/`, with a `summary.json` envelope the index notebook can pick up.

## Getting started

This catalog follows the [vignette-catalog-skills](https://github.com/carpenter-singh-lab/vignette-catalog-skills) pattern.
The skill stores are gitignored, so a fresh clone has only `skills-lock.json`; restore the on-disk skill content first.
`npx skills update` won't create Claude Code's `.claude/skills/` on a fresh clone, so add the skills explicitly:

```bash
npx skills add carpenter-singh-lab/vignette-catalog-skills --agent claude-code -y
npx skills add marimo-team/marimo-pair --agent claude-code -y
```

This writes `.claude/skills/` as well as `.agents/skills/` (both gitignored); if a session is already open, run `/reload-skills` or restart it so the new skills register.

Then open Claude Code in this repo and ask to *get started* - the `vignette-catalog-setup` skill installs prereqs ([uv](https://docs.astral.sh/uv/) and the [marimo-pair](https://github.com/marimo-team/marimo-pair) skill), launches `nb01_orientation` in a live marimo kernel, and hands off to `vignette-catalog-compose-notebook` for the actual analysis.

To run setup by hand:

```bash
uv --version  # or: curl -LsSf https://astral.sh/uv/install.sh | sh
npx skills add marimo-team/marimo-pair --agent claude-code -y
uvx marimo edit --sandbox notebooks/nb01_orientation.py
```

Always launch notebooks with `--sandbox` so the PEP 723 inline dependencies are provisioned.

## Related

Built on the [vignette-catalog-skills](https://github.com/carpenter-singh-lab/vignette-catalog-skills) pattern: a thin skill plus a catalog of parameterized marimo notebooks an agent composes from in a live kernel.
Sibling catalogs of the same pattern: [jx](https://github.com/broadinstitute/jx) for JUMP Cell Painting, [fgx](https://github.com/broadinstitute/fgx) for FinnGenie human genetics, and [prx](https://github.com/broadinstitute/prx) for PROSPECT chemical genetics.
