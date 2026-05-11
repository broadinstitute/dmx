# dmx — DepMap eXplore

An experiment in agent-driven scientific data exploration, built around the public [DepMap](https://depmap.org/portal/) Breadbox API — genome-scale perturbation screens, omics, drug sensitivity, and cell-line metadata for cancer-dependency analysis.

dmx is a curated catalog of [marimo](https://marimo.io) notebooks for cancer-dependency analysis, plus a thin skill that lets an agent compose new analyses from them.
Each notebook is both a runnable demonstration and a source of pure functions other notebooks can [import and reuse](https://docs.marimo.io/guides/reusing_functions/) directly.
Given a new cancer-dependency question, the agent picks relevant notebooks, composes their functions into a new notebook, executes it in a live kernel, and hands back a self-contained, re-runnable result.

Breadbox is the REST API underlying the DepMap Portal: find datasets, search genes and compounds, retrieve matrix values, define model contexts, run two-class comparisons, and inspect precomputed associations. By default dmx targets `https://depmap.org/portal/breadbox`, which does not require an API key for these read-only examples.

## The catalog

Each notebook ships with a committed session snapshot under [`notebooks/__marimo__/session/`](notebooks/__marimo__/session/) so the molab preview renders cell outputs without re-executing.

| Notebook | Role | Preview |
|---|---|---|
| [`nb01_orientation.py`](notebooks/nb01_orientation.py) | Landing page: what Breadbox exposes and what questions this catalog should answer first | [![Open in molab](https://marimo.io/molab-shield.svg)](https://molab.marimo.io/github/broadinstitute/dmx/blob/main/notebooks/nb01_orientation.py) |
| [`nb02_dataset_discovery.py`](notebooks/nb02_dataset_discovery.py) | Breadbox primitives: `bb_get`, `bb_post`, dataset discovery, dimension search, feature listing | [![Open in molab](https://marimo.io/molab-shield.svg)](https://molab.marimo.io/github/broadinstitute/dmx/blob/main/notebooks/nb02_dataset_discovery.py) |
| [`nb03_gene_dependency_profile.py`](notebooks/nb03_gene_dependency_profile.py) | Gene -> Chronos dependency profile across cell lines, joined to model metadata | [![Open in molab](https://marimo.io/molab-shield.svg)](https://molab.marimo.io/github/broadinstitute/dmx/blob/main/notebooks/nb03_gene_dependency_profile.py) |
| [`nb04_context_comparison.py`](notebooks/nb04_context_comparison.py) | Lineage / mutation-defined context -> compare dependency distributions inside vs outside | [![Open in molab](https://marimo.io/molab-shield.svg)](https://molab.marimo.io/github/broadinstitute/dmx/blob/main/notebooks/nb04_context_comparison.py) |
| [`nb05_association_query.py`](notebooks/nb05_association_query.py) | Gene dependency or drug sensitivity slice -> precomputed Breadbox association hits across modalities | [![Open in molab](https://marimo.io/molab-shield.svg)](https://molab.marimo.io/github/broadinstitute/dmx/blob/main/notebooks/nb05_association_query.py) |

The agent-facing catalog table in `.claude/skills/compose-notebook/SKILL.md` is the detailed contract: it lists reusable helpers, import patterns, and current gotchas.

Related public catalogs of the same pattern: [jx](https://github.com/broadinstitute/jx) for JUMP Cell Painting, [fgx](https://github.com/broadinstitute/fgx) for FinnGenie human genetics, and [prx](https://github.com/broadinstitute/prx) for PROSPECT chemical genetics.

## Getting started

Clone this repo, open Claude Code inside it, and ask: *help me get started*.
The `getting-started` skill installs prereqs ([uv](https://docs.astral.sh/uv/) and the [marimo-pair](https://github.com/marimo-team/marimo-pair) skill), launches `nb01_orientation` in a live marimo kernel, and hands off to the `compose-notebook` skill for the actual analysis.

If you prefer to run setup by hand:

```bash
uv --version  # or: curl -LsSf https://astral.sh/uv/install.sh | sh
AGENT=claude-code  # or: codex
npx skills add marimo-team/marimo-pair -g --agent "$AGENT" -y
uvx marimo edit --sandbox notebooks/nb01_orientation.py
```

## License

BSD 3-Clause — see [LICENSE](LICENSE).
