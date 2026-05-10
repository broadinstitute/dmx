# dmx -- DepMap eXplore

A marimo notebook catalog for agent-composable analysis against the public [DepMap](https://depmap.org/portal/) Breadbox API.

The pattern mirrors [jx](https://github.com/broadinstitute/jx), [fgx](https://github.com/broadinstitute/fgx), and [prx](https://github.com/broadinstitute/prx): a small catalog of numbered marimo notebooks, each a runnable demonstration and a source of pure helpers that later notebooks can import.
The data surface is Breadbox, the REST API underlying the DepMap Portal.

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

## Getting started

Clone this repo, open Codex inside it, and ask: *help me get started*.
The `getting-started` skill launches `nb01_orientation` in a live marimo kernel and hands off to `compose-notebook` for analysis.

If you prefer to run setup by hand:

```bash
uv --version
AGENT=codex  # or: claude-code
npx skills add marimo-team/skills -g --agent "$AGENT" -y
npx skills add marimo-team/marimo-pair -g --agent "$AGENT" -y
uvx marimo edit --sandbox notebooks/nb01_orientation.py
```

## What this is for

DepMap combines genome-scale perturbation screens, omics, drug sensitivity, and cell-line metadata.
Breadbox gives a stable programmatic API for common moves: find datasets, search genes and compounds, retrieve matrix values, define model contexts, run two-class comparisons, and inspect precomputed associations.
By default dmx targets the public Breadbox API at `https://depmap.org/portal/breadbox`, which does not require an API key for these read-only examples.

dmx is the catalog version of those moves.
The first goal is simple but useful: given a gene, find cancer contexts where its dependency is selective, annotate the relevant cell lines, and surface correlated genomic or drug features as follow-up hypotheses.

Related public catalogs: [jx](https://github.com/broadinstitute/jx) for JUMP Cell Painting, [fgx](https://github.com/broadinstitute/fgx) for FinnGenie human genetics, and [prx](https://github.com/broadinstitute/prx) for PROSPECT chemical genetics.

## License

BSD 3-Clause - see [LICENSE](LICENSE).
