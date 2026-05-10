# dmx -- DepMap eXplore

A marimo notebook catalog for agent-composable analysis against the public [DepMap](https://depmap.org/portal/) Breadbox API.

The pattern mirrors `jx`, `fgx`, and `prx`: a small catalog of numbered marimo notebooks, each a runnable demonstration and a source of pure helpers that later notebooks can import. The data surface is Breadbox, the REST API underlying the DepMap Portal.

## The catalog

- `notebooks/nb01_orientation.py` - landing page: what Breadbox exposes and what questions this catalog should answer first
- `notebooks/nb02_dataset_discovery.py` - Breadbox primitives: `bb_get`, `bb_post`, dataset discovery, dimension search, and feature listing
- `notebooks/nb03_gene_dependency_profile.py` - gene -> Chronos dependency profile across cell lines, joined to model metadata
- `notebooks/nb04_context_comparison.py` - lineage / mutation-defined context -> compare dependency distributions inside vs outside the context
- `notebooks/nb05_association_query.py` - gene dependency or drug sensitivity slice -> precomputed Breadbox association hits across modalities

The agent-facing catalog table in `.claude/skills/compose-notebook/SKILL.md` is the detailed contract: it lists reusable helpers, import patterns, and current gotchas.

## Getting started

Clone this repo, open Codex inside it, and ask: *help me get started*. The `getting-started` skill launches `nb01_orientation` in a live marimo kernel and hands off to `compose-notebook` for analysis.

If you prefer to run setup by hand:

```bash
uv --version
cp .env.example .env  # optional; public Breadbox reads usually work without auth
npx skills add marimo-team/skills -g --agent codex -y
npx skills add marimo-team/marimo-pair -g --agent codex -y
uvx marimo edit --sandbox notebooks/nb01_orientation.py
```

## What this is for

DepMap combines genome-scale perturbation screens, omics, drug sensitivity, and cell-line metadata. Breadbox gives a stable programmatic API for common moves: find datasets, search genes and compounds, retrieve matrix values, define model contexts, run two-class comparisons, and inspect precomputed associations. By default dmx targets the public Breadbox API at `https://depmap.org/portal/breadbox`; override with `BREADBOX_BASE_URL` if needed. `DEPMAP_MCP_TOKEN` and `API_PROXY_PASSWORD` are optional fallbacks for authenticated gateways.

dmx is the catalog version of those moves. The first goal is simple but useful: given a gene, find cancer contexts where its dependency is selective, annotate the relevant cell lines, and surface correlated genomic or drug features as follow-up hypotheses.

## License

BSD 3-Clause - see [LICENSE](LICENSE).
