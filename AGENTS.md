# AGENTS.md - dmx

Project-specific guidance for agents working in this repository. This is the
public, runnable catalog of marimo notebooks for DepMap Breadbox analysis.
Planning and cross-instance coordination live in the primary
[`jx`](https://github.com/broadinstitute/jx) repo.

`README.md` is the human entry point. The skills under `.claude/skills/` are
the operational entry points: `getting-started` for first-run setup and
`compose-notebook` for adding a new analysis.

## Validation Rule

After composing or editing any notebook in `notebooks/`, launch it in a
marimo sandbox kernel and run all cells before reporting the task complete.
Static checks do not catch wrong endpoints, empty tables, stale API
assumptions, or broken plots.

Minimal launch:

```bash
PORT=$(python -c "import socket; s=socket.socket(); s.bind(('127.0.0.1',0)); print(s.getsockname()[1])")
env -u PYTHONPATH uvx marimo edit --sandbox --headless --no-token --port $PORT notebooks/nbNN_*.py
```

Then run static checks:

```bash
uvx ruff check notebooks/
uvx ruff format notebooks/
uvx marimo check notebooks/*.py
```

**Then, last, refresh the molab session snapshot** for any notebook whose source
changed in this task:

```bash
env -u PYTHONPATH uvx marimo export session --sandbox notebooks/nbNN_*.py
```

Order matters. Session snapshots store a `code_hash` per cell, and molab
attaches the stored output only when the snapshot hash matches the source
cell. Any later edit to the notebook source - including a `ruff format`
whitespace pass - shifts every `code_hash` and silently strips outputs in
the public molab preview. Always regenerate snapshots **after** the final
formatter / source edit, and commit the regenerated `.json` files in the
same change that touched the `.py` files.

## Architecture

- Catalog over library. Helpers live as `@app.function` cells in numbered
  notebooks. Later notebooks import from earlier notebooks by adding
  `notebooks/` to `sys.path`.
- Breadbox access is direct REST via `requests`; public read-only examples need
  no API key, and the existing DepMap MCP tools are thin wrappers around the
  same endpoints.
- Keep helpers close to API primitives: `requests`, `polars`, and small parsing
  functions.
- Raw API responses should be summarized before printing. Association and
  matrix endpoints can produce large payloads.
- Do not add a Python package until repeated cross-notebook imports make the
  notebook-as-library pattern painful.

## When the Question Fits the Catalog

Almost every first DepMap request should compose existing helpers:

- "What datasets exist?" -> `nb02_dataset_discovery`
- "Show KRAS dependency by lineage" -> `nb03_gene_dependency_profile`
- "Is gene X selective in KRAS-mutant lung?" -> `nb04_context_comparison`
- "What correlates with gene X dependency?" -> `nb05_association_query`

Read `.claude/skills/compose-notebook/SKILL.md` before writing new analysis
code.
