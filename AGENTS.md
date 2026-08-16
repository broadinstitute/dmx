# AGENTS.md - dmx

Project-specific guidance for agents working in this catalog.
`README.md` is the human entry point.
This catalog uses the shared `vignette-catalog-compose-notebook` skill for setup, execution, and composition; its specifics live in `catalog.toml`.

## Skills (restore after clone)

The catalog skills are installed via `npx skills add`, recorded in the tracked `skills-lock.json`, but **not vendored** - the install stores (`.agents/`, `.claude/skills/*`) are gitignored.
A fresh clone has only the lock, so the on-disk skill content (and the `validate-notebook.sh` the rule below depends on) is missing until you restore it.

Do **not** use `npx skills update` for this.
`update` only refreshes agent stores that already exist on disk and has no `--agent` flag (only `-g`/`-p`/`-y`); on a fresh clone it materializes the Universal store (`.agents/`) but **not** `.claude/skills/`, where Claude Code discovers project skills - so the skills look installed yet Claude Code sees nothing.
Add them explicitly instead, from the repo root, with `add --agent` (the only command that targets a specific store):

```bash
npx skills@1.5.20 add carpenter-singh-lab/vignette-catalog-skills -s vignette-catalog-compose-notebook -s vignette-catalog-scaffold -a claude-code -a codex -y
npx skills@1.5.20 add marimo-team/skills -s marimo-notebook -a claude-code -a codex -y
npx skills@1.5.20 add marimo-team/marimo-pair -s marimo-pair -a claude-code -a codex -y
```

Repository instructions, `catalog.toml`, and the compose notebook contract override generic `marimo-notebook` advice when they are more specific.
This writes `.claude/skills/` and `.agents/`; both are gitignored.
Do it before relying on the skills or the validation rule, and run `/reload-skills` (or restart an already-open session) so they register.
(This instruction lives here, in a tracked file, on purpose: a skill cannot bootstrap its own install.)

## Launching notebooks

Always use `--sandbox` so PEP 723 inline metadata is provisioned:

```bash
uvx marimo edit --sandbox notebooks/nbNN_*.py
```

Do not improvise alternative launch commands.

## Validation rule

After composing or editing any notebook, run the `validate-notebook.sh` bundled with the installed `vignette-catalog-compose-notebook` skill with `--write` followed by the notebook path, then open it and look at the outputs.
Static checks do not catch wrong outputs, empty tables, stale endpoints, broken plots, or sign-convention mistakes.

## Architecture

- Catalog over library.
  Helpers are top-level `@app.function` cells in numbered notebooks; later notebooks import them via `sys.path`.
- Data surface is `rest`: this catalog reaches DepMap data over the public Breadbox REST API (`https://depmap.org/portal/breadbox`) via `requests`, no auth for these read-only examples.
- All HTTP goes through `bb_get` / `bb_post` in `nb02_dataset_discovery.py`.
  Import those rather than calling `requests` directly, so every call inherits the two hardening behaviors below.
- Do not add a Python package until repeated cross-notebook imports make it painful.

## Breadbox quirks (learned the hard way - keep `bb_get`/`bb_post` as the single choke point)

- **A non-default `User-Agent` is required.** A bare `requests` call (default `python-requests/x.y` UA) gets a `403 Forbidden` from the portal nginx, returning an HTML error page, not JSON.
  `breadbox_headers()` sends an identifying UA; do not drop it.
- **Breadbox 5xx's intermittently, and can go fully down.** Healthy reads still return transient `504 Gateway Time-out`, and during this catalog's first build the portal app tier was unreachable for ~an hour (edge up, app tier hanging - a portal outage, not a client or network problem; the separate Cloud Run DepMap MCP server stayed up).
  `bb_request` retries on `5xx` and connection errors with linear backoff; `4xx` is a real client error and raises immediately.
  If a whole composing session is timing out, check `https://depmap.org/portal/breadbox/types/dimensions` directly before assuming your code is wrong.
- **Dataset ids are stable enough to hardcode for examples** (`Chronos_Combined`, `depmap_model_metadata`, `mutations_hotspot`), but they can drift across DepMap quarterly releases.
  When a matrix/tabular/context call 404s, re-discover ids with `nb02`'s `list_datasets()` / `search_dimensions()` rather than guessing.
- **The DepMap release is pinned documentarily**, not in the request - Breadbox always serves the current release.
  The release the catalog was built against is recorded in `catalog.toml` `[data].version` (`26Q1`) and `nb01`'s `DEPMAP_RELEASE`; nb01 reads the live release off the `Chronos_Combined` dataset name and warns on drift.
  When DepMap ships a new quarter, bump both deliberately and re-run - the numbers re-anchor and the committed outputs no longer match.

## Conventions

Semantic line breaks in markdown.
ASCII-only.
Conventional Commits.
`ruff line-length = 120` is Python only.

## Shared contract (read before editing)

Read the installed `vignette-catalog-compose-notebook` skill and its notebook contract before authoring or editing a notebook.
Keep the shared workflow there and the DMX-specific rules here or in `catalog.toml`.

## When the question fits the catalog

The notebook-to-question routing lives in the `[[vignette]]` table in `catalog.toml` - each notebook, its helpers, and what it does - which is the single source the `vignette-catalog-compose-notebook` skill reads.
Do not mirror that table here; point at it.
