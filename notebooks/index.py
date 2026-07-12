# /// script
# requires-python = ">=3.12,<3.14"
# dependencies = [
#     "marimo",
#     "polars",
# ]
# ///

import marimo

__generated_with = "0.23.10"
app = marimo.App(width="medium")

with app.setup:
    import json
    from pathlib import Path

    import marimo as mo
    import polars as pl

    REPO = Path(__file__).resolve().parent.parent
    PROCESSED = REPO / "data" / "processed"
    REQUIRED_KEYS = {"description", "numbers", "files"}


@app.function
def discover_envelopes() -> list[dict]:
    """Glob data/processed/**/summary.json and keep the ones shaped like an envelope.

    Superset check on the required keys, so envelopes with extra keys (e.g. an
    optional `status`) are picked up too. The index computes nothing of its own -
    every number here is read straight from the node that produced it.
    """
    envelopes = []
    for path in sorted(PROCESSED.glob("**/summary.json")):
        try:
            data = json.loads(path.read_text())
        except (json.JSONDecodeError, OSError):
            continue
        if isinstance(data, dict) and REQUIRED_KEYS <= set(data):
            data["_analysis"] = path.parent.relative_to(PROCESSED).as_posix()
            data["_summary_path"] = path.relative_to(REPO).as_posix()
            envelopes.append(data)
    return envelopes


@app.function
def normalize_files(files: object) -> list[dict]:
    """A `files` entry is either a path string or a {path, caption} dict; return dicts."""
    rows = []
    for f in files or []:
        if isinstance(f, str):
            rows.append({"path": f, "caption": ""})
        elif isinstance(f, dict):
            rows.append({"path": str(f.get("path", "")), "caption": str(f.get("caption", ""))})
    return rows


@app.function
def status_badge(env: dict) -> str:
    """Render the optional status object as an inline badge; empty string when complete."""
    status = env.get("status")
    if not isinstance(status, dict):
        return ""
    state = status.get("state", "complete")
    if state == "complete":
        return ""
    note = status.get("note", "")
    return f"  `[{state}]`{f' - {note}' if note else ''}"


@app.cell(hide_code=True)
def _():
    mo.md(
        r"""
        # Catalog index

        Every `data/processed/**/summary.json` envelope this catalog has produced, one block
        each: its description, the numbers behind it, and the files it wrote. The index reads
        each number from the node that canonically produced it and computes nothing of its own.

        Dual-mode: in a checkout it discovers via the repo and writes `data/processed/index.{json,csv}`;
        standalone (molab) it globs whatever `data/processed/` is present and only renders.
        """
    )
    return


@app.cell
def _():
    envelopes = discover_envelopes()
    mo.md(
        f"**{len(envelopes)}** envelope(s) discovered under `data/processed/`."
        if envelopes
        else "No `summary.json` envelopes found under `data/processed/`. Run a notebook that writes one, then re-run this index."
    )
    return (envelopes,)


@app.cell
def _(envelopes):
    blocks = []
    for env in envelopes:
        numbers = env.get("numbers") or {}
        numbers_tbl = (
            pl.DataFrame({"metric": list(numbers.keys()), "value": [str(v) for v in numbers.values()]})
            if numbers
            else pl.DataFrame({"metric": [], "value": []})
        )
        files = normalize_files(env.get("files"))
        files_md = "\n".join(f"- `{r['path']}`{f' - {r['caption']}' if r['caption'] else ''}" for r in files) or "_(none)_"
        blocks.append(
            mo.vstack(
                [
                    mo.md(f"### `{env['_analysis']}`{status_badge(env)}"),
                    mo.md(env.get("description", "")),
                    mo.ui.table(numbers_tbl, selection=None, pagination=False),
                    mo.md(f"**Files** ({len(files)}):\n\n{files_md}"),
                ]
            )
        )
    mo.vstack(blocks) if blocks else mo.md("")
    return


@app.cell
def _(envelopes):
    # Collate one row per envelope and, when the checkout is writable, persist the artifact.
    index_rows = [
        {
            "analysis": env["_analysis"],
            "description": env.get("description", ""),
            "state": (env.get("status") or {}).get("state", "complete"),
            "n_numbers": len(env.get("numbers") or {}),
            "n_files": len(normalize_files(env.get("files"))),
            "summary_path": env["_summary_path"],
        }
        for env in envelopes
    ]
    index_df = pl.DataFrame(index_rows)

    written = []
    if index_rows:
        try:
            (PROCESSED / "index.json").write_text(json.dumps(index_rows, indent=2))
            index_df.write_csv(PROCESSED / "index.csv")
            written = ["data/processed/index.json", "data/processed/index.csv"]
        except OSError:
            written = []  # read-only datastore (e.g. molab): render only, no write

    mo.vstack(
        [
            mo.md("## Consolidated index"),
            mo.ui.table(index_df, selection=None, pagination=False),
            mo.md(
                "Wrote " + ", ".join(f"`{w}`" for w in written)
                if written
                else "_Render-only: the datastore is read-only, so no `index.{json,csv}` was written._"
            ),
        ]
    )
    return


if __name__ == "__main__":
    app.run()
