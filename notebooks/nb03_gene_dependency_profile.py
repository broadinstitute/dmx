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
    import sys
    from pathlib import Path

    import altair as alt
    import marimo as mo
    import polars as pl

    NOTEBOOK_DIR = Path(__file__).resolve().parent
    if str(NOTEBOOK_DIR) not in sys.path:
        sys.path.insert(0, str(NOTEBOOK_DIR))

    from nb02_dataset_discovery import bb_post  # noqa: E402

    DEFAULT_GENE = "KRAS"
    DEFAULT_DEPENDENCY_DATASET = "Chronos_Combined"
    DEFAULT_METADATA_COLUMNS = [
        "CellLineName",
        "OncotreeLineage",
        "OncotreePrimaryDisease",
    ]


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    # nb03 - gene dependency profile

    Pull one gene's dependency scores from a Breadbox matrix dataset and join
    them to cell-line metadata. This is the canonical DepMap first move:
    gene -> dependency profile -> lineages and diseases where the score is
    most negative.
    """)
    return


@app.function
def column_mapping_to_frame(raw: dict, index_name: str = "depmap_id") -> pl.DataFrame:
    """Convert a Breadbox column-oriented mapping to a row-oriented DataFrame."""
    ids: set[str] = set()
    for values in raw.values():
        if isinstance(values, dict):
            ids.update(values.keys())
    rows = []
    for item_id in sorted(ids):
        row = {index_name: item_id}
        for column, values in raw.items():
            if isinstance(values, dict):
                row[column] = values.get(item_id)
        rows.append(row)
    return pl.DataFrame(rows, infer_schema_length=10_000) if rows else pl.DataFrame()


@app.function
def tabular_metadata(columns: list[str]) -> pl.DataFrame:
    """Fetch cell-line metadata columns from Breadbox."""
    raw = bb_post("datasets/tabular/depmap_model_metadata", {"columns": columns})
    return column_mapping_to_frame(raw, index_name="depmap_id")


@app.function
def matrix_feature(dataset_id: str, feature_label: str, value_name: str = "value") -> pl.DataFrame:
    """Fetch one feature from a gene x model matrix dataset."""
    raw = bb_post(
        f"datasets/matrix/{dataset_id}",
        {"features": [feature_label], "feature_identifier": "label"},
    )
    values = raw.get(feature_label)
    if values is None and isinstance(raw, dict) and raw:
        values = next(iter(raw.values()))
    if not isinstance(values, dict):
        return pl.DataFrame()
    return pl.DataFrame(
        {
            "depmap_id": list(values.keys()),
            value_name: list(values.values()),
        }
    )


@app.function
def gene_dependency_profile(
    gene: str,
    dataset_id: str = DEFAULT_DEPENDENCY_DATASET,
    metadata_columns: list[str] = DEFAULT_METADATA_COLUMNS,
) -> pl.DataFrame:
    """Return dependency scores for one gene joined to model metadata."""
    scores = matrix_feature(dataset_id, gene, value_name="dependency")
    metadata = tabular_metadata(metadata_columns)
    return scores.join(metadata, on="depmap_id", how="left").sort("dependency")


@app.cell
def _():
    profile = gene_dependency_profile(DEFAULT_GENE)
    mo.md(f"## {DEFAULT_GENE} dependency profile")
    mo.ui.table(profile.head(20), page_size=20)
    return (profile,)


@app.cell
def _(profile):
    lineage_summary = (
        profile.drop_nulls(["OncotreeLineage", "dependency"])
        .group_by("OncotreeLineage")
        .agg(
            pl.col("dependency").mean().alias("mean_dependency"),
            pl.col("dependency").median().alias("median_dependency"),
            pl.len().alias("n_models"),
        )
        .filter(pl.col("n_models") >= 5)
        .sort("mean_dependency")
    )
    mo.md("## Lineage summary")
    mo.ui.table(lineage_summary.head(15), page_size=15)
    return (lineage_summary,)


@app.cell
def _(lineage_summary):
    _chart = (
        alt.Chart(lineage_summary.head(20))
        .mark_bar()
        .encode(
            x=alt.X("mean_dependency:Q", title="mean Chronos dependency"),
            y=alt.Y("OncotreeLineage:N", sort="x", title=None),
            tooltip=[
                "OncotreeLineage",
                "mean_dependency",
                "median_dependency",
                "n_models",
            ],
        )
        .properties(width=560, height=420)
    )
    mo.ui.altair_chart(_chart)
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ## To extend

    - Swap `DEFAULT_GENE` for another dependency (e.g. `BRAF`, `MYC`) and compare which lineages light up.
    - Restrict to a single lineage and rank the most dependent cell lines for a gene (feeds nb04).
    - Pull a second dataset (RNAi vs CRISPR) for the same gene and check concordance.
    """)
    return


if __name__ == "__main__":
    app.run()
