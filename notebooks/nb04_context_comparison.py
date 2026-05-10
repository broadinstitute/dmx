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

__generated_with = "0.23.5"
app = marimo.App(width="medium")

with app.setup:
    import sys
    from pathlib import Path

    import altair as alt
    import marimo as mo
    import polars as pl

    NOTEBOOK_DIR = Path(__file__).parent
    if str(NOTEBOOK_DIR) not in sys.path:
        sys.path.insert(0, str(NOTEBOOK_DIR))

    from nb02_dataset_discovery import bb_post
    from nb03_gene_dependency_profile import (
        DEFAULT_DEPENDENCY_DATASET,
        matrix_feature,
        tabular_metadata,
    )

    DEFAULT_GENE = "KRAS"
    DEFAULT_LINEAGE = "Lung"


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    # nb04 - context comparison

    Breadbox context expressions define subsets of DepMap models using
    metadata, mutation calls, expression, or other matrix slices. This
    notebook starts with a simple lineage context and compares a gene's
    dependency inside vs outside the context.
    """)
    return


@app.function
def lineage_context(lineage: str) -> dict:
    """Breadbox context expression for a single Oncotree lineage."""
    return {
        "name": lineage,
        "dimension_type": "depmap_model",
        "expr": {"==": [{"var": "lineage"}, lineage]},
        "vars": {
            "lineage": {
                "dataset_id": "depmap_model_metadata",
                "identifier": "OncotreeLineage",
                "identifier_type": "column",
            }
        },
    }


@app.function
def mutation_lineage_context(gene: str, lineage: str) -> dict:
    """Breadbox context expression for mutant models within one lineage."""
    return {
        "name": f"{gene}-mutant {lineage}",
        "dimension_type": "depmap_model",
        "expr": {
            "and": [
                {"==": [{"var": "lineage"}, lineage]},
                {">": [{"var": "mutation"}, 0]},
            ]
        },
        "vars": {
            "lineage": {
                "dataset_id": "depmap_model_metadata",
                "identifier": "OncotreeLineage",
                "identifier_type": "column",
            },
            "mutation": {
                "dataset_id": "mutations_hotspot",
                "identifier": gene,
                "identifier_type": "feature_label",
            },
        },
    }


@app.function
def resolve_context(context: dict) -> dict:
    """Resolve a Breadbox context expression to model IDs."""
    return bb_post("temp/context", context)


@app.function
def dependency_context_table(
    gene: str,
    context: dict,
    dataset_id: str = DEFAULT_DEPENDENCY_DATASET,
) -> pl.DataFrame:
    """Return dependency scores labeled by context membership."""
    resolved = resolve_context(context)
    in_ids = set(resolved.get("ids", []))
    scores = matrix_feature(dataset_id, gene, value_name="dependency")
    metadata = tabular_metadata(
        ["CellLineName", "OncotreeLineage", "OncotreePrimaryDisease"]
    )
    return (
        scores.join(metadata, on="depmap_id", how="left")
        .with_columns(pl.col("depmap_id").is_in(list(in_ids)).alias("in_context"))
        .sort(["in_context", "dependency"], descending=[True, False])
    )


@app.cell
def _():
    context = lineage_context(DEFAULT_LINEAGE)
    resolved_context = resolve_context(context)
    mo.md(
        f"## Context: {DEFAULT_LINEAGE} ({len(resolved_context.get('ids', []))} models)"
    )
    return (context,)


@app.cell
def _(context):
    context_table = dependency_context_table(DEFAULT_GENE, context)
    mo.ui.table(context_table.head(20), page_size=20)
    return (context_table,)


@app.cell
def _(context_table):
    context_summary = (
        context_table.drop_nulls("dependency")
        .group_by("in_context")
        .agg(
            pl.col("dependency").mean().alias("mean_dependency"),
            pl.col("dependency").median().alias("median_dependency"),
            pl.col("dependency").quantile(0.25).alias("q25_dependency"),
            pl.col("dependency").quantile(0.75).alias("q75_dependency"),
            pl.len().alias("n_models"),
        )
        .sort("in_context", descending=True)
    )
    mo.md(f"## {DEFAULT_GENE}: in-context vs out-of-context")
    mo.ui.table(context_summary, page_size=5)
    return


@app.cell
def _(context_table):
    _plot_frame = context_table.drop_nulls(["dependency", "in_context"]).with_columns(
        pl.col("in_context").cast(pl.String).alias("context")
    )
    _chart = (
        alt.Chart(_plot_frame)
        .mark_boxplot(extent="min-max")
        .encode(
            x=alt.X("context:N", title=f"in {DEFAULT_LINEAGE} context?"),
            y=alt.Y("dependency:Q", title=f"{DEFAULT_GENE} Chronos dependency"),
            color=alt.Color("context:N", legend=None),
        )
        .properties(width=360, height=360)
    )
    mo.ui.altair_chart(_chart)
    return


if __name__ == "__main__":
    app.run()
