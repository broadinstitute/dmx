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

    from nb02_dataset_discovery import bb_post, records_frame  # noqa: E402

    DEFAULT_DATASET = "Chronos_Combined"
    DEFAULT_IDENTIFIER = "KRAS"
    DEFAULT_IDENTIFIER_TYPE = "feature_label"


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    # nb05 - association query

    Breadbox exposes precomputed associations for a dataset slice. For a gene
    dependency slice, this answers: what expression, mutation, dependency, drug
    sensitivity, or subtype features correlate with the dependency profile?
    """)
    return


@app.function
def association_query(dataset_id: str, identifier: str, identifier_type: str = "feature_label") -> dict:
    """Query Breadbox precomputed associations for one dataset slice."""
    return bb_post(
        "temp/associations/query-slice",
        {
            "slice_query": {
                "dataset_id": dataset_id,
                "identifier": identifier,
                "identifier_type": identifier_type,
            }
        },
    )


@app.function
def association_table(response: dict) -> pl.DataFrame:
    """Convert a Breadbox association response to a sorted DataFrame."""
    associations = response.get("associated_dimensions", [])
    df = records_frame(associations)
    if df.is_empty():
        return df
    sort_col = "log10qvalue" if "log10qvalue" in df.columns else df.columns[0]
    return df.sort(sort_col)


@app.function
def top_associations(
    dataset_id: str,
    identifier: str,
    identifier_type: str = "feature_label",
    n: int = 50,
) -> pl.DataFrame:
    """Return the top association rows for a dataset slice."""
    response = association_query(dataset_id, identifier, identifier_type)
    return association_table(response).head(n)


@app.cell
def _():
    associations = top_associations(DEFAULT_DATASET, DEFAULT_IDENTIFIER, DEFAULT_IDENTIFIER_TYPE, n=50)
    mo.md(f"## Top associations: {DEFAULT_IDENTIFIER} in {DEFAULT_DATASET}")
    mo.ui.table(associations, page_size=15)
    return (associations,)


@app.cell
def _(associations):
    if "other_dataset_given_id" in associations.columns:
        modality_counts = associations.group_by("other_dataset_given_id").len().sort("len", descending=True)
    else:
        modality_counts = pl.DataFrame()
    mo.md("## Returned modalities")
    mo.ui.table(modality_counts, page_size=10)
    return


@app.cell
def _(associations):
    required = {
        "other_dimension_label",
        "correlation",
        "log10qvalue",
        "other_dataset_given_id",
    }
    if required.issubset(set(associations.columns)):
        plot_frame = associations.head(20).with_columns(pl.col("log10qvalue").clip(-50, 0).alias("capped_log10qvalue"))
        chart_output = mo.ui.altair_chart(
            alt.Chart(plot_frame)
            .mark_bar()
            .encode(
                x=alt.X("capped_log10qvalue:Q", title="log10 q-value (capped)"),
                y=alt.Y("other_dimension_label:N", sort="x", title=None),
                color=alt.Color("other_dataset_given_id:N", title="dataset"),
                tooltip=[
                    "other_dimension_label",
                    "other_dataset_given_id",
                    "correlation",
                    "log10qvalue",
                ],
            )
            .properties(width=620, height=460)
        )
    else:
        chart_output = mo.md("Association response did not include the expected plotting columns.")
    chart_output
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ## To extend

    - Swap the slice for a drug-sensitivity dataset identifier to find genes correlated with a compound response.
    - Filter associations to a single modality (e.g. expression only) and rank by correlation sign.
    - Chain from nb03/nb04: take a selectively-essential gene and ask what co-dependencies travel with it.
    """)
    return


if __name__ == "__main__":
    app.run()
