# /// script
# requires-python = ">=3.12,<3.14"
# dependencies = [
#     "marimo",
#     "polars==1.40.1",
#     "requests==2.32.5",
# ]
# ///

import marimo

__generated_with = "0.23.5"
app = marimo.App(width="medium")

with app.setup:
    import marimo as mo
    import polars as pl
    import requests

    BASE_URL = "https://depmap.org/portal/breadbox"


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    # nb02 - dataset discovery

    Breadbox is easiest to compose against when the agent can discover the
    available datasets and dimensions first. This notebook owns the shared
    HTTP helpers for the dmx catalog and demonstrates the first discovery
    calls.
    """)
    return


@app.function
def breadbox_headers() -> dict[str, str]:
    """Return headers for Breadbox requests."""
    return {
        "Accept": "application/json",
        "User-Agent": "dmx/0.1 (+https://github.com/broadinstitute/dmx)",
    }


@app.function
def breadbox_url(endpoint: str) -> str:
    """Return the full Breadbox URL for an endpoint."""
    return f"{BASE_URL}/{endpoint.lstrip('/')}"


@app.function
def bb_get(endpoint: str, params: dict | None = None) -> object:
    """GET from Breadbox and return decoded JSON."""
    response = requests.get(
        breadbox_url(endpoint),
        params=params,
        headers=breadbox_headers(),
        timeout=60,
    )
    response.raise_for_status()
    return response.json()


@app.function
def bb_post(endpoint: str, body: dict) -> object:
    """POST JSON to Breadbox and return decoded JSON."""
    response = requests.post(
        breadbox_url(endpoint),
        json=body,
        headers=breadbox_headers(),
        timeout=120,
    )
    response.raise_for_status()
    return response.json()


@app.function
def records_frame(records: object) -> pl.DataFrame:
    """Convert a list of JSON records into a polars DataFrame."""
    if records is None:
        return pl.DataFrame()
    if isinstance(records, list):
        return (
            pl.DataFrame(records, infer_schema_length=10_000)
            if records
            else pl.DataFrame()
        )
    if isinstance(records, dict):
        return pl.DataFrame([records], infer_schema_length=10_000)
    return pl.DataFrame({"value": [records]})


@app.function
def list_data_types() -> pl.DataFrame:
    """List Breadbox data types derived from the dataset catalog."""
    datasets = list_datasets()
    if "data_type" not in datasets.columns:
        return pl.DataFrame()
    return datasets.select("data_type").drop_nulls().unique().sort("data_type")


@app.function
def list_dimension_types() -> pl.DataFrame:
    """List Breadbox dimension types."""
    return records_frame(bb_get("types/dimensions"))


@app.function
def list_datasets(**filters: object) -> pl.DataFrame:
    """List Breadbox datasets, optionally filtered by feature/sample type."""
    params = {key: value for key, value in filters.items() if value is not None}
    return records_frame(bb_get("datasets/", params=params))


@app.function
def search_dimensions(
    substring: str, type_name: str | None = None, limit: int = 10
) -> pl.DataFrame:
    """Search genes, compounds, or cell lines by substring."""
    params: dict[str, object] = {"substring": substring, "limit": limit}
    if type_name is not None:
        params["type_name"] = type_name
    return records_frame(bb_get("datasets/dimensions/", params=params))


@app.function
def dataset_features(dataset_id: str) -> pl.DataFrame:
    """List features available in a matrix dataset."""
    return records_frame(bb_get(f"datasets/features/{dataset_id}"))


@app.cell
def _():
    data_types = list_data_types()
    mo.md("## Data types")
    mo.ui.table(data_types, page_size=10)
    return


@app.cell
def _():
    dimension_types = list_dimension_types()
    mo.md("## Dimension types")
    mo.ui.table(dimension_types, page_size=10)
    return


@app.cell
def _():
    gene_model_datasets = list_datasets(feature_type="gene", sample_type="depmap_model")
    mo.md("## Gene x DepMap model datasets")
    mo.ui.table(gene_model_datasets, page_size=12)
    return


@app.cell
def _():
    kras_hits = search_dimensions("KRAS", type_name="gene", limit=5)
    mo.md("## Search example: KRAS")
    mo.ui.table(kras_hits, page_size=5)
    return


if __name__ == "__main__":
    app.run()
