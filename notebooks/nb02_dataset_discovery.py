# /// script
# requires-python = ">=3.12,<3.14"
# dependencies = [
#     "marimo",
#     "polars==1.40.1",
#     "requests==2.32.5",
# ]
# ///

import marimo

__generated_with = "0.23.10"
app = marimo.App(width="medium")

with app.setup:
    import time

    import marimo as mo
    import polars as pl
    import requests

    BASE_URL = "https://depmap.org/portal/breadbox"
    # The portal nginx 403s the default python-requests User-Agent, so every call
    # sends an identifying one. Breadbox also intermittently 504s (gateway
    # timeout); bb_get / bb_post retry on 5xx and connection errors.
    USER_AGENT = "dmx/0.1 (+https://github.com/broadinstitute/dmx)"
    MAX_RETRIES = 6
    RETRY_BACKOFF = 3.0


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    # nb02 - dataset discovery

    Breadbox is easiest to compose against when the agent can discover the
    available datasets and dimensions first. This notebook owns the shared HTTP
    helpers for the dmx catalog and demonstrates the first discovery calls.

    The helpers are hardened against two portal quirks: a non-default
    `User-Agent` (the default is 403'd) and transient `5xx` gateway timeouts
    (retried with backoff).
    """)
    return


@app.function
def breadbox_headers() -> dict[str, str]:
    """Return headers for Breadbox requests (identifying User-Agent required)."""
    return {
        "Accept": "application/json",
        "User-Agent": USER_AGENT,
    }


@app.function
def breadbox_url(endpoint: str) -> str:
    """Return the full Breadbox URL for an endpoint."""
    return f"{BASE_URL}/{endpoint.lstrip('/')}"


@app.function
def bb_request(method: str, endpoint: str, *, params: dict | None = None, json: dict | None = None) -> object:
    """Call Breadbox with retries on 5xx and connection errors, return decoded JSON.

    The DepMap portal occasionally returns 504 Gateway Time-out on otherwise
    valid reads, so retry the same request a few times with linear backoff. 4xx
    is a real client error and is raised immediately.
    """
    url = breadbox_url(endpoint)
    timeout = 120 if method == "POST" else 60
    last_err: str | None = None
    for attempt in range(MAX_RETRIES):
        try:
            response = requests.request(
                method,
                url,
                params=params,
                json=json,
                headers=breadbox_headers(),
                timeout=timeout,
            )
            if response.status_code < 500:
                response.raise_for_status()
                return response.json()
            last_err = f"HTTP {response.status_code}"
        except requests.RequestException as exc:
            last_err = str(exc)
        time.sleep(RETRY_BACKOFF * (attempt + 1))
    raise RuntimeError(f"Breadbox {method} {endpoint} failed after {MAX_RETRIES} retries: {last_err}")


@app.function
def bb_get(endpoint: str, params: dict | None = None) -> object:
    """GET from Breadbox and return decoded JSON."""
    return bb_request("GET", endpoint, params=params)


@app.function
def bb_post(endpoint: str, body: dict) -> object:
    """POST JSON to Breadbox and return decoded JSON."""
    return bb_request("POST", endpoint, json=body)


@app.function
def records_frame(records: object) -> pl.DataFrame:
    """Convert a list of JSON records into a polars DataFrame."""
    if records is None:
        return pl.DataFrame()
    if isinstance(records, list):
        return pl.DataFrame(records, infer_schema_length=10_000) if records else pl.DataFrame()
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
def search_dimensions(substring: str, type_name: str | None = None, limit: int = 10) -> pl.DataFrame:
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


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ## To extend

    - Pull one gene's dependency profile from a `Chronos_Combined`-style matrix dataset and join to model metadata (nb03).
    - Define a lineage or mutation context with `temp/context` and compare dependency inside vs outside (nb04).
    - Query `temp/associations/query-slice` for a dependency slice to surface correlated features (nb05).
    """)
    return


if __name__ == "__main__":
    app.run()
