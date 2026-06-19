# /// script
# requires-python = ">=3.12,<3.14"
# dependencies = [
#     "marimo",
#     "requests==2.32.5",
# ]
# ///

import marimo

__generated_with = "0.23.10"
app = marimo.App(width="medium")

with app.setup:
    import marimo as mo
    import requests

    # Breadbox is the public REST layer behind the DepMap Portal. These reads are
    # public and need no token. Later notebooks own the shared HTTP helpers; nb01
    # keeps one minimal GET so orientation is live, not just prose.
    # A non-default User-Agent is required: the portal's nginx 403s the bare
    # python-requests UA, so send an identifying one.
    BASE_URL = "https://depmap.org/portal/breadbox"
    REQUEST_TIMEOUT = 60
    USER_AGENT = "dmx/0.1 (+https://github.com/broadinstitute/dmx)"


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    # nb01 - orientation

    dmx is a marimo catalog for the public **DepMap Breadbox API**. Breadbox is the
    REST layer behind the [DepMap Portal](https://depmap.org/portal/): datasets,
    dimensions, matrix values, tabular metadata, context expressions, two-class
    comparisons, and precomputed associations - genome-scale perturbation screens,
    omics, drug sensitivity, and cell-line metadata for cancer-dependency analysis.

    The catalog is intentionally small. Each notebook teaches one move and exposes
    the helpers it introduces as `@app.function`s so later notebooks can import
    them. By default dmx targets `https://depmap.org/portal/breadbox`, which needs
    no API key for these read-only examples.
    """)
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ## Catalog map

    | Notebook | Move |
    |---|---|
    | `nb02_dataset_discovery.py` | Breadbox HTTP helpers (`bb_get`, `bb_post`), dataset discovery, entity search |
    | `nb03_gene_dependency_profile.py` | Gene -> Chronos dependency profile across cell lines, joined to model metadata |
    | `nb04_context_comparison.py` | Context expression -> in/out dependency comparison |
    | `nb05_association_query.py` | Query precomputed associations for a dependency or drug slice |

    First composition target: given a gene, find cancer contexts where its
    dependency is selective, annotate the relevant cell lines, and surface
    correlated genomic or drug features.
    """)
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ## Breadbox surfaces this catalog starts with

    - Dataset discovery: `data_types/`, `datasets/`, `types/dimensions`
    - Entity search: `datasets/dimensions/`
    - Matrix retrieval: `datasets/matrix/{dataset_id}`
    - Tabular metadata: `datasets/tabular/{dataset_id}`
    - Context expressions: `temp/context`
    - Associations: `temp/associations/query-slice`

    Public reads do not require authentication. The notebooks call Breadbox
    directly with `requests`.
    """)
    return


@app.function
def reach_surface() -> list[dict]:
    """Hit Breadbox once and return its dimension types.

    GETs `types/dimensions`, the small enumeration that names the kinds of things
    Breadbox indexes (gene, compound, depmap_model, ...). It is the cheapest call
    that proves the surface is reachable; nb02 generalizes this into `bb_get`.
    """
    response = requests.get(
        f"{BASE_URL}/types/dimensions",
        headers={"Accept": "application/json", "User-Agent": USER_AGENT},
        timeout=REQUEST_TIMEOUT,
    )
    response.raise_for_status()
    return response.json()


@app.cell
def _():
    _dimension_types = reach_surface()
    _names = ", ".join(sorted(str(row.get("name", "?")) for row in _dimension_types))
    mo.md(
        f"""
        ## Live check

        Breadbox returned **{len(_dimension_types)}** dimension types: {_names}.
        These are the entity kinds the catalog searches, profiles, and compares.
        """
    )
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ## To extend

    - Discover which gene x DepMap-model matrix datasets exist (e.g. `Chronos_Combined`) and search for a gene of interest (nb02).
    - Pull one gene's dependency profile across cell lines and rank lineages by mean Chronos score (nb03).
    - Define a lineage or mutation context and compare a gene's dependency inside vs outside it (nb04).
    - Query precomputed associations for a dependency slice to surface correlated features across modalities (nb05).
    """)
    return


if __name__ == "__main__":
    app.run()
