# /// script
# requires-python = ">=3.12,<3.14"
# dependencies = [
#     "marimo",
# ]
# ///

import marimo

__generated_with = "0.23.5"
app = marimo.App(width="medium")

with app.setup:
    import marimo as mo


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    # nb01 - orientation

    dmx is a marimo catalog for the public DepMap Breadbox API. Breadbox is the
    REST layer behind the DepMap Portal: datasets, dimensions, matrix values,
    tabular metadata, context expressions, two-class comparisons, and
    precomputed associations.

    The catalog is intentionally small. Each notebook teaches one move and
    exposes the helpers it introduces as `@app.function`s so later notebooks
    can import them.
    """)
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ## Catalog map

    | Notebook | Move |
    |---|---|
    | `nb02_dataset_discovery.py` | Breadbox HTTP helpers, dataset discovery, entity search |
    | `nb03_gene_dependency_profile.py` | Gene -> Chronos dependency profile across cell lines |
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


if __name__ == "__main__":
    app.run()
