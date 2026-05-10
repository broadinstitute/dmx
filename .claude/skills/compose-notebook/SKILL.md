---
name: compose-notebook
description: >-
  Compose a new marimo notebook in the dmx repo by reusing @app.function
  helpers from the existing notebook catalog to answer a DepMap Breadbox
  question end-to-end. Trigger for notebooks, analyses, figures, or vignettes
  involving DepMap CRISPR dependency, expression, mutations, drug sensitivity,
  cell-line metadata, Breadbox context expressions, two-class comparisons, or
  precomputed associations.
---

# Compose a new marimo notebook from the dmx catalog

## What this skill is for

dmx is a catalog of marimo notebooks against the public Breadbox REST API at
`https://depmap.org/portal/breadbox`. There is no Python SDK and no MCP server
in the analysis path. Existing notebooks own the helpers they introduce, and
later notebooks import those helpers directly.

## The catalog at a glance

| Module | Reusable functions / globals | What they do |
|---|---|---|
| `nb01_orientation` | no helpers | Human/agent landing page: Breadbox surfaces, first questions, and catalog map. |
| `nb02_dataset_discovery` | `BASE_URL`, `bb_get(endpoint, params=None)`, `bb_post(endpoint, body)`, `records_frame(records)`, `list_data_types()`, `list_dimension_types()`, `list_datasets(**filters)`, `search_dimensions(substring, type_name=None, limit=10)`, `dataset_features(dataset_id)` | Breadbox core primitives. Start here for "what data exists?", gene/compound/model search, and dataset IDs. |
| `nb03_gene_dependency_profile` | `tabular_metadata(columns)`, `matrix_feature(dataset_id, feature_label)`, `gene_dependency_profile(gene, dataset_id, metadata_columns)` | Gene -> Chronos dependency profile across DepMap models, joined to cell-line metadata. |
| `nb04_context_comparison` | `lineage_context(lineage)`, `mutation_lineage_context(gene, lineage)`, `resolve_context(context)`, `dependency_context_table(gene, context)` | Build a Breadbox context, fetch dependency scores, and label in-context vs out-of-context models. |
| `nb05_association_query` | `association_query(dataset_id, identifier, identifier_type)`, `association_table(response)`, `top_associations(dataset_id, identifier, identifier_type, n)` | Query Breadbox precomputed associations for a gene dependency or drug sensitivity slice. |

When the question is not obviously answered by a row, read the closest
notebook before writing new code.

## Cross-notebook import recipe

Use plain Python imports from the `notebooks/` directory:

```python
with app.setup:
    import sys
    from pathlib import Path

    NOTEBOOK_DIR = Path(__file__).parent
    if str(NOTEBOOK_DIR) not in sys.path:
        sys.path.insert(0, str(NOTEBOOK_DIR))

    from nb02_dataset_discovery import bb_get, bb_post, records_frame
```

Match transitive dependencies in the PEP 723 header. Importing a helper from
another notebook executes that notebook's setup block.

## Process for a new composition

1. Map the user question to existing helpers.
2. Validate the dataset ID and identifiers early with `search_dimensions` and
   `dataset_features`.
3. Create `notebooks/nbNN_<topic>.py` with a PEP 723 header.
4. Reuse helpers by import; do not duplicate Breadbox plumbing.
5. Keep expensive API calls behind clear cells and summarize large responses.
6. Run the notebook in a live marimo sandbox kernel and inspect headline
   tables/plots before reporting completion.
7. Add a row to this catalog table if the new notebook introduces stable
   helpers.

## Gotchas

- Matrix endpoints return column-oriented mappings: `{feature: {sample: value}}`.
  Convert deliberately before joining metadata.
- Tabular endpoints are also column-oriented. Build a `depmap_id` column from
  the union of column keys before joining.
- Association results include multiple modalities. Filter by
  `other_dataset_given_id` when the user wants expression-only, CRISPR-only,
  mutation-only, or drug-only hits.
- `log10qvalue` can include sentinel values for self-hits. Cap or filter
  before plotting.
- Context expressions use JSON-logic syntax. Start with simple lineage or
  mutation+lineage contexts before composing more elaborate filters.
