# BRENDA files in the scenarios

GECKO and geckopy read the same BRENDA files from a database folder: `kcat.tsv`, `mw.tsv` and `sa.tsv`. Each file has a header row and explicit columns. EC codes are bare (`1.1.1.1`) and the organism column holds the plain name.

| File | Columns |
| --- | --- |
| `kcat.tsv` | `ec_code`, `substrate`, `organism`, `kcat_max`, `kcat_median`, `n`, `references` |
| `mw.tsv` | `ec_code`, `substrate`, `organism`, `mw`, `n`, `references` |
| `sa.tsv` | `ec_code`, `substrate`, `organism`, `sa_max`, `sa_median`, `n`, `references` |

GECKO resolves the folder through the adapter's `params.path` (the real release lives in `databases/brenda-<release>`), geckopy takes a folder for `load_brenda_data`. Specific activity is joined with molecular weight on EC code and organism, case-insensitively, and gives derived kcats in 1/s.

## Hand-built fixtures

`kcat_chain_ectestgem`, `fuzzy_kcat_coefficient_ectestgem` and `kcat_merge_edge_cases_ectestgem` each carry a small BRENDA folder in their `data/` folder. They share the following conventions:

- `mw.tsv` holds one inert row (EC `9.9.9.9`, organism `nonexistent organism`). MATLAB's `loadBRENDAdata` builds a `containers.Map` from `mw.tsv` and refuses an empty key list, so the file cannot hold the header only.
- `sa.tsv` holds the header only. The hand-built fixtures therefore do not exercise the specific-activity derived kcat path.
- A comment header on `mw.tsv` and `sa.tsv` is exactly one line. A longer header shifts the real header row far enough down that a loader misreads it as a data row.

The real release folder is read by `ec_model_full_yeastgem`, where the specific-activity join has data.
