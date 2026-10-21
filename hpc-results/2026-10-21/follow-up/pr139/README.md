# Scenarios of raven-gecko-parity PR 139

Scenario files from `fix/min-nr-fluxes-solver` at `9326b3f`; RAVEN at `875e45b` (PR 11), raven-toolbox at `8fbaed4` (PR 181); MATLAB R2024b, GLPK, Vera. MATLAB against Python seeds 0, 1 and 2. Both checkouts were restored to `develop` afterwards.

| Scenario | Result |
| --- | --- |
| `remove_low_score_genes_ties_smallyeast` | `MATCH` (keeps `G2` for `PGI`, `C3`+`C4` for `PFK`, `H2` for `HXK`) |
| `task_list_edge_cases` | `MATCH` (4 tasks; `NaN` for the empty `CHANGED` bounds: `lb` of `R1` and `R2` in `T1`, `ub` of `R3` in `T2`, both bounds of `R4`) |
| `yaml_float_digits_smallyeast` | `MATCH` (bound lines `1.5e-07`, `2.5`, `0.30000000000000004`, `0.3333333333333333`, `-0.6666666666666666`, `1e+16`) |
| `task_list_parsing` | `MATCH` |
| `replace_metabolite_toy` | `MATCH` |

`run.sh` is the job script; `status/` holds the verdicts, `logs/` the compare output, `matlab/` and `python/` the result files.
