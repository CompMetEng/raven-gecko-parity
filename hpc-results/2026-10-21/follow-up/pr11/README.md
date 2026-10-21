# RAVEN PR 11 and raven-toolbox PR 181

RAVEN `fix/10-yaml-booleans` at `875e45b` (base `79a7d8c`) and raven-toolbox `fix/parity-decisions` at `8fbaed4` (base `4128e2a`), MATLAB R2024b, GLPK, Vera. Both checkouts were restored to `develop` afterwards.

| Test class | Result |
| --- | --- |
| `tTasks` | 7 passed, 0 failed, 0 incomplete |
| `tIO` | 24 passed, 0 failed, 0 incomplete |
| `tManipulation` | 85 passed, 0 failed, 0 incomplete |
| `tSyntax` | 1 passed |
| `tINIT` | 31 passed, 0 failed, 0 incomplete |

`tTasks/parseTaskListReadsFlagsListsAndHeadersLeniently` and `tIO/writeYAMLmodelWritesShortestRoundTripFloats` both passed (rows in `matlab-tests/`).

Scenarios, MATLAB against Python seeds 0, 1 and 2, all `MATCH`: `remove_low_score_genes_smallyeast`, `replace_metabolite_toy`, `task_checking_smallyeast`, `task_list_parsing`, `yaml_roundtrip_smallyeast`. A grep of `scenarios/` for `parseTaskList`, `checkTasks`, `fitTasks`, `writeYAMLmodel`, `exportYAML` and `removeLowScoreGenes` found the first, third, fourth and fifth; `replace_metabolite_toy` was added by request.

`run.sh` is the job script. Per-test rows are in `matlab-tests/`, logs in `logs/`, verdicts in `status/`.
