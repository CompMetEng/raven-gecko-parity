# Generic functions in GECKO and geckopy that belong in RAVEN

Functions that need no enzyme concept live in RAVEN (MATLAB) and raven-toolbox (Python), and
GECKO and geckopy import them. This record lists each candidate, the decision, and where the
MATLAB and Python sides stand. Both sides get the same function in the same release, and the
ledgers (`ledgers/raven.yml`, `ledgers/gecko.yml`) carry the rows once the functions exist,
because `parity check` fails on a row that names a function that is not there.

The MATLAB functions of the experiments package are in `src/experiments` of RAVEN. In GECKO, `constrainFluxData`, `loadEvotuneData` and `evotuneScore` are wrappers of them, and `loadFluxData`, `fillCarbonNum` and `parseGECKOargs` are removed (`addCarbonNum` and `getcarbonnum` are stubs that raise `GECKO:removedFunction`).

Status values: `done` (implemented and tested on a branch), `in progress`, `planned`, `kept`
(stays where it is, with the reason).

## Moves to RAVEN

| Function (today) | Target | Python | MATLAB | Note |
|---|---|---|---|---|
| Flux data and experiment files: `load_flux_data`, `FluxData` / `loadFluxData` | `raven_toolbox.experiments`, `src/experiments` | done | in progress | The columns are `Condition Ptot grRate <met> (<rxn>) ... weight source`. |
| `load_experiment_data` / `loadEvotuneData` | `raven_toolbox.experiments`, `src/experiments` as `loadExperimentData` | done | in progress | Takes a folder. Files `fluxData.tsv`, `maxGrowth.tsv`, `zeroFlux.tsv`. |
| `simulate_conditions` / part of `evotuneScore` | `raven_toolbox.experiments`, `simulateConditions` | done | in progress | The anaerobic switch, the protein-content change and the solver call are function arguments. |
| `experiment_distance`, `dataset_rmse` / part of `evotuneScore` | `raven_toolbox.experiments`, `experimentDistance`, `datasetRmse` | done | in progress | Constants: biomass carbon 41, infeasible penalty 99. |
| `compute_carbon_weights` / `fillCarbonNum` (`model.excarbon`) | `raven_toolbox.experiments`, `computeCarbonWeights` | done | in progress | Returns a map keyed by reaction id and does not write a model field. |
| `apply_flux_data_constraints` / `constrainFluxData` | `raven_toolbox.experiments`, `applyFluxDataConstraints` | done | done | Takes the biomass reaction and carbon source as arguments. GECKO and geckopy keep a wrapper that reads them from the adapter. |
| `solver_support.py` (`optimize_model`, scipy and HiGHS repair, shadow prices) | `raven_toolbox.utils.solver` | done | not needed | MATLAB uses `solveLP`. |
| `calculate_mw`, `find_met_smiles` | raven-toolbox (`utils.protein`, `annotation.smiles`) | done | to decide | raven-toolbox had no equivalent. `find_met_smiles` needs `cache_path`. |
| Maintenance fit (Huber loss, study weights, bootstrap, uptake grid, NGAM from retentostat) | `raven_toolbox.biomass` | done (yeast-GEM draft PR 433) | planned | Reuses `set_gam` and `fit_parameters`. yeast-GEM keeps data, ids and a thin wrapper. |
| `parseGECKOargs` | `parseRAVENargs` | not needed | done | RAVEN's function accepts every call form that GECKO's does. GECKO's is removed. |

## Stays in GECKO and geckopy

| Function | Reason |
|---|---|
| `ec_fva`, `pfba_enzymes`, `get_enzyme_bottlenecks` | Built on enzyme usage reactions and enzyme constraints. cobrapy's `pfba` covers the generic case. |
| KEGG per-organism REST download, UniProt download and loader | Per-organism protein data, different from RAVEN's bulk KEGG dump, and tied to the adapter. |
| `findGECKOroot` | Needed at install time, before RAVEN is on the path. |
| `buildGECKOdoc` | Legacy m2html build. The documentation is built from `gecko-docs`. Candidate for removal. |
| Evotune CMA-ES search | The search is generic, its parameters are per-kcat (mask, sigma, tying, source groups). |
| `gather_kcats`, `get_enzyme_data`, `kcat_tuning`, `ec_model`, the rest of `limit_proteins`, the adapter, `upgrade.py`, OKP, DLKcat and BRENDA readers | Use kcats, enzymes or the protein pool. |

## Differences between the two APIs

* The MATLAB functions work on plain structs with the GECKO field names (`conds`, `Ptot`,
  `grRate`, `exchFluxes`, `exchMets`, `exchRxnIDs`, `weight`, `source`).
* `applyFluxDataConstraints(model, fluxData, bioRxnId, ...)` takes the biomass reaction as a
  positional argument, and its `'condition'` option counts from 1. In Python `bio_rxn_id` is
  a keyword and the condition counts from 0.
* `computeCarbonWeights` returns a `containers.Map` and `simulateConditions` a struct array.
  The hooks return the model, and `optimize` defaults to `solveLP`.
* A condition is feasible only when the solve is optimal, and a biomass reaction missing from
  the model raises an error.

## Ledger rows

`parity check` passes for the raven pair (two `unreviewed` warnings: `get_gam` and `fit_maintenance`, until it is decided whether RAVEN gets a MATLAB version) and `parity ledger-audit` is clean for both pairs.

* `raven.yml`: `loadFluxData`, `loadExperimentData`, `simulateConditions`, `computeCarbonWeights`, `datasetRmse`, `experimentDistance` and `applyFluxDataConstraints` are `parity` rows with the Python names under `raven_toolbox.experiments`. The result types and constants of the package, `find_met_smiles`, `met_smiles` and the maintenance functions are `python-only`, `internal` or `unreviewed` rows.
* `gecko.yml`: `fillCarbonNum`, `loadFluxData` (the row now names `geckopy.load_flux_data` as `internal`) and `parseGECKOargs` are gone. `evotuneScore` and `loadEvotuneData` are `subsumed` (wrappers of the RAVEN functions; geckopy calls `raven_toolbox.experiments` inside the tuning). `addCarbonNum` and `getcarbonnum` are `matlab-only` stubs. The seven evotune result types and constants of geckopy and the re-exports in `geckopy.databases` are updated. The notes that named `geckopy.solver_support` name `raven_toolbox.utils.solver`.
* The scenario entries `loadFluxData` are removed from the `gecko` scenarios, because the function belongs to the raven pair.

Eight findings of `parity check` for the gecko pair are not related to this move: `fillKcatsFromIsozymes`, `showAdapterParams`, `format_adapter_params`, `resolve_gr_exp`, `KcatEvidence`, `kcat_evidence`, `kcats_for_reactions` and `show_adapter_params` have no row in this branch of the ledger.

## Parity runs on the moved code

The scenarios that touch the moved code run Python (hash seeds 0, 1 and 2) against MATLAB, with RAVEN `feat/experiments` (e431f96a) and GECKO `refactor/parse-raven-args` (f84240c1) on the MATLAB side, and geckopy `feat/experiment-data` with raven-toolbox `feat/experiments` on the Python side. The results are in `/cephyr/NOBACKUP/groups/compmeteng/eduardk/pw/out/20261120` (`status/<id>.txt`, `logs/<id>.log`).

| Scenario | Result |
|---|---|
| `flexibilize_concs_ectestgem` | match (`loadFluxData` and `constrainFluxData`, loose and percentage mode) |
| `flux_data_edge_cases_ectestgem` | match: `loadFluxData` raises on both sides for a ragged row, a file with only the fixed columns and a missing file |
| `evotune_edge_cases_ectestgem` | match |
| `okp_kcat_ectestgem` | match |
| `kcat_io_edge_cases_ectestgem` | match |
| `molecular_weight_sequences` | match |
| `ec_analysis_light_humangem` | match |
| `ec_fva_ecyeastgem` | running (genome scale) |
| `ec_proteomics_ecyeastgem` | differs in `flexibilize` only (75 values; Python seeds 0, 1 and 2 agree with each other). The parity plan lists the optimal vertex under the proteomics constraints as a documented difference of this scenario. The flux data steps match. |

The scenarios changed for the move: the MATLAB side calls `loadFluxData(path)`, the Python side imports from `raven_toolbox.experiments`, `raven_toolbox.annotation` and `raven_toolbox.utils.solver`, and the evotune scenario uses the experiment-data names. Two further changes in `evotune_edge_cases_ectestgem` follow other changes in the toolboxes: geckopy takes `sigma0_log` inside each source group, and GECKO takes the biomass reaction from `modelAdapter.params.bioRxn` and the group of each source from `kcatGroups`, so the scenario sets both.

## Open points

* geckopy's pin on raven-toolbox has to name the first release with `experiments` and the
  moved solver module.
* `ec_fva_ecyeastgem` has to finish; the table above gets its result.
* `ec_proteomics_ecyeastgem` has no run on the code from before the move, so that the `flexibilize` difference is known to be the documented one and not a new one.
* The seeded evotune run against the stored RMSE trace (`~/ecmodel-runs/b2/full`, seed 0) has
  to give the same result on the moved code.
* GECKO needs a RAVEN version requirement for the new functions.

## Pull requests

| Repository | Pull request | Waits for |
|---|---|---|
| raven-toolbox | [#204](https://github.com/CompMetEng/raven-toolbox/pull/204) (draft) | |
| RAVEN | [#28](https://github.com/CompMetEng/RAVEN/pull/28) (draft) | |
| geckopy | [#98](https://github.com/CompMetEng/geckopy/pull/98) (draft) | the raven-toolbox release, and a raised `raven-toolbox` requirement |
| GECKO | [#38](https://github.com/CompMetEng/GECKO/pull/38) (draft) | the RAVEN release with `src/experiments` |
| yeast-GEM | [#433](https://github.com/SysBioChalmers/yeast-GEM/pull/433) (draft) | the raven-toolbox release |

