# Generic functions in GECKO and geckopy that belong in RAVEN

Functions that need no enzyme concept live in RAVEN (MATLAB) and raven-toolbox (Python), and
GECKO and geckopy import them. This record lists each candidate, the decision, and where the
MATLAB and Python sides stand. Both sides get the same function in the same release, and the
ledgers (`ledgers/raven.yml`, `ledgers/gecko.yml`) carry the rows once the functions exist,
because `parity check` fails on a row that names a function that is not there.

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
| `apply_flux_data_constraints` / `constrainFluxData` | `raven_toolbox.experiments`, `applyFluxDataConstraints` | done | in progress | Takes the biomass reaction and carbon source as arguments. GECKO and geckopy keep a wrapper that reads them from the adapter. |
| `solver_support.py` (`optimize_model`, scipy and HiGHS repair, shadow prices) | `raven_toolbox.utils.solver` | in progress | not needed | MATLAB uses `solveLP`. |
| `calculate_mw`, `find_met_smiles` | raven-toolbox | in progress | to decide | Checked against what raven-toolbox already has. |
| Maintenance fit (Huber loss, study weights, bootstrap, uptake grid, NGAM from retentostat) | `raven_toolbox.biomass` | in progress | planned | Reuses `set_gam` and `fit_parameters`. yeast-GEM keeps data, ids and a thin wrapper. |
| `parseGECKOargs` | `parseRAVENargs` | not needed | in progress | RAVEN's function accepts every call form that GECKO's does. GECKO's is removed. |

## Stays in GECKO and geckopy

| Function | Reason |
|---|---|
| `ec_fva`, `pfba_enzymes`, `get_enzyme_bottlenecks` | Built on enzyme usage reactions and enzyme constraints. cobrapy's `pfba` covers the generic case. |
| KEGG per-organism REST download, UniProt download and loader | Per-organism protein data, different from RAVEN's bulk KEGG dump, and tied to the adapter. |
| `findGECKOroot` | Needed at install time, before RAVEN is on the path. |
| `buildGECKOdoc` | Legacy m2html build. The documentation is built from `gecko-docs`. Candidate for removal. |
| Evotune CMA-ES search | The search is generic, its parameters are per-kcat (mask, sigma, tying, source groups). |
| `gather_kcats`, `get_enzyme_data`, `kcat_tuning`, `ec_model`, the rest of `limit_proteins`, the adapter, `upgrade.py`, OKP, DLKcat and BRENDA readers | Use kcats, enzymes or the protein pool. |

## Ledger rows to add when the functions exist

* `raven.yml`: `loadFluxData`, `loadExperimentData`, `simulateConditions`,
  `computeCarbonWeights`, `experimentDistance`, `datasetRmse`, `applyFluxDataConstraints`, each
  `parity` with the Python name under `raven_toolbox.experiments`.
* `gecko.yml`: `fillCarbonNum` is replaced by `computeCarbonWeights`. `evotuneScore` and
  `loadEvotuneData` are wrappers around the RAVEN functions, with the Python names of the
  wrappers in geckopy. `constrainFluxData` is a wrapper of `applyFluxDataConstraints`.
* `parseGECKOargs` is removed from `gecko.yml`.

## Parity runs to repeat after the move

The `sensitivity_tuning_ectestgem` and `flexibilize_concs_ectestgem` scenarios call
`loadFluxData` and `constrainFluxData`. They match exactly in every mode today, so they are
rerun on the moved functions to show the match is kept.

## Open points

* geckopy's pin on raven-toolbox has to name the first release with `experiments` and the
  moved solver module.
* The seeded evotune run against the stored RMSE trace (`~/ecmodel-runs/b2/full`, seed 0) has
  to give the same result on the moved code.
* GECKO needs a RAVEN version requirement for the new functions.
