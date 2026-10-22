# HPC results 2026-10-22: merged develop, probes, seeds, solvers

Vera, MATLAB R2024b, GLPK unless a scenario names Gurobi. Gurobi 13.0.1 (token server; `versioninfo.license` is 2), gurobipy 13.0.3, Python 3.14.2, cobra 0.32.1.

Revisions (all `develop`): RAVEN `7bc599d` (contains `beadbc0`), raven-toolbox `12f04bf`, raven-gecko-parity `83d260e`, GECKO `52e46f7`, geckopy `277b74a` (main). `MANIFEST.txt` in `full-run/` records them. Its `dirty=68` for raven-gecko-parity counts untracked files (the per-scenario `.run.lock` files and the `results/` folder); `git status` lists no modified tracked file.

Batch layout: `tier-b/` (the new scenarios, below), `full-run/` (4.0), `probes/` (scripts, one JSON per probe, logs, the YAML files written by Python), `hashseeds/` (4), `python-hang/` (the `gap_report_toy` hang), `scip-tests/` (3.5).

## 4.0 Full run on merged develop

67 scenarios against Python seeds 0, 1 and 2: 61 `MATCH`, 4 `DIFFER`, 2 `ERROR`. 23 RAVEN test classes: 20 `PASS`, 3 `INCOMPLETE`, none failed. The first batch had 57 of 64 `MATCH`. The new `MATCH` are `replace_metabolite_toy` and `min_nr_fluxes_toy`; there are 3 new scenarios, all `MATCH`.

| Scenario | Verdict | Cause |
| --- | --- | --- |
| `model_edits_ectestgem`, `ec_fseof_ectestgem`, `ec_model_io_ectestgem` | `DIFFER` | the documented intended differences of the GECKO pair; not looked at again |
| `fuzzy_kcat_coefficient_ectestgem` | `DIFFER` | geckopy `main` is behind `develop`, as in the handover; not looked at again |
| `ec_model_full_yeastgem` | `ERROR` | `YeastGEMAdapter` was removed from GECKO; a rewrite against `TomlAdapter` exists locally and is not part of this batch |
| `gap_report_toy` | `ERROR` | the Python side hangs, see below. With cobra's FVA run in-process the scenario gives `MATCH` |

Test classes: `tReconstruction` has 3 skipped (`getKEGGModelForOrganismNeedsData`, `getModelFromKEGGNeedsData`, `getPhylDistNeedsData`, data missing), `tSolver` 1 (`optimizeProbSolvesMILPWithCobraGlpk`), `tUtils` 2 (`parallelWorkersRAVENTrueNoPCTReturnsZero`, `parallelWorkersRAVENTruePCTReturnsInf`, no Parallel Computing Toolbox cases). `tINIT` 31, `tGapfilling` 40 and `tTasks` 7 passed with Gurobi on the path (`GUROBI_OK=1` in the logs of all 23 classes).

### The `gap_report_toy` hang (Python side)

In `gap_report` (`raven_toolbox/gapfilling/report.py:106`, `_blocked_ids`) the call `cobra.flux_analysis.find_blocked_reactions` starts a multiprocessing pool and never returns. The stack dumps in `python-hang/` show the main thread in `Pool.join` inside `flux_variability_analysis` (and the pool's handler threads waiting for results). In the full run seed 1 stalled for 11 minutes at 0% CPU. Run again on a compute node with faulthandler set to 150 s, all 3 serial runs (seeds 1, 2, 3) and all 8 concurrent runs (seeds 0, 0, 0, 0, 4, 8, 12, 16) hung. Seed 0 had completed once in the full run, so the hang is not tied to a seed; it was reproducible on every later attempt. The Python start method here is the Python 3.14 default (forkserver). With `cobra.Configuration().processes = 1` set by a wrapper (`probes/venv-p1`), the scenario gives `MATCH` and 20 of 20 hash seeds match. `certify.py` in the same package already runs FVA with `processes=1` for a similar reason (its comment). Support: strong that the hang is in the cobra pool (stacks); weak for why the pool never finishes (the cause was not isolated). Nothing was changed in raven-toolbox.

## 2 RAVEN behaviour probes (GLPK, smallYeast, outputs only)

One JSON per probe in `probes/` and the script `probes/p2_*.m`. Signatures needing a correction: 2.4 (see there).

**2.1 `removeGenes`** (`2.1_removeGenes.json`): removing `YFR053C` gives 60 genes, 52 mets, 53 reactions with all four option variants (`removeUnusedMets`, `removeBlockedRxns` and `standardizeRules false` change nothing here); no gene is left outside the rules. A single-gene reaction (22 exist; first is `PGI`, gene `YBR196C`): after removing the gene `PGI` stays with `lb = ub = 0` (before -1000, 1000) and an empty `grRule`; genes 61 to 60, mets and reactions unchanged. A gene used by no reaction (`UNUSED1`, added with a zero column in `rxnGeneMat` and filled gene fields) is removed (`genes` 62 to 61). No gene field has a wrong length afterwards; the field names are the same before and after.

**2.2 `scoreModel`** (`2.2_scoreModel.json`), gene scores for G1 to G6 (`YAL054C`, `YBR117C`, `YBR196C`, `YBR218C`, `YCL040W`, `YCR012W`):

| Call | G1 | G2 | G3 | G4 | G5 | G6 |
| --- | --- | --- | --- | --- | --- | --- |
| tissue A | -5 | NaN | 0 | 3.4657 | 10 | NaN |
| tissue B | -5 | 8.0472 | 3.4657 | 3.4657 | -5 | NaN |
| G3 listed twice (levels 1 2 3 and 3 2 1), tissue A | -5 | NaN | 0 | 3.4657 | 10 | NaN |
| column A all NaN, tissue A | NaN | NaN | NaN | NaN | NaN | NaN |

The reaction score sums are -89.534 (A), -80.974 (B), -89.534 (G3 twice) and -106 (all NaN, 53 reactions at the default -2). The second G3 row is ignored (the first listing is used). HPA path: not run, no HPA fixture in `tests/function_tests/test_data`.

**2.3 `checkTasks`** (`2.3_checkTasks.json`, with the task file text): A (unknown `CHANGED RXN` `NOSUCHRXN`) `ok = true`. B (same metabolite in IN and OUT with nonzero LB) error `RAVEN:badInput`, "The IN LB and OUT LB in "[B] B" cannot be nonzero for the same metabolite". C (`a[s]; a[s]` twice) error `RAVEN:badInput`, "The constraints on some input(s) in "[C] C" are defined more than one time". D (`CHANGED LB 5`, `UB 1`) error with no identifier, "Invalid set of bounds for reaction(s): R2.". E (`CHANGED RXN` `R2; R3` with `CHANGED LB` `-5; 10`) error `MATLAB:innerdim` (matrix multiplication dimensions); the test file is in the JSON, so the cause can be checked. F (`EQU a[c] => b[c]` with `LB 7`, `UB 3`) error "Invalid set of bounds for reaction(s): TEMPORARY_1.".

**2.4 `FSEOF`** (`2.4_FSEOF.json`): the biomass reaction is `biomassOUT`. In smallYeast the nutrient inputs `glcIN` and `o2IN` have `ub = 0`, so the model does not grow as it is; the probe opens both to 1000 first. No reaction carries negative flux at the optimum, so as instructed the glucose exchange `glcIN` was taken and reversed. Result for the requested steps (`iterations` 10 and 20, `corrThreshold` 0.9 and 0.5, target `glcIN` reversed): no target, no nonzero slope (inference: a reversed target has no positive flux to enforce). After `reverseRxns` back, 44 slopes are nonzero but none passes the thresholds. The structure has `logical`, `amplify`, `knockdown`, `knockout` and `slope`, not `rxns`, so the probe reads ids from `model.rxns`. Added to the request, because `glcIN` gives no sign information: the product exchanges `ethOUT`, `acOUT` and `glyOUT` as they are (flux positive) and reversed (flux negative), `corrThreshold` 0.9 and 0.5. As they are: `ethOUT` at 0.9 gives targets `ethOUT, GPM, ENO, PYK, PDC, ADH1` with slopes 1, 0.069, 0.069, 0.117, 0.562, 1; `acOUT` and `glyOUT` give 17 and 11 targets (`extra` in the JSON). Reversed, all three give no target and no nonzero slope. So in RAVEN the sign of the target flux decides whether anything is found at all; the slopes of the targets are all positive.

**2.5 `applyCondition`** (`2.5_applyCondition.json`): `reset_exchanges` values and the reactions whose bounds changed (`glcIN`, `o2IN`, `ethIN`, all with `lb 0`, `ub` 0 to 1000):

| Value | Changed |
| --- | --- |
| `in`, `both`, `blocked`, `all` | `glcIN`, `o2IN`, `ethIN` (3) |
| `out`, `uptake`, `excrete`, `reverse` | none |
| `true`, `false` | error `MATLAB:invalidConversion`, "Conversion to char from logical is not possible." |
| `none`, `0`, `1` | error "Invalid reactionType specified" |

**2.6 `convertToIrrev`** (`2.6_convertToIrrev.json`): A (as it is) 71 reactions. B (`PGI`, `rev = 1`, set to `lb 0`, `ub 1000`): still split, `PGI` (0, 1000) and `PGI_REV` (0, 0); 71 reactions. C (`acOUT`, `rev = 0`, set to `lb -1000`): not split, `acOUT` (-1000, 1000); 71 reactions. So the split follows `rev`, not the bounds. The new reversed reactions are appended at the end in the order of the model (`..., PGI_REV, FBA_REV`). `convertToIrrev(model, {'PGI','FBA'})` and `{'FBA','PGI'}` give the same order of ids and the same `matchRev` (position 10 and 13 pointing at 54 and 55).

## 3 Solvers

`3.1_solvers.json`: `setRavenSolver('gurobi')` ok (Gurobi 13.0.1, tiny LP `OPTIMAL`, licence from the token server `gurobi.lic.chalmers.se`, port 27690). `scip`: "SCIP not found. RAVEN only provides the precompiled SCIP MEX binary for Windows." `cobra`: "COBRA toolbox is not initialized". `glpk` ok. Gurobi MATLAB is on the path only if `GUROBI_MATLAB` is added: `restoredefaultpath` drops `MATLABPATH`, so the test runner and the probe add it explicitly. The tiny-MILP part of this probe failed with an error in the probe's own call (`optimizeProb` argument form), not a solver result; the MILP tests of `tINIT` and `tGapfilling` with Gurobi are the evidence for 3.4.

3.2 Scenarios (the verdicts are from the full run; the scenario files name a solver only where shown, otherwise MATLAB uses the stored preference, which is GLPK):

| Scenario | Verdict | `matlab_solver` | `python_solver` |
| --- | --- | --- | --- |
| min_nr_fluxes_toy | MATCH | gurobi | gurobi |
| gapfill_connect_smallyeast | MATCH | gurobi | gurobi |
| find_leak_metabolite_toy | MATCH | not set (glpk) | not set |
| task_checking_smallyeast | MATCH | gurobi | gurobi |
| fit_parameters_toy | MATCH | not set (glpk) | not set |
| gap_report_toy | ERROR | not set (glpk) | not set |
| remove_bad_reactions_toy | MATCH | not set (glpk) | not set |
| init_merge_linear_smallyeast | MATCH | not set (glpk) | not set |
| init_scores_smallyeast | MATCH | not set (glpk) | not set |
| apply_condition_smallyeast | MATCH | gurobi | gurobi |
| close_model_smallyeast | MATCH | gurobi | gurobi |
| compare_fluxes_smallyeast | MATCH | not set (glpk) | not set |
| elemental_balance_smallyeast | MATCH | not set (glpk) | not set |
| model_manipulation_smallyeast | MATCH | not set (glpk) | not set |
| model_simplification_smallyeast | MATCH | not set (glpk) | not set |
| trace_flux_path_toy | MATCH | not set (glpk) | not set |
| walk_fluxes_smallyeast | MATCH | not set (glpk) | not set |
| gpr_dnf_rules | MATCH | not set (glpk) | not set |
| ec_fseof_ectestgem | DIFFER | gurobi | gurobi |
| ec_fva_ectestgem | MATCH | gurobi | gurobi |
| ec_model_expansion_ectestgem | MATCH | not set (glpk) | not set |
| ec_model_io_ectestgem | DIFFER | not set (glpk) | not set |
| enzyme_annotation_ectestgem | MATCH | not set (glpk) | not set |
| enzyme_usage_ectestgem | MATCH | gurobi | gurobi |
| flexibilize_concs_ectestgem | MATCH | gurobi | gurobi |
| fuzzy_kcat_coefficient_ectestgem | DIFFER | not set (glpk) | not set |
| kcat_chain_ectestgem | MATCH | not set (glpk) | not set |
| model_edits_ectestgem | DIFFER | not set (glpk) | not set |
| okp_kcat_ectestgem | MATCH | not set (glpk) | not set |
| protein_pool_ectestgem | MATCH | not set (glpk) | not set |
| sensitivity_tuning_ectestgem | MATCH | gurobi | gurobi |

The scenarios the handover lists as "MILP" mostly do not name Gurobi in `scenario.yml`; the 11 that do (among all 67) give `MATCH` except `ec_fseof_ectestgem` (documented) and `ec_model_full_yeastgem` (stale adapter). Python used Gurobi 13.0.3 (gurobipy) where `python_solver` says gurobi.

3.4 `tINIT` 31, `tGapfilling` 40, `tTasks` 7, all passed, none incomplete, with Gurobi. 3.5 SCIP is not available on Vera (error above), so no SCIP run was possible. Without Gurobi on the path the same three classes (`scip-tests/`) fall back to GLPK and pass with the same counts (31, 40, 7; status line `NO-GUROBI`), so none of those tests depends on Gurobi.

## 4 Hash seeds

For each scenario the Python side was run with `PYTHONHASHSEED` 0 to 19, each seed twice, and compared with the single MATLAB result of the full run (`hashseeds/<scenario>/summary.txt`). All 19 scenarios: 20 of 20 seeds `MATCH`, and the two runs of each seed are byte-identical (20 of 20). The 19 are the 18 handover scenarios in the full run plus `gap_report_toy` (with `processes=1`, see above): `merge_models_smallyeast`, `add_reactions_from_equations_smallyeast`, `add_reactions_from_model_smallyeast`, `change_gene_rules_smallyeast`, `standardize_gr_rules_smallyeast`, `gpr_dnf_rules`, `remove_genes_smallyeast`, `replace_metabolite_toy`, `remove_low_score_genes_smallyeast`, `remove_low_score_genes_ties_smallyeast`, `model_simplification_smallyeast`, `model_manipulation_smallyeast`, `homology_chain`, `get_gene_data_fixture`, `init_scores_smallyeast`, `init_merge_linear_smallyeast`, `gap_report_toy`, `yaml_roundtrip_smallyeast`, `task_list_edge_cases`.

## 5 RAVEN reading `name: ""`

Python (raven-toolbox `12f04bf`) wrote four variants of smallYeast to `probes/p5-files/`: empty metabolite name (`s_0001` does not exist in this model, so the first metabolite `MAL_m` was used), empty reaction name (first reaction), empty gene name (first gene) and a single space as metabolite name. Python warns only for the metabolite: "1 metabolite(s) have an empty name (MAL_m); the YAML is written anyway, but RAVEN's writeYAMLmodel refuses such a model.".

| File | `readYAMLmodel` | Result |
| --- | --- | --- |
| empty metabolite name | loads, no error | `metNames(1)` is empty; the first 5 `mets`, `metNames` and `metFormulas` are not shifted (only the name differs from the original); `checkModelStruct(model,false)` warns "The following metabolites have empty names: MAL_m" |
| empty reaction name | loads, no error | `rxnNames(1)` is empty; `checkModelStruct` gives no message |
| empty gene name | loads, no error | the file has no `name:` key for that gene (cobra drops it); `geneShortNames(1)` is empty; no message |
| single-space metabolite name | loads | `metNames(1)` is a space |

`writeYAMLmodel` with a single space as metabolite name writes `    - name: " "`.

## Changes to the handover signatures

2.4 as above (inputs opened, field names, extra product targets). 2.2 uses `arrayData.threshold = 1` as given. 2.3 task files were built from `taskTestModel` (copied into the probe) and the header of `test_tasks.txt`, one row per case. 3.1 added `addpath(GUROBI_MATLAB)`. The hash-seed run of `gap_report_toy` uses a wrapper that sets `cobra.Configuration().processes = 1`; nothing in RAVEN or raven-toolbox was edited, and RAVEN and raven-toolbox are back on `develop`.

Not checked: the HPA path of 2.2 (no fixture); SCIP and the cobra solver (not installed); the cause of the pool hang; whether `fuzzy_kcat_coefficient_ectestgem` and the three documented `DIFFER` scenarios still differ for the stated reasons (verdicts only); tier B probes (not built yet).

## Tier B: new scenarios (`tier-b/`)

Eight new scenarios on hand-written inputs, run on RAVEN `f667229` (branch `fix/export-for-git`, which only changes `exportForGit`) and raven-toolbox `c7af46a` (branch with the `exportForGit` changes and in-process FVA), seeds 0, 1 and 2. `DRAFT-inputs.md` lists the inputs as shown before the run. Each folder holds the compare output, the MATLAB result and the Python seed 0 result. RAVEN is not changed by the other scenarios, so their verdicts do not depend on the `exportForGit` branch.

| Scenario | Verdict | Differences (Python against MATLAB) |
| --- | --- | --- |
| `merge_linear_toy` | `MATCH` | none; in the three-reaction case both sides merge `R1` and `R2` (survivors `R0`, `R1`), also with the reactions in reverse order |
| `export_for_git_toy` | `MATCH` | on the branches: header, bounds with 1/3 and 0.30000000000000004, objective, files written; before the branches the header, the numbers and the default formats differed |
| `gap_fill_toy` | `DIFFER` | see below |
| `gap_fill_milp_toy` | `DIFFER` | see below |
| `homology_toy` | `DIFFER` | see below |
| `curate_tables_toy` | `DIFFER` | see below |
| `score_model_toy` | `DIFFER` | see below |
| `edge_model_toy` | `DIFFER` | RAVEN refuses to write two metabolites with one name in one compartment (Python writes them); `removeGenes` with defaults keeps `R3` (`G1 and G2`, rule emptied) where Python deletes it; the YAML header: Python writes `defaultLB`/`defaultUB` and an empty reaction name as `name: ""`, RAVEN writes neither. `convertToIrrev`/`convert_to_irreversible` agree on this model |

**`gap_fill_toy`** (`fillGaps` and `connect_blocked_reactions`, defaults): `fill_a` (`EX_A` upper bound 0.5): Python cannot connect `EX_D`, `EX_F`, `R2`, `R3`; RAVEN also lists `EX_A`, `EX_B`, `R1`. `fill_b` (`EX_A` 1.5, two blocked chains): Python raises "Gap-filling is infeasible"; RAVEN adds `T1`, `T2` and connects `EX_D`, `EX_F`, `R2`, `R3`. `fill_c` (net production allowed, `D` has no exit): Python cannot connect `R2`; RAVEN adds `T1` and connects it. `fill_d` with scores `{T1:-1, T2:-1, T3:0.5}` and `{T1:-2, T2:-1, T3:0}`: Python adds `T1, T3` and `T2, T3`; RAVEN adds `T1` and `T2`. `fill_d` without scores agrees. Python's default `eps` is 1.0; RAVEN uses 0.1 internally.

**`gap_fill_milp_toy`** (`gapFillMILP`/`fill_gaps_milp`, `gapFillFastLP`/`fill_gaps_fast_lp`): `milp_a` (`EX_A` and `BIO` bounded at 5000, minimum growth 1500, `bigM` 1000): Python adds `T1`; RAVEN does not solve it. `milp_b` (objective coefficient 2, minimum growth 15, supply 10): Python does not solve it; RAVEN adds `T1`. `milp_c` (objective `BIO + BIO2`, minimum growth 5): Python adds `T1`, `T2`; RAVEN adds `T2`. `milp_d` (template reaction `T2` with lower bound 1): Python adds `T1`, `T2`; RAVEN adds `T1`. `milp_e` (a template reaction with the id of a draft reaction): RAVEN adds it as `Rname_template`, Python adds nothing (MILP not solved; LP lists `R1`, `Rname` as unconnectable).

**`homology_toy`** (`getModelFromHomology`): `hom_a1_flag` (`(A or B) and C`, only `C` mapped): RAVEN keeps the unmapped subunit as `OLD_T_B and nC`, Python drops it and keeps `nC`. `hom_a2_drop`, `hom_a3_drop`: Python keeps the reaction with the mapped genes, RAVEN drops the reaction. `hom_b_t1t2` (preferred order `T1`, `T2`, `X` also hit by `T1`'s gene mapping): Python keeps `R_T1` and `R_T2`, RAVEN only `R_T1`. `hom_c_t2t1` (two templates with one reaction id, preferred order given): RAVEN fails with "Input A of class double and input B of class cell must be cell arrays of character vectors"; Python renames the second reaction (`R1_T1`). `makeFakeBlastStructure` needs at least 10 pairs, so the MATLAB side pads the list with pairs of genes that no template has.

**`curate_tables_toy`** (`curateModelFromTables`, two empty reactions plus the review's tables with `genes.tsv` extended to the genes in the rules): RAVEN links the two table reactions to the two existing reactions by `rxnIdx` (2 reactions); Python adds two new reactions (4). RAVEN reads the `metFormulas` header and the gene short names; Python does not. New ids: Python `M_1`, `R_0011`; RAVEN `M_0001`. Default bounds: RAVEN (-Inf, Inf), Python (-1000, 1000). With unpadded identifiers (`R_9`, `R_10`) RAVEN fails: `generateNewIds` sorts the identifiers as text, takes `9` as the largest and mints `R_10`, which exists.

**`score_model_toy`** (`scoreModel`, six genes, threshold 1 and 0): threshold 0 (reference level 0): RAVEN gives +10 where Python gives -5, and a zero level gives no score on the RAVEN side; reaction scores: Python returns NaN for reactions whose genes all lack a level where RAVEN returns the default -2 (3 reactions at threshold 1, 7 when all levels are missing).

**`export_for_git_toy`**: see above. The branches also give `Rxn id` as header and write an equation without a coefficient of 1.

Not checked: ftINIT and `prepINITModel` (the two sides take different inputs: RAVEN expression data, Python reaction scores; an interface decision comes first), `setGAM`, `scaleBiomassPseudoreaction`, `reporterMetabolites`, `getUniProtScores`, T2 as a fixture shared by several scenarios (`edge_model_toy` covers the cases for YAML, `convertToIrrev` and `removeGenes`), T4 and T6.
