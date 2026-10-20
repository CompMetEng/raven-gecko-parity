# Follow-up runs for the RAVEN vs raven-toolbox review

MATLAB R2024b on Vera, RAVEN `develop` at `79a7d8c`, GLPK. `probes/p1.m` and `probes/p4.m` are the scripts, `probes/p1.out` and `probes/p4.out` their output (`PROBE <n> <key>=<value>` lines), `hang/` the ten runs of `add_reactions_from_model_smallyeast`. Nothing under RAVEN or GECKO was changed.

## 1. `replace_metabolite_toy`

`replaceMets(model,'oxygen','o2')` renames `a` to `o2`, adds the row of `b` (name `o2`, the same name and compartment) to the row of `a`, and deletes `b`. The last two lines then run `repIdx = find(ismember(model.mets, repMetIds))` with `repMetIds = {'b'}`, which no longer exists, so `repIdx` is `zeros(0,1)`. `contractModel(model,'mets',[])` sees an empty `p.mets` and runs without the restriction (the extra unique column is never built), which is the same as `contractModel(model)` (step d): `unique` rows give `I=[1 3]`, `J=[1 1 2 2]`, so `unrelated1` and `unrelated2` merge. With the index of the surviving metabolite (`a`) the result keeps `unrelated2`.

## 2. `min_nr_fluxes_toy` on `e103c8c`

`MATCH` against all three Python seeds (run from a separate worktree).

## 3. The hang in `add_reactions_from_model_smallyeast`

Ten runs, one job each, one after another: 9 `MATCH` in about 40 s, run 8 hung (`hang/run8.txt`, `hang/run8.log`). MATLAB had printed the output of both `addRxnsGenesMets` calls in full, then sat at 0% CPU; the final `[matlab] -> ...` line never appeared. The stack from `kill -SEGV` shows only the main thread waiting in `mcr_run_main` (no MATLAB call stack). The node was not quiet (load 10 to 13).

## 4. Probes (`probes/p4.out`)

| Probe | Observed |
| --- | --- |
| P1 `mergeModels` | orphan metabolite and gene survive (53 mets, 62 genes); merged objective has 2 non-zero entries (model2 was an unmodified `smallYeast` and kept its own) |
| P2 `removeGenes` | removing `YBR196C` (only gene of `PGI`): mets 52 to 52, rxns 53 to 53, genes 61 to 60; `PGI` keeps its place with an empty `grRule`; nothing else pruned |
| P3 `addTransport` | `smallYeast` has compartments `c` and `m` (no `e`); `CO2` c to m gives `tr_0001`, `subSystems` = `Transport between c and m` |
| P4 `removeLowScoreGenes` | rules `G1`, `G3 or G4`, `G6 or G7`; removed `G2`, `G5` (the tie keeps `G1`; the order within `G3 or G4` and `G6 or G7` is changed) |
| P5 `parseTaskList` | `SHOULD FAIL` of `0` or `"0"`: `shouldFail=true` (only empty is false); `IN 'A[c];;B[c]'`: 3 inputs `{A[c], '', B[c]}`; `IN 'A[c]; B[c]'`: 2 inputs, the second `' B[c]'` with a leading space; `CHANGED RXN` with empty bounds: `LBrxn=[]`, `UBrxn=[]`, no error |
| P6 `writeYAMLmodel` | an empty metabolite name is refused (`RAVEN:modelError`, "The following metabolites have empty names"); a lower bound of `1/3` is written as `0.333333333333333`; `compartments:` follows `model.comps` order |
| P7 `exportForGit` | `'prefix'` and `'path'` are name-value arguments and `out` is not assigned; files: mat, txt, xlsx, xml, yml and `dependencies.txt` |
| P8 `isDnfGrRule` | `(A and B) or C` true; `A && B` true; `(A or B) and C` false; `a OR b` true; `A and` error "grRule ends after an operator: A and"; empty string true |
| P9 `checkTasks` | two metabolites named `A` (reaction on `a2`): `ok=false`, so only the first match is constrained |

Calls changed from the request: P3 uses `m` instead of `e`; P7 passes `prefix` and `path` by name.
