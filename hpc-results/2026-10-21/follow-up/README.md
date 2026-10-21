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

## 5. The hang, second round (`hang-hunt/`)

With a trace in `parity_run.m`, `finish.m` and the scenario file (marks with timestamps), the scenario ran 40 times in a row (`sequential-40-summary.txt`) and 8 workers x 15 times at once on one node, each with its own copy of the checkout (`parallel-8x15-summary.txt`): 160 of 160 `MATCH`, no hang. A complete trace takes about 0.4 s from `scenario start` to `finish.m reached` (`example-trace.txt`).

An earlier version of the parallel job had a setup error: the checkout copies had no sibling `RAVEN` folder, so MATLAB stopped at `readYAMLmodel` with an error. In 3 of 79 of those runs the MATLAB process was still alive 4 minutes after `finish.m` had run (`HANG-w*.txt`, `hang-trace-w*.txt`; node load about 20). All threads were waiting on locks or events at 0% CPU. So a MATLAB process can stall while exiting on a busy node. The hang of run 8 matches that picture (output finished, process idle), but it happened in a run that succeeded, so the match is partial. Support: moderate. A time limit per MATLAB attempt with one retry turns such a stall into a clean result.

## 6. `replaceMets` fix (`replaceMets-fix/`)

Patch against RAVEN `79a7d8c` (two files), run on MATLAB R2024b with GLPK:

| Step | Result |
| --- | --- |
| A: new tests only, `replaceMets.m` unpatched | `tManipulation` 84 passed, 1 failed, 0 incomplete; `replaceMetsByNameLeavesUnrelatedDuplicatesAlone` fails, the other two new tests pass |
| B: full patch | 85 passed, 0 failed, 0 incomplete; all three new tests pass |
| C: `replace_metabolite_toy`, MATLAB against Python seed 0 only | `MATCH` |

Applying needed `git apply --recount --ignore-whitespace`; `=&gt;` in two comments of `duplicateModel` was written as `=>`. Nothing else was changed. The RAVEN clone was restored to `79a7d8c` afterwards. Per-test rows: `tManipulation-A-tests-only.csv`, `tManipulation-B-full-patch.csv`.
