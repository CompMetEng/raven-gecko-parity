# Tier B probes: draft inputs (not run)

Each case becomes a parity scenario (`scenario.yml`, `run.py`, `<id>.m`), so the same toy input runs on both sides and the verdict is `MATCH` or `DIFFER`. Models are written as lists of `(id, equation, lb, ub)`; both sides build them with the same equations. Metabolite names equal ids, compartment `c`.

## fillGaps / connect_blocked_reactions (review gapfill F1, F2, F4)

All cases: draft `EX_A: -> A [0, ub]`, `R1: A -> B`, `EX_B: B ->`, `R2: C -> D`, `EX_D: D ->`, `R3: E -> F`, `EX_F: F ->` (all `[0, 1000]` unless said); template `T1: A -> C`, `T2: B -> E` or `T2: A -> E` as stated.

| Case | Change | Recorded | Review expectation |
| --- | --- | --- | --- |
| fill_a | `EX_A` ub 0.5, template T1, T2: B -> E, threshold 1.0 and 0.1 | added, newly connected, cannot connect | RAVEN treats max flux <= 1 as blocked, Python about 0 |
| fill_b | `EX_A` ub 1.5, template T1, T2: A -> E, threshold 1.0 and 0.1 | same | forced flux per reaction 0.1 against 1.0 |
| fill_c | draft `EX_A`, `R1`, `EX_B`, `R2: C -> D` (D has no exit), template `T1: A -> C`, net production allowed | same | net-production path |
| fill_d | draft as fill_a with `EX_A` ub 1000 plus `EX_D`; template `T1: A -> C`, `T2: B -> C`, `T3: D -> B`; scores none, `{T1:-1, T2:-1, T3:0.5}`, `{T1:-2, T2:-1, T3:0}` | same | positive scores not clamped |

## gapFillMILP / fill_gaps_milp (review gapfill F6 to F8, plus a renamed id)

| Case | Draft | Template | Setting | Review expectation |
| --- | --- | --- | --- | --- |
| milp_a | `EX_A: -> A [0,5000]`, `BIO: B -> [0,5000]`, objective BIO | `T1: A -> B [0,5000]` | min growth 1500, bigM 1000 | bigM ignored for finite bounds on the Python side |
| milp_b | `EX_A [0,10]`, `BIO` objective coefficient 2 | `T1: A -> B` | min growth 15 | growth floor on c'v (RAVEN) or on v |
| milp_c | `EX_A [0,10]`, `BIO: B ->`, `BIO2: C ->`, objective BIO + BIO2 | `T1: A -> B`, `T2: A -> C` | min growth 5 | same |
| milp_d | `EX_A [0,10]`, `BIO`, `EX_C: C -> [0,1000]` | `T1: A -> B`, `T2: A -> C [1,1000]` | min growth 5 | template reaction with lb > 0 is always added |
| milp_e | `R1: a -> b`, `Rname: b ->`, objective R1 | `Rname: -> a` (same id as a draft reaction), also run on fastLP and connect | none | renamed-id handling |

## mergeLinear (review INIT F16, F17)

| Case | Input | Recorded | Review expectation |
| --- | --- | --- | --- |
| ml_a | `R0: <=> M0`, `R1: M1 + M2 <=> M0`, `R2: M1 <=> M0 + M2`, all `[-1000,1000]`, metabolite order M0, M1, M2 | surviving reaction ids and merged groups | a Python re-implementation of `mergeLinear.m` merges R1 and R2 (survivors R0, R1); `merge_linear` merges all three (survivor R0) |
| ml_b | `smallYeast` after the task-free reduction | counts and survivors | no difference expected |

The expected MATLAB column rests on the reviewer's Python simulation of `mergeLinear.m`, not on a MATLAB run.

## prepINITModel / ftINIT output model (review INIT F1, F2)

`smallYeast.yml` with `glcIN` and `o2IN` upper bound 1000, the tasks of `task_checking_smallyeast`, with the infeasible task `GLC_AER_TIGHT` renamed `GLC_AER` (the id of a feasible one). Scores: reactions with a gene rule, alternately +1 and -1 in model order. Recorded: ids and descriptions of the kept tasks, merged essential reactions, number of reference-model reactions whose bounds or stoichiometry differ from the template, and bounds and equation of `MDH1`, `PGK`, `ENO` in the `ftINIT` output (no gap filling). Review expectation: the Python output model is built from the oriented, tightened model (16 bounds, 1 stoichiometry differ), and the same-id task is dropped.

## getModelFromHomology (review recon F2, F3, F4, F6)

Template reactions `a_c -> b_c` with gene rules as in the review scripts; ortholog hits as listed.

| Case | Rules and hits | Settings | Recorded |
| --- | --- | --- | --- |
| hom_a | `(A or B) and C` with only C mapped; `A or (B and C)` with B unmapped; `(A and B) or C` with A and B unmapped | complex policy flag, keep, drop | resulting rule per reaction |
| hom_b | T1: R_T1 gene Y; T2: R_T2 gene X; hits Y -> N1 (T1), X -> N1 and N2 (T2) | preferred order none, T1 T2, T2 T1 | reactions, rules, gene map |
| hom_c | two templates, same reaction id R1 (rules Y and X) | preferred order none and T2 T1 | reaction rule and model name (merge order) |

The MATLAB side needs a `blast` structure built from the same hits; the exact fields will follow `getBlastFromExcel` and `getModelFromHomology` before this is run.

## curateModelFromTables (review recon F15 to F17)

Toy model with two empty reactions `R_9`, `R_10`; tables `mets.tsv`, `rxns.tsv`, `coeffs.tsv` as in the review folder `cur/`, with header aliases in the header row, reactions linked by `rxnIdx` where the name differs, and new reactions without an id. Recorded: metabolite names, formulas, charges, reaction equations, bounds, rules, and the new ids.

## Not drafted yet

`setGAM` and `scaleBiomassPseudoreaction` (review io F22), `reporterMetabolites` ties, `getUniProtScores` column order, `exportForGit` (pending your decision).
