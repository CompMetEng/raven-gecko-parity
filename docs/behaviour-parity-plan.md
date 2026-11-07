# Behavioural parity for RAVEN and raven-toolbox

The RAVEN ledger (`ledgers/raven.yml`) has 314 rows, of which 94 are marked `parity`: both sides implement the function and are intended to behave the same. Matching names prove nothing, as the README says, so a scenario runs both sides on the same input and compares the results. The repository holds 88 raven-pair scenarios (138 scenario folders in total, 50 of them for the gecko pair), and 85 of the 94 `parity` rows list at least one scenario. In the nightly report of 2026-11-04, 87 of the 88 scenarios match. The one that differs is `edge_model_toy`, in two details of the written YAML: whether a duplicated name is written, and the order of two genes. This document describes what exists, how it runs and what is still open. The gecko pair is described in [gecko-behaviour-parity-plan.md](gecko-behaviour-parity-plan.md).

## Scenarios by area

The full list is the folder `scenarios/`, and the rows that each scenario covers are in the `scenarios:` field of the ledger. The table names examples.

| Area | Examples |
|---|---|
| Editing and combining models | `add_reactions_from_equations_smallyeast`, `add_reactions_from_model_smallyeast`, `add_transport_reactions_smallyeast`, `change_gene_rules_smallyeast`, `change_reaction_equations_smallyeast`, `remove_genes_smallyeast`, `remove_metabolites_smallyeast`, `copy_to_compartment_smallyeast`, `merge_compartments_smallyeast`, `merge_models_smallyeast`, `contract_model_toy`, `sort_model_toy`, `scale_biomass_toy` |
| Formats | `yaml_roundtrip_smallyeast`, `yaml_float_digits_smallyeast`, `export_to_excel_smallyeast`, `delta_g_csv_smallyeast`, `save_delta_g_csv_smallyeast`, `export_for_git_toy`, `edge_model_toy` |
| Gene rules | `gpr_dnf_rules`, `standardize_gr_rules_smallyeast` |
| Reconstruction from homology | `homology_chain`, `homology_toy`, `hits_from_supporth_toy`, `hpa_scores_toy`, `uniprot_scores_toy`, `process_protein_fasta_fixture`, `remove_low_score_genes_smallyeast` |
| ftINIT, tasks and gap-filling | `init_scores_smallyeast`, `init_merge_linear_smallyeast`, `init_steps_toy`, `min_nr_fluxes_toy`, `ftinit_toy`, `init_output_model_smallyeast`, `ftinit_task_metabolites_toy`, `task_list_parsing`, `task_checking_smallyeast`, `gapfill_connect_smallyeast`, `gapfill_topology_smallyeast`, `gap_fill_toy`, `gap_fill_milp_toy` |
| Analysis and checks | `fseof_toy`, `reporter_metabolites_toy`, `random_sampling_toy`, `differential_sampling_toy`, `walk_fluxes_smallyeast`, `elemental_balance_smallyeast`, `check_model_struct_smallyeast`, `close_model_smallyeast`, `apply_condition_smallyeast`, `model_simplification_smallyeast` |

## How the comparison is built

**Both sides run the same binaries.** RAVEN and raven-toolbox both download BLAST+ 2.17.0, DIAMOND 2.1.17 and HMMER 3.4.0 from the release assets of `CompMetEng/toolbox-assets` (`downloadRavenBinaries.m`, `binaries.py`). On Windows both sides use HMMER 3.3.2, because no native Windows build of 3.4 exists. Identical binaries on identical inputs give identical hits, so the aligner steps are compared directly.

**Intermediate stages are compared as well as the result.** A whole-pipeline difference says that two models differ and not where. A chain emits named checkpoints in one run per implementation, and `parity compare` walks the results structurally, so the first stage that diverges is visible. `homology_chain` has two checkpoints (the BLAST hits and the draft model built from a synthetic ortholog mapping), and `yaml_roundtrip_smallyeast` has checkpoints as well. The other scenarios cover one function or one stage each. The ftINIT stages are separate scenarios (`init_scores_smallyeast`, `init_merge_linear_smallyeast`, `init_steps_toy`, `min_nr_fluxes_toy`, `ftinit_toy`).

**Solver scenarios name the solver.** Eleven raven-pair scenarios name Gurobi on both sides in their `scenario.yml`, so that a difference cannot be blamed on a solver preference. They do not run without a licence.

## Differences between the two tools

**Resolved.** Each of these was found by a scenario, tracked as an issue, and closed with the scenario matching. The ledger rows hold the current statement of each behaviour.

| Issue | Function | Subject |
|---|---|---|
| [#6](https://github.com/CompMetEng/raven-gecko-parity/issues/6) | `writeYAMLmodel` | layout of the written file, and the reader (see [yaml-reconciliation.md](yaml-reconciliation.md)) |
| [#7](https://github.com/CompMetEng/raven-gecko-parity/issues/7) | `checkTasks` | the linear program that checks a task |
| [#8](https://github.com/CompMetEng/raven-gecko-parity/issues/8) | `contractModel` | expand and contract as a round trip |
| [#9](https://github.com/CompMetEng/raven-gecko-parity/issues/9) | `addRxns`, `changeRxns` | names of metabolites created from an equation |
| [#10](https://github.com/CompMetEng/raven-gecko-parity/issues/10) | `addTransport` | the default call |
| [#12](https://github.com/CompMetEng/raven-gecko-parity/issues/12) | `changeGrRules` | appending a rule to a reaction without one |
| [#15](https://github.com/CompMetEng/raven-gecko-parity/issues/15) | `applyCondition` | the direction of `reset_exchanges` |
| [#16](https://github.com/CompMetEng/raven-gecko-parity/issues/16) | `deltaGCSV` | the "no measurement" sentinel of yeast-GEM |
| [#17](https://github.com/CompMetEng/raven-gecko-parity/issues/17) | `exportToExcelFormat` | hiding of default bounds |
| [#104](https://github.com/CompMetEng/raven-gecko-parity/issues/104) | `ftINIT` | pinning of tied optima with `resolveTies` |

**Remaining.**
- `predictLocalization` and `predict_localization` differ by design: RAVEN places reactions by a random search, and Python solves a MILP. The ledger row says so, and no scenario covers it.
- `contractModel`: the ledger row lists the remaining differences, which are the definition of a reversible reaction (`lb<0<ub` against `lb<0`) and the mirroring of bounds of opposite direction. `contract_model_toy` matches.
- `edge_model_toy` differs in the two YAML details named above.
- The ledger row of `FSEOF` says that `fseof_toy` gives DIFFER on three seeds, and the nightly report has `fseof_toy` as MATCH. The row and the scenario need to be reconciled.

## How it runs

### Nightly

`.github/workflows/nightly.yml` runs at 05:00 daily as a matrix over the pairs `raven` and `gecko`, one at a time.
1. **It skips when nothing moved.** `nightly/state.json` (for the gecko pair `nightly/gecko-state.json`) holds the revisions of the last comparison. The job resolves the tracked branches with `git ls-remote` and stops when the pair is unchanged. `workflow_dispatch` takes a `force` input and a `scenario` input.
2. **It uses the branches that `parity.toml` names.** `parity refs` reads them. All four repositories default to `main` on GitHub, and `parity.toml` tracks `develop` in all four. The raven entry checks out RAVEN and raven-toolbox. The gecko entry checks out all four, because geckopy installs raven-toolbox from git.
3. **It writes a report.** `nightly/report.md` (gecko: `nightly/gecko-report.md`) is regenerated and committed each time the comparison runs, so drift has a history.
4. **It separates a crash from agreement.** A scenario without a result file is recorded as `ERROR`, never as a match, and the job fails. One failing scenario does not stop the others.
5. **It provides a Gurobi licence.** The licence comes from the organisation secret `GUROBI_EDUARD` (a licence file) or from the three fields `GUROBI_WLSACCESSID`, `GUROBI_WLSSECRET` and `GUROBI_LICENSEID`. When neither is configured, the job says so and the solver-free scenarios still run. Only a WLS licence works on a hosted runner.
6. **It runs the Python side on several Python versions.** A `python-versions` job runs Python 3.11, 3.13, 3.14 and 3.15 against the MATLAB results of the compare job.

MATLAB comes from `matlab-actions/setup-matlab` (R2024b). A weekly job, `.github/workflows/gecko-genome-scale.yml` (Sundays 06:00), runs the genome-scale scenario `ec_model_full_yeastgem` of the gecko pair.

### On demand

```bash
parity run <scenario> --impl python
parity run <scenario> --impl matlab
parity compare results/python/<scenario>.json results/matlab/<scenario>.json --scenario <scenario>
```

`parity.local.toml` points a run at any checkout, and `parity run` takes `--output` and `--matlab`. `parity prepare`, `parity collect` and `parity seeds` (the Python side under several `PYTHONHASHSEED` values) support runs on a cluster, which [hpc-handover.md](hpc-handover.md) describes. `docs/workflow.md` describes the day-to-day workflow with the ledger.

## Open work

- **KEGG chain.** There is no scenario. `getKEGGModelForOrganism` and `getPhylDist` are `parity` rows without a scenario. `getPhylDist` needs a pinned KEGG artefact that the harness does not hold.
- **Compartment assignment.** `assignCompartments` has no scenario. The Python side certifies and gap-fills, and a scenario written before the two algorithms are diffed would measure the difference between two algorithms. The diff comes first, recorded on the ledger row.
- **Homology chain.** The checkpoint `orthologs` (the ortholog map after the cutoffs) is not in `homology_chain`.
- **Rows without a scenario.** Nine `parity` rows list none. Examples are `getDiamond`, `runSuppOrth`, `fitTasks` and `ftINITInternalAlg` (the ftINIT stages are covered through `ftinit_toy` and the other ftINIT scenarios, and the rows do not list them).
- **Pending rows.** 13 rows are `matlab-pending`.
- **Branch-targeted runs.** `--python-ref` and `--matlab-ref` and a single-command `parity verify` do not exist. A run against another branch uses `parity.local.toml`.

## Open questions

- **Where do the MATLAB results live between runs?** The nightly commits only the report and the state file, and uploads the results as a workflow artifact. Committing the results would let CI compare without MATLAB, and a stale committed result could become the reference. Recording the revisions in the result and warning on a mismatch avoids that.
- **Does the ledger need a `divergence` status?** The statuses are `parity`, `python-pending`, `matlab-pending`, `matlab-only`, `python-only`, `via-dependency`, `subsumed`, `internal` and `unreviewed`. A known difference is recorded in `notes` and in the `differs` field, which `parity ledger-audit` checks with the option `strict_differs`. No raven row uses `differs`.
- **How large should the fixtures be?** Small enough for a nightly, and large enough that a difference has somewhere to hide. The gecko pair has a weekly genome-scale job, and the raven pair has none. The genome-scale comparison behind the published Jaccard figures takes about one hour per side.
