# Plan: behavioural parity for GECKO ↔ geckopy

The gecko ledger carries **78 rows marked `parity`**, 61 of them linked to at least one scenario. A scenario runs both sides on the same input and compares the results; a row without one is a claim about matching names. This plan describes the scenarios, the fixture they share and what each reports. It is the companion to [behaviour-parity-plan.md](behaviour-parity-plan.md), which does the same job for RAVEN.

## What exists now

48 scenarios, 25 that match and 23 that report differences. A scenario reports DIFFER while any case differs; the differences are described in its `scenario.yml` and in the ledger notes of the functions concerned, and a deliberate difference is said to be deliberate there.

| Scenario | Covers | Verdict |
|---|---|---|
| `ec_fseof_ectestgem` | `ecFSEOF` | match |
| `ec_fva_ectestgem` | `ecFVA` | match |
| `ec_model_build_edge_cases_ectestgem` | `makeEcModel` | match |
| `ec_model_expansion_ectestgem` | `makeEcModel`, `loadConventionalGEM`, `ModelAdapter` | match |
| `ec_model_full_yeastgem` | `makeEcModel`, `getECfromGEM`, `fuzzyKcatMatching`, `assignKcatValues`, `applyKcatConstraints`, `setProtPoolSize` | match |
| `enzyme_annotation_ectestgem` | `applyComplexData`, `getECfromGEM`, `getECfromDatabase` | match |
| `enzyme_usage_ectestgem` | `getEnzymeUsage`, `reportEnzymeUsage`, `getConcControlCoeffs` | match |
| `flexibilize_concs_ectestgem` | `loadFluxData`, `constrainFluxData`, `flexibilizeEnzConcs` | match |
| `fuzzy_kcat_coefficient_ectestgem` | `fuzzyKcatMatching`, `loadBRENDAdata` | match |
| `kcat_chain_ectestgem` | `fuzzyKcatMatching`, `loadBRENDAdata`, `writeDLKcatInput`, `readDLKcatOutput`, `mergeKcats`, `assignKcatValues`, `applyCustomKcats`, `getKcatAcrossIsozymes`, `getStandardKcat`, `removeStandardKcat`, `applyKcatConstraints` | match |
| `model_edits_ectestgem` | `getReactionsFromEnzyme`, `setKcatForReactions`, `applyECcodes`, `mapRxnsToConv`, `addNewRxnsToEC` | match |
| `molecular_weight_sequences` | `calculateMW` | match |
| `okp_kcat_ectestgem` | `findMetSmiles`, `writeOpenKineticsPredictorInput`, `readOpenKineticsPredictorOutput` | match |
| `protein_pool_ectestgem` | `setProtPoolSize`, `calculateFfactor`, `loadProtData`, `fillEnzConcs`, `constrainEnzConcs` | match |
| `sensitivity_tuning_ectestgem` | `sensitivityTuning`, `sigmaFitter`, `truncateValues` | match |
| `sigma_fitter_edge_cases_ectestgem` | `sigmaFitter` | match |
| `ec_model_io_crossload_ectestgem` | `loadEcModel`, `saveEcModel` | match |
| `ec_model_load_ecyeastgem` | `loadEcModel` | match |
| `ec_model_load_ecyeastgem911` | `loadEcModel` | match |
| `ec_model_load_light_humangem` | `loadEcModel` | match |
| `ec_model_build_light_humangem` | `makeEcModel`, `applyComplexData`, `getECfromGEM`, `fuzzyKcatMatching`, `readDLKcatOutput`, `mergeDLKcatAndFuzzyKcats`, `assignKcatValues`, `getStandardKcat`, `applyKcatConstraints`, `setProtPoolSize` | match |
| `ec_tuning_ecyeastgem` | `sensitivityTuning`, `sigmaFitter` | match |
| `ec_subset_genes_ectestgem` | `getSubsetEcModel` | match |
| `standard_kcat_subsystem_ectestgem` | `getStandardKcat` | match |
| `complex_data_edge_cases_ectestgem` | `applyComplexData`, `calculateFfactor`, `fillEnzConcs` | differs: partial stoichiometry `[2, 0]`, `calculateFfactor` with a missing sample, a zero total and an empty enzyme list, a negative abundance in `fillEnzConcs` |
| `ec_analysis_ecyeastgem` | `getEnzymeUsage`, `reportEnzymeUsage`, `getEnzymeBottlenecks`, `getPfbaEnzymes`, `getConcControlCoeffs`, `mapRxnsToConv`, `ecFSEOF`, `applyKcatConstraints`, `setKcatForReactions` | differs: the tied bottlenecks (none of the top 20 enzymes shared) and the order of 27 of the 372 ecFSEOF targets, which have the same target set |
| `ec_analysis_ecyeastgem911` | the same without `mapRxnsToConv` | differs: the tied bottlenecks and the order of 10 of the 363 ecFSEOF targets, which have the same target set |
| `ec_analysis_light_humangem` | `mapRxnsToConv`, `getSubsetEcModel`, `setProtPoolSize` | match |
| `ec_model_build_yarrowia` | `makeEcModel`, `getECfromGEM`, `getECfromDatabase`, `fuzzyKcatMatching`, `assignKcatValues`, `applyCustomKcats`, `getStandardKcat`, `applyKcatConstraints`, `setProtPoolSize` | differs: `getECfromGEM` for a model without EC codes |
| `ec_model_io_eciyali5` | `loadEcModel`, `saveEcModel`, `upgradeEcModel` | differs: the number of reactions `upgradeEcModel` turns (895 and 0) |
| `ec_proteomics_ecyeastgem` | `loadProtData`, `fillEnzConcs`, `constrainEnzConcs`, `calculateFfactor`, `setProtPoolSize`, `loadFluxData`, `constrainFluxData`, `flexibilizeEnzConcs`, `getConcControlCoeffs` | differs: the replicate-count check of `loadProtData`, the mean of `calculateFfactor`, the optimal vertex under the proteomics constraints |
| `ec_lookup_edge_cases_ectestgem` | `getECfromDatabase`, `findECInDB`, `applyECcodes`, `getSubsetEcModel` | differs: conflicting ECs for one gene |
| `ec_model_edit_edge_cases_ectestgem` | `addNewRxnsToEC`, `mapRxnsToConv`, `ecFVA` | differs: `addNewRxnsToEC` when no new enzyme is needed |
| `ec_model_io_ectestgem` | `saveEcModel`, `loadEcModel` | differs: content of the written files: `metaData` description, charge, SBO terms, `subsystem`, gene names |
| `enzyme_usage_edge_cases_ectestgem` | `getEnzymeUsage`, `reportEnzymeUsage`, `getEnzymeBottlenecks`, `getPfbaEnzymes`, `getConcControlCoeffs` | differs: the margin of a given pFBA fraction, tied bottlenecks, `zero=false` in the report, a missing usage reaction, a negative `topAbsUsage` |
| `evotune_edge_cases_ectestgem` | `screenKcatLeverage`, `selectTunableMask`, `cmaesKcatTuning`, `nChanged`, `parsimonyFrontier`, `bestParsimonious` | differs: order of tied rows, source label, a model without formulas, a flat objective, a NaN RMSE in `bestParsimonious` |
| `flexibilize_edge_cases_ectestgem` | `flexibilizeEnzConcs`, `relaxProteomicsGreedy`, `constrainEnzConcs` | differs: `constrainEnzConcs` with a prior lower bound or a negative concentration, an empty `enzymeSet`, the first of two tied enzymes |
| `flux_data_edge_cases_ectestgem` | `loadFluxData`, `constrainFluxData` | differs: argument handling of `loadFluxData` and `constrainFluxData` |
| `kcat_apply_edge_cases_ectestgem` | `applyKcatConstraints`, `setKcatForReactions`, `getKcatAcrossIsozymes`, `applyCustomKcats` | differs: the id vocabulary of the update list in `applyKcatConstraints` |
| `kcat_assignment_pipeline_ectestgem` | `kcatAssignmentConfig`, `collectKcats`, `assignKcats`, `ecTreeEstimate`, `kcatsFromECTree` | differs: config validation, doubled candidates and the `ec.source` text in `collectKcats`, NaN `k` and subsets in the EC tree |
| `kcat_io_edge_cases_ectestgem` | `findMetSmiles`, `writeOpenKineticsPredictorInput`, `saveKcatSource`, `loadKcatSource`, `getStandardKcat` | differs: `findMetSmiles` duplicates and unnamed metabolites, invalid rows in `loadKcatSource`, `getStandardKcat` |
| `kcat_merge_edge_cases_ectestgem` | `readDLKcatOutput`, `writeDLKcatInput`, `loadBRENDAdata`, `fuzzyKcatMatching`, `mergeKcats`, `mergeDLKcatAndFuzzyKcats`, `assignKcatValues` | differs: `OKP-` sources in `mergeKcats`, `readDLKcatOutput`, `writeDLKcatInput` |
| `kcat_review_stats_ectestgem` | `ecStatsFromBrenda`, `reviewKcatAssignment`, `loadBrendaMeasurements` | differs: `ecStatsFromBrenda`, the case of the source in `reviewKcatAssignment`, `measurements.tsv.xz` |
| `protein_data_filters_ectestgem` | `loadProtData` | differs: column-count validation of `loadProtData` |
| `sensitivity_tuning_edge_cases_ectestgem` | `sensitivityTuning`, `truncateValues` | differs: `Inf` in `truncateValues` |
| `upgrade_gecko3_ectestgem` | `upgradeEcModel`, `upgradeGECKOproject`, `convertCustomKcats` | differs: how the readers turn reverse protein reactions, the project report, a project with an adapter class, DLKcat output with foreign substrates |

The scenarios that match cover the base of the pipeline on GECKO's own unit-test model: the expansion into an ecModel, the annotation of its `ec` substructure, the protein budget, the kcat chain, OKP, enzyme usage, ecFVA, flexibilization, sensitivity tuning, ecFSEOF and calculateMW, and `ec_model_full_yeastgem` runs the full tutorial yeast-GEM through `makeEcModel` and a BRENDA kcat assignment. The `*_edge_cases_*` scenarios and `upgrade_gecko3_ectestgem`, `kcat_review_stats_ectestgem`, `kcat_assignment_pipeline_ectestgem`, `protein_data_filters_ectestgem` and `ec_model_io_crossload_ectestgem` feed the same functions inputs that the first fixtures skip, one input per case, so a difference names its cause. `ecTestGEM` is a seven-reaction model, so a scenario that needs a solver declares it explicitly (Gurobi on both sides).

## The fixture

`test/unit_tests/ecTestGEM` in the GECKO checkout: 7 reactions, 5 genes, and a hand-written
UniProt / KEGG / Complex Portal / BRENDA / proteomics snapshot. It is the model MATLAB's
`geckoCoreFunctionTests.m` pins its expected values against, and geckopy carries a port of that
suite (`tests/test_gecko_matlab_parity.py`) against its own copy under `examples/ecTestGEM`.

Both sides of a scenario read the **GECKO** copy. geckopy's adapter is loaded from its own
`model_adapter.toml` --- that file is one of the things being compared --- and only its `path`
and `conv_gem` are repointed at the GECKO fixture, so the parameter values stay geckopy's while
the data bytes are shared. Where the two would otherwise read different files, the scenario
declaration says so.

Three consequences worth stating up front:

**The SBML readers get tested for free.** MATLAB's own tests build the test model in code
(`getGeckoTestModel.m`) rather than reading `models/testModel.xml`. A scenario cannot do that
--- the two sides have to start from the same bytes --- so both read the SBML file, one through
RAVEN's `importModel` and one through cobrapy. That is why `ec_model_expansion_ectestgem` emits
a `conventional` checkpoint before it emits anything about `makeEcModel`: a divergence there is
the readers, not GECKO. Today they agree on every reaction, bound, metabolite, gene association
and stoichiometric coefficient.

**The adapters get compared, not assumed.** The same scenario emits an `adapter` checkpoint, so
a drift between GECKO's `TestGEMAdapter.m` and geckopy's `model_adapter.toml` is reported
rather than silently changing what the rest of the scenario measures.

**BRENDA is read from the same files on both sides.** GECKO's `databases/brenda-<release>` folder and geckopy's `brenda-refresh` output hold `kcat.tsv`, `mw.tsv` and `sa.tsv`, and both loaders read them. The scenarios that need BRENDA point both sides at one folder (`brenda_dir` for geckopy, the adapter's `params.path` for MATLAB). Confirmed equivalent, by execution on both sides: [docs/brenda-reconciliation.md](brenda-reconciliation.md).

## The scenario programme

Ordered by what each one needs, cheapest first. Every tier below is one scenario per row of the
table unless it says otherwise, and each emits checkpoints rather than a single verdict, so the
first stage that diverges is the one to fix.

### Tier 1 --- deterministic, offline, solver-free

| Scenario | Covers | Status |
|---|---|---|
| `ec_model_expansion_ectestgem` | `makeEcModel`, `loadConventionalGEM`, `ModelAdapter` | done, matching |
| `enzyme_annotation_ectestgem` | `applyComplexData`, `getECfromGEM`, `getECfromDatabase` | done, matching |
| `protein_pool_ectestgem` | `setProtPoolSize`, `calculateFfactor`, `loadProtData`, `fillEnzConcs`, `constrainEnzConcs` | done, matching |
| `molecular_weight_sequences` | `calculateMW` | done, matching |
| `ec_model_io_ectestgem` | `saveEcModel`, `loadEcModel` | done, round trips match, written files differ in content (see below) |
| `ec_model_io_crossload_ectestgem` | `loadEcModel`, `saveEcModel` | done, one `subsystem` difference in MATLAB's reader (see below) |
| `model_edits_ectestgem` | `addNewRxnsToEC`, `copyECtoGEM`, `setKcatForReactions`, `getReactionsFromEnzyme`, `mapRxnsToConv` | done, matching; the unknown-enzyme case is asserted on each side as a deliberate difference |

`molecular_weight_sequences` runs `calculateMW` and `calculate_mw` on eighteen sequences, each isolating one variable: the water-mass constant (`18.01528`), `X` ("any residue") and `B` ("D or N") as the mean and computed values of geckopy's table, lowercase letters (counted as residues after upper-casing), and a sequence with no recognised residue (`NaN` on both sides, so an empty sequence cannot read as an 18 Da protein). The twenty standard residues, `J`, `O`, `U` and `Z` agree between the tables to the gram-mole. All eighteen sequences agree, so the scenario is a regression fixture: a delta that reappears is a finding.

`ec_model_io_ectestgem` writes one freshly built full ecModel (17 reactions after expansion)
with both writers and reads each file back with its own reader. `saveEcModel` /
`save_ec_model` are thin wrappers around `writeYAMLmodel` / `write_yaml_model`, and
`loadEcModel` / `load_ec_model` around the matching readers. The model before any write and
the model after write-then-read agree on every field on both sides --- reactions, bounds, gene
associations, stoichiometry, and every `ec.*` array down to the coupling matrix --- so each
side's own round trip is lossless. The written *files* have the same layout (2-space indent,
the same key order and quoting, identical `ec-rxns`, `ec-enzymes` and `gecko_light` blocks)
and differ in content: 347 lines / 6774 characters from `save_ec_model` against 297 / 5496
from `saveEcModel`. None of the differences touches the `ec` data. They are:

- `metaData`: both writers record `gecko_format` and a `generator` line (GECKO in the MATLAB
  file, geckopy in the Python file); `save_ec_model` also writes a `description`, which
  `saveEcModel` cannot, because `writeYAMLmodel` has no slot for a model description.
- `charge: 0.0` on four metabolites in the Python file: the SBML fixture sets no charge and
  cobra's SBML import fills in 0; the MATLAB import leaves it unset and the writer omits it.
- SBO terms on metabolites, reactions and genes in the Python file: cobra keeps the `sboTerm`
  of every SBML entity (and `make_ec_model` adds `SBO:0000252` to the protein metabolites);
  the MATLAB import drops an SBO term that is the same for every entity of one type, as in
  this model.
- `subsystem: Protein usage` on `prot_pool_exchange` and the `usage_prot_*` reactions in the
  Python file: `make_ec_model` always sets it, `makeEcModel` only when the model has
  subsystems.
- `name` on the five genes in the MATLAB file: the MATLAB import reads the SBML gene label as
  the short name; cobra's import gives an empty name, which the writer omits.

A single `ec-code` is written as a list by both writers.

`ec_model_io_crossload_ectestgem` reads the two files with both readers. Charge, SBO terms, subsystems, gene names, EC codes, the `ec` data and the `metaData` description all load identically on both sides, from either file, including the file in which only some reactions carry a `subsystem` entry (as in the Python file). That needs the `readYAMLmodel` of RAVEN with CompMetEng/RAVEN#14, which keeps each subsystem with its reaction.

ecTestGEM ships no pre-existing `ecModel.yml`, so the `direct` checkpoint of this scenario is a freshly built ecModel. Reading ecModels that a GECKO or geckopy run wrote is the job of the real-model scenarios of tier 4: `ec_model_load_ecyeastgem`, `ec_model_load_ecyeastgem911`, `ec_model_load_light_humangem` and `ec_model_io_eciyali5`, which carry the legacy conventions of GECKO 3 files (reverse-direction `usage_prot_*` and `prot_pool_exchange` reactions, empty `enzymes` entries, notes with doubled quotes).

`model_edits_ectestgem` covers `addNewRxnsToEC` with more than one new reaction that needs isozyme (OR-rule) expansion in a single call, including one reversible reaction with an isozyme rule on its own, since the reversibility split duplicates that rule onto a second entry before the isozyme split runs. Both sides agree on the reactions added, the enzymes added and every field of the resulting model: reactions, gene associations, stoichiometry and every `ec.*` array.

Every checkpoint of the scenario matches. The one remaining behavioural difference is deliberate and decided in
raven-gecko-parity#70: for a protein id that is not in `ec.enzymes`, `getReactionsFromEnzyme` returns empty outputs and
geckopy's `get_reactions_from_enzyme` raises `ValueError`. Each side asserts its own behaviour inside the scenario (the run
fails if either changes) and reports the same empty result, so the checkpoint compares what both stand for.

The other checkpoints agree: `getReactionsFromEnzyme` matches protein ids case-insensitively on both sides;
`setKcatForReactions` gives the same kcat values and the same `ec.source` string for a suffixed id, a base name with one
scalar and a base name with a per-isozyme list; `copyECtoGEM` fills an empty eccode, leaves a populated one alone under
`overwrite=false` and treats an empty `ec.eccodes` entry under `overwrite=true` the same way; `mapRxnsToConv` matches on a
synthetic per-reaction flux vector keyed by id (the alphabetical rank of each ecModel reaction's own id, so the comparison
never assumes the two sides agree on row order, only on the set of ids), keeping the whole scenario solver-free.

`protein_pool_ectestgem` matches on `setProtPoolSize`, `calculateFfactor`, `loadProtData`, `fillEnzConcs` and `constrainEnzConcs`. MATLAB's `calculateFfactor` needs a `modelAdapter` in every call, to read `params.path`, where geckopy's `calculate_f_factor` needs one only when there is no `prot_data` at all: resolving that path does no I/O of its own. `loadProtData` is called with the filters switched on on both sides.

### Tier 2 --- the kcat chain

`kcat_chain_ectestgem` is one flowing ecModel with one checkpoint per stage, on ecTestGEM plus two reactions: `R2a` (the same gene rule shape as `R2` but no EC code, so it goes through DLKcat and the isozyme fill) and `R6` (no gene rule, the one shape the stock fixture lacks, needed for the new-standard-pseudo-protein path). It reports MATCH at every checkpoint.

| Checkpoint | Covers |
|---|---|
| `fuzzy` | `fuzzyKcatMatching`, `loadBRENDAdata` |
| `standard` | `getStandardKcat`, `removeStandardKcat` |
| `dlkcat` | `writeDLKcatInput`, `readDLKcatOutput`, `mergeDLKcatAndFuzzyKcats`, `mergeKcats` |
| `selection` | `assignKcatValues`, `applyCustomKcats`, `getKcatAcrossIsozymes` |
| `constraints` | `applyKcatConstraints` |

The behaviours the checkpoints pin down:

- `applyCustomKcats`, mode A (a reaction id and no protein): both sides write the new `ec.kcat`, set `ec.source` to `custom` and record the note.
- `assignKcatValues` with `criteria='median'` and `'mean'`: for one reaction with three candidates (kcat 100 from `custom` first, then 200 and 180 from `dlkcat`) both sides write 180 and 160 and record the source of the first row.
- `applyKcatConstraints`, light formulation, partial isozyme coverage: when only some isozymes of a reaction have a usable kcat, both sides skip the others before taking the cheapest. The shared ecTestGEM fixture cannot reach this case, so a small purpose-built light model inside the scenario does.
- `getStandardKcat` aggregates a subsystem's kcat as the median on both sides. The ecTestGEM fixture has no `subSystems` data, so the scenario exercises the global fallback only.

The BRENDA files are read by both loaders and are described in [docs/brenda-reconciliation.md](brenda-reconciliation.md). geckopy's BRENDA search order and its `origin` output match MATLAB's ([docs/brenda-search-order.md](brenda-search-order.md)). `runDLKcat` is out of scope because it needs the DLKcat container: the scenario reads a stored `DLKcat.tsv`.

`okp_kcat_ectestgem` covers OpenKineticsPredictor, the REST-submitted alternative to DLKcat: `findMetSmiles` feeding `writeOpenKineticsPredictorInput` (both `onlyWithSmiles` branches), then `readOpenKineticsPredictorOutput` parsing a result file into a kcat list. It stops at the CSV boundary on both sides. Both sides prefix each entry's source with `OKP-` (for example `OKP-CataPro`), and the scenario compares the source verbatim. It reports MATCH.

The edge cases of the chain are in `kcat_apply_edge_cases_ectestgem`, `kcat_merge_edge_cases_ectestgem`, `kcat_io_edge_cases_ectestgem`, `kcat_assignment_pipeline_ectestgem`, `kcat_review_stats_ectestgem` and `upgrade_gecko3_ectestgem`. Each reports DIFFER while a case differs, and its `scenario.yml` and the ledger notes of the functions concerned describe each difference.

### Tier 3 --- the solver tier

These scenarios need Gurobi, which `nightly.yml` provisions for the RAVEN pair. They compare at set level and on magnitudes, not on the exact vertex an LP happened to land on.

| Scenario | Covers | Verdict |
|---|---|---|
| `enzyme_usage_ectestgem` | `getEnzymeUsage`, `reportEnzymeUsage`, `getConcControlCoeffs` | match |
| `ec_fva_ectestgem` | `ecFVA` | match |
| `flexibilize_concs_ectestgem` | `loadFluxData`, `constrainFluxData`, `flexibilizeEnzConcs` | match |
| `sensitivity_tuning_ectestgem` | `sensitivityTuning`, `sigmaFitter`, `truncateValues` | match |
| `ec_fseof_ectestgem` | `ecFSEOF` | match |

`enzyme_usage_ectestgem`: with both of ecTestGEM's parallel m1 to m2 routes (R2, R4) open, MATLAB and geckopy solve to different, equally optimal vertices of the same LP once proteomics narrows every protein's cap and leaves the protein pool with slack. The scenario blocks R2 and R4 and leaves R3, ecTestGEM's only single-gene, single-isozyme reaction, as the only route, which removes the degeneracy structurally. P5's measured concentration, the tightest of the five, is then the only binding constraint on growth. `getConcControlCoeffs` reads the reduced cost of each usage bound in geckopy and perturbs each bound in MATLAB, and both give the same coefficients here, because nothing else binds within the perturbation range.

`ec_fva_ectestgem`: FVA needs no unique optimum, so R2 and R4 stay open, and R2's two isozymes get different kcats to give the solver a trade-off between them. Both sides read each canonical reaction's bounds off the solve of its own group, so the isozyme-split reaction agrees.

`flexibilize_concs_ectestgem`: `loadFluxData` and `constrainFluxData` match in both bound-setting modes (RAVEN's `setParam('var', ...)` sets `lb` and `ub` with the same sign-aware arithmetic as geckopy). `flexibilizeEnzConcs` reuses the R2/R4-blocked fixture of `enzyme_usage_ectestgem`, which has exactly one growth-limiting enzyme (P5): both sides flex it once from the same starting concentration and agree on the flexed concentration and the final growth. A harness note: MATLAB's `jsonencode` collapses a one-element numeric array to a bare scalar, so the MATLAB side converts such arrays with `num2cell` (see `docs/scenarios.md`).

`sensitivity_tuning_ectestgem`: `sensitivityTuning` bumps the same reaction and enzyme on both sides (`R3`, `P4`) with the same old and new kcat. `sigmaFitter` and `fit_sigma` return the same sigma and the same protein pool bound for the returned model. `truncateValues` agrees for one value per order-of-magnitude regime, a negative value, and a value small enough to round away to zero. The edge cases are in `sigma_fitter_edge_cases_ectestgem` (match) and `sensitivity_tuning_edge_cases_ectestgem` (differs).

**`ec_fseof_ectestgem`.** `geckopy.ec_fseof` is a thin wrapper over `raven_toolbox.analysis.fseof`, which classifies a reaction as a target by the correlation of its absolute flux with the enforced production flux, as RAVEN's `FSEOF` does, and `ecFSEOF` does the same. The scenario compares the enforced levels, the target set and the ranked order of the targets, on ecTestGEM:

1. The enforced production-flux levels. Both sides start at the biomass-optimal production flux
   (`0` here, since the cheaper of two isozyme routes covers all of it) and end at 90% of the
   theoretical maximum of the pFBA solution that maximises the target. Both fix the objective at
   0.999999 of its optimum in the pFBA step (`solveLP(model,1)` in RAVEN, `fseof` in raven_toolbox), so the levels are
   `[0, 259.19974, 518.39948, 777.59922, 1036.79896, 1295.998704]` on both sides. `ec_fseof` passes
   `equality=True`, which fixes the target flux at each level, as `ecFSEOF` does.
2. The candidate search space. Both sides consider reactions with a gene rule only, drop the
   `usage_prot_*` reactions from the report and keep the enforced production reaction itself
   (`R3`, whose flux follows the enforced level by definition). Uncatalysed and exchange
   reactions (`R1`, `R4`, `S1`, `S2`) are not targets on either side.
3. The selection criterion. Both sides take a reaction as a target when some absolute flux exceeds 1e-10, its spread is nonzero and the
   absolute Pearson correlation with the enforced flux reaches `corrThreshold` (default 0.9),
   and rank by slope times correlation. On this fixture both sides report `R3` and `R5`.

Only the levels, the target set and its ranked order are compared. The ranked order lists the reaction targets first and the transport targets after them, each group ranked by slope times correlation (the `transport` column of `ec_fseof` marks the second group; on ecYeastGEM 9.1.1 the 32 transport targets are the last 32 of 363 rows on both sides). `ecFSEOF.m` additionally computes a per-gene
essentiality column (one knockout LP per gene) and reports the absolute regression slope; `ec_fseof` does neither, pointing callers
at `cobra.flux_analysis.single_gene_deletion` and a subsystem/metabolite check instead.

### Tier 4 --- a real model

ecTestGEM proves the functions agree on every behaviour MATLAB itself pins down. It does not
prove they agree on yeast-GEM: 7 reactions leave a divergence nowhere to hide. The GECKO
tutorial ships `yeast-GEM.yml`, so a slow job --- separate from the nightly, as the RAVEN plan
says for its genome-scale comparison --- compares a full `makeEcModel` and a full BRENDA-based
kcat assignment against it.

`ec_model_full_yeastgem` runs with the other scenarios in the nightly workflow and weekly on its own (`.github/workflows/gecko-genome-scale.yml`); both download the BRENDA release into the GECKO checkout first (`scripts/ci/stage_brenda.py`). DLKcat is out of scope: it is a separate neural-network tool with its own PyTorch environment, and BRENDA fuzzy matching alone is a substantial genome-scale test of the shared pipeline. Two pieces of scenario-side data preparation are needed that no ecTestGEM-scale scenario has needed: `cobra.io.load_yaml_model` (what `load_conventional_gem` calls for a `.yml` file) does not populate `reaction.annotation['ec-code']` from yeast-GEM.yml's own `eccodes` reaction field, so the Python side reads the file a second time with a plain YAML parse to recover it; and `fuzzyKcatMatching`'s organism-closeness escalation needs a real KEGG taxonomy artefact this repo does not vendor or download, so a small two-organism `PhylDist.mat` stand-in (`YeastGEMAdapterParity.m`) exercises the same code path without it.

The scenario matches: both sides build the same 4671 ec reactions from yeast-GEM's 4102 conventional reactions, the same number of reactions with a kcat, the same protein pool bound and the same kcat for every reaction. Both loaders drop a BRENDA row whose kcat is more than four orders of magnitude below the median of its EC when five or more values exist (`filterOutlierKcats` and `filter_outlier_kcats`), so the 0.002 1/s kcat derived from a specific activity for EC 1.1.1.101 is dropped before matching on both sides. The organism-closeness tier of fuzzy matching is exercised only through the two-organism `PhylDist.mat` stand-in (`YeastGEMAdapterParity.m`), so both sides take the same wildcard-organism fallback for almost every query; the EC-code and substrate tiers are exercised at full genome scale. See `scenarios/ec_model_full_yeastgem/scenario.yml`.

The other real-data scenarios run the same pairs of functions on models that GECKO or geckopy wrote or built, so that a difference that needs many enzymes, a subsystem, a rule that real genes have or a file that an older version wrote has somewhere to show. Each needs a solver or a checkout of its data; the Vera paths are in the `inputs` of each `scenario.yml`.

| Scenario | Model and data | Verdict |
|---|---|---|
| `ec_model_load_ecyeastgem`, `ec_model_load_ecyeastgem911` | the tutorial ecYeastGEM written by GECKO, and ecYeastGEM 9.1.1 written by geckopy | match: every field of the loaded model, 8001 and 8010 reactions |
| `ec_model_load_light_humangem` | the light ecModel of the `light_ecModel` tutorial (16021 reactions, 30494 ec reactions), whose file has `enzymes: !!omap` entries without content | match: both readers load the file as it is |
| `ec_model_io_eciyali5` | eciYali5-GEM written by GECKO 3 (reverse protein reactions, three notes with doubled quotes, empty `enzymes` entries), read, upgraded, written and read back | differs in one value: `upgradeEcModel` turns 895 reactions in MATLAB and 0 in geckopy, whose reader has turned them already; the upgraded model agrees |
| `ec_model_build_light_humangem` | `makeEcModel` to `setProtPoolSize` on Human-GEM with the tutorial data (2421 enzymes, 29321 ec reactions, a 22 MB DLKcat file) | match |
| `ec_model_build_yarrowia` | the same pipeline on iYali5-GEM with the eciYali5-GEM data, once with the data as it is and once without double quotes | differs in one value per run: `getECfromGEM` raises for a model without EC codes and `fill_eccodes_from_gem` returns |
| `ec_analysis_ecyeastgem`, `ec_analysis_ecyeastgem911` | the usage, bottleneck, pFBA, control coefficient, mapping, ecFSEOF and kcat edit functions on ecYeastGEM | differs: tied bottlenecks and the order of the ecFSEOF targets, each from the choice of an optimal vertex |
| `ec_analysis_light_humangem` | growth, mapping and `getSubsetEcModel` to HT29-GEM on the light Human-GEM ecModel | match: `get_subset_ec_model` removes the genes and metabolites that no kept reaction uses, as `getSubsetEcModel` does |
| `ec_fva_ecyeastgem`, `ec_fva_light_humangem` | `ecFVA` and `ec_fva` on ecYeastGEM with flux data and on the light Human-GEM ecModel | `ec_fva_ecyeastgem` differs only in `fixed_cobra_tolerance`, the Python model at cobrapy's default tolerance (1944 of 4102 reactions differ), and matches with RAVEN's tolerance of 1e-9; `ec_fva_light_humangem` differs in the NaN records: RAVEN's Gurobi presolve reports the fixed-growth LP infeasible for 44 maxima and 164 minima of 300 sampled reactions (geckopy: 6 and 6), and of the 391 values both solve one differs (the minimum of MAR07157, 0 in geckopy and 2.45e-4 in MATLAB) |
| `ec_proteomics_ecyeastgem` | the proteomics and flux data functions on ecYeastGEM with the tutorial data and an edited file | differs: the replicate-count check of `loadProtData`, the average of `calculateFfactor` and the optimal vertex under the proteomics constraints |
| `ec_tuning_ecyeastgem` | `sensitivityTuning` and `sigmaFitter` on ecYeastGEM before its tuning (7, 12 and 5 tuned kcats) | match with Gurobi on both sides, with GLPK on both sides and with GLPK in MATLAB and scipy in geckopy, with the python models at RAVEN's tolerance of 1e-9 |
| `ec_subset_genes_ectestgem`, `standard_kcat_subsystem_ectestgem` | `getSubsetEcModel` with a context model that lacks a reaction and a gene of its rule, and `getStandardKcat` with subsystems that real models have, on ecTestGEM | match |

Fixes in the two toolboxes show in these scenarios, and each is the reason one of them matches. `sigmaFitter` takes the first sigma whose growth is within 1e-9 (relative) of the best growth, so the result no longer depends on the rounding noise of GLPK on a plateau (`sigma_fitter_edge_cases_ectestgem`). `getStandardKcat` and `assign_standard_kcat` treat subsystem names that differ by letter case as one subsystem and accept a reaction without a subsystem, `getSubsetEcModel` accepts a context model that is the full model and a light ecModel with the `001_` prefix, and geckopy reads the cells of `uniprot.tsv` like `textscan %q`, loads files with empty `enzymes` entries and doubled quotes, and works with the scipy solver (`geckopy.solver_support`, which has no MATLAB counterpart). The infeasible start of `flexibilizeEnzConcs` gives the same result on Gurobi, GLPK and scipy.

What stays different is a convention or follows the solver. RAVEN's `solveLP` gives Gurobi tolerances of 1e-9 and cobrapy leaves them at 1e-6, so the Python models of the FVA and analysis scenarios are set to RAVEN's value. An LP with many optimal vertices or a degenerate dual gives the sides different fluxes, shadow prices and control coefficients, which is why the usage functions are given one flux vector and why the tied bottlenecks, the target set of `ecFSEOF` and the enzymes analysed by `getConcControlCoeffs` differ. The remaining differences are in the ledger notes of the functions: the replicate-count check of `loadProtData`, the average of `calculateFfactor` (`nanmean` in geckopy), the KEGG fallback of `make_ec_model`, which is used when `kegg_db` is passed, and the reactions and metabolites that `getSubsetEcModel` removes with its unused genes.

Solver combinations. `PARITY_MATLAB_SOLVER` and `PARITY_PYTHON_SOLVER` replace the declared solvers of the scenarios that declare one. The default is Gurobi on both sides. Seed 0 was run with GLPK in MATLAB and scipy in geckopy and with Gurobi in MATLAB and scipy in geckopy on the solver-declaring scenarios that call `get_subset_ec_model`, `apply_kcat_constraints` or `switch_to_scipy` and on `ec_tuning_ecyeastgem` and `ec_proteomics_ecyeastgem`, and with GLPK on both sides on those two. `ec_fva_ecyeastgem` did not finish with GLPK in MATLAB (the MATLAB step passed a limit of one hour twice), so its GLPK combinations were not run, and `ec_fva_light_humangem` ran with Gurobi in MATLAB, whose step takes about 3 hours, and with Gurobi in MATLAB and scipy in geckopy. The scenarios that declare no solver (`ec_model_build_edge_cases_ectestgem`, `kcat_chain_ectestgem`, `ec_subset_genes_ectestgem`, `ec_model_load_light_humangem` and the other ecTestGEM scenarios without a solver input) ran once, with the default; the scipy path of `get_subset_ec_model` on ecTestGEM is covered by a unit test of geckopy.

- Matching on every combination that ran: `ec_fseof_ectestgem`, `ec_fva_ectestgem`, `enzyme_usage_ectestgem`, `ec_model_full_yeastgem`, `ec_model_build_light_humangem`, `ec_analysis_light_humangem`, `flexibilize_concs_ectestgem`, `sensitivity_tuning_ectestgem` and `sigma_fitter_edge_cases_ectestgem` (the sigma and the pool bound are the same on GLPK, Gurobi and scipy). `ec_tuning_ecyeastgem` matches with Gurobi on both sides, with GLPK on both sides and with GLPK in MATLAB and scipy in geckopy: the python models are set to RAVEN's solver tolerance of 1e-9 (`python_tolerance`), where cobrapy's GLPK tolerance of 1e-7 gave growth rates that differ in the sixth decimal, and with GLPK on both sides the growth rates before and after the tuning differ by 1.1e-10 to 6.6e-10 (absolute).
- The differing groups of the default run recur with every combination in the edge-case scenarios and in `ec_proteomics_ecyeastgem`, `ec_model_build_yarrowia` and `ec_model_edit_edge_cases_ectestgem`. With scipy geckopy lands on the vertex MATLAB reports for the tied rows of `getEnzymeBottlenecks` and `relaxProteomicsGreedy` in the ecTestGEM scenarios, so `tied_bottlenecks` and `relax_tied_enzymes` match there. The flux-data LP of `ec_proteomics_ecyeastgem` with 10 percent bounds is infeasible (the pyruvate exchange `r_2033` is bound to [7.657e-5, 8.463e-5] and pyruvate cannot be produced). Gurobi on both sides and Gurobi in MATLAB with scipy in geckopy report it infeasible. GLPK accepts it because its optimum violates the mass balance of `s_1400` by 7.657e-5, within its tolerance: MATLAB and geckopy with GLPK both report a growth rate of 0.0913870060941, and the case matches on GLPK on both sides, while HiGHS on scipy reports it infeasible, so it differs for GLPK in MATLAB with scipy in geckopy. This is a solver difference and not a defect of geckopy; the sides agree when both use the same solver.
- Differences on ecYeastGEM by solver pair, tutorial model and in brackets 9.1.1 (`ec_analysis_ecyeastgem`, `ec_analysis_ecyeastgem911`):

| MATLAB/geckopy | `bottlenecks` order and rows | `fseof` targets in MATLAB and geckopy | `fseof` ranks that differ | `min_uptake` | `pfba` |
|---|---|---|---|---|---|
| Gurobi/Gurobi | 20 and 49 (20 and 50) | 372 and 372 (363 and 363) | 27 (10) | none | none |
| GLPK/GLPK | 19 and 51 (16 and 33) | 358 and 357 (363 and 363) | 207 (188) | none (9.1.1: the pool flux is 125 in MATLAB and 113.7 in geckopy at the same minimal glucose uptake) | none |
| GLPK/scipy | 20 and 58 (20 and 42) | 358 and 380 (363 and 363) | 204 (188) | as GLPK/GLPK | none |
| Gurobi/scipy | 20 and 55 (20 and 56) | 372 and 380 (363 and 363) | 147 (11) | as GLPK/GLPK | none |

  The ranked FSEOF tables differ by ties (reactions of one pathway with equal scores to 1e-12) and near-ties, and the pFBA vertex of GLPK and HiGHS carries flux through other reactions of the same lipid pathway than the Gurobi vertex does. The enforced levels agree on every pair.
- scipy. `get_subset_ec_model` and `apply_kcat_constraints` edit the model on glpk when it is on the scipy solver and move it back to the sparse scipy problem in one pass (`geckopy.solver_support`), because optlang rewrites every constraint for each removed variable and builds the scipy matrix densely. The Python step of `ec_model_build_light_humangem` takes 6.8 minutes with GLPK in MATLAB and 6.9 minutes with Gurobi in MATLAB on scipy, against 2.1 minutes on Gurobi; the Python step of `ec_analysis_light_humangem`, which includes `get_subset_ec_model` on the light ecModel, takes 5.3 and 5.4 minutes on scipy against 5.8 minutes on Gurobi. Both scenarios match on scipy. Genome-scale models are still slow in the FVA: `ec_fva_ecyeastgem` on scipy (Gurobi in MATLAB) agrees with growth free, differs for one bound with growth fixed (the minimum of r_2720, -173.615 against -173.3495, 0.15 percent relative) and for 2257 of 8204 bounds in the run at cobrapy's tolerance (2258 of 24612 bounds against geckopy on Gurobi), because of the accuracy of HiGHS on the ill-scaled LPs of ecYeastGEM; that is a limit of the solver and not a parity statement, and the Python part takes about 100 minutes. `ec_fva_light_humangem` with scipy in geckopy (Gurobi in MATLAB) raises an OptimizationError (HiGHS status 'aborted') after 171 seconds, so it has no scipy result. `geckopy.solver_support` has no MATLAB counterpart.

A whole-model diff on top of an unchecked base reports many differences with one cause, so tiers 1 and 2 come first.

## The GECKO nightly

`nightly.yml` runs both pairs as two entries in one matrix, one at a time rather
than in parallel --- both entries commit `nightly/` back to this repo, and a concurrent push
from the other entry would race it, which running sequentially avoids without needing retry
logic. Each writes its own state and report (`nightly/state.json` / `nightly/report.md` for
raven, unchanged; `nightly/gecko-state.json` / `nightly/gecko-report.md` for gecko), so
neither entry's history overwrites the other's.

**A gecko run is a four-repo job.** geckopy depends on raven-toolbox (`git+...@develop`, not yet
on PyPI), and GECKO depends on RAVEN on the MATLAB side. Installing geckopy from its own
metadata would pull raven-toolbox's `develop` HEAD, which is *not* necessarily the revision
`parity.toml` tracks, and the MATLAB side needs a RAVEN checkout on the path regardless. So the
gecko entry checks out all four repos and installs raven-toolbox from the checkout before
geckopy (`pip install --no-deps` on geckopy itself, then its other dependencies read out of its
own `pyproject.toml` rather than duplicated by hand, so the workflow cannot drift out of sync
with them), or it would be comparing against something it did not record. `parity.toml` already
resolves the MATLAB half of this: the gecko pair's `setup` names `{raven_matlab}`, so RAVEN goes
on the MATLAB path from wherever that pair is configured, including a `parity.local.toml`
override. RAVEN and raven-toolbox are checked out unconditionally, in both matrix entries, so
the gecko entry always has them available even though raven is not the pair under test there.

`parity refs --format shell` prints a `<pair>_<side>_repo=` line next to its `_ref=` line, so the workflow can check out or `git ls-remote` a pair without a second, hand-maintained copy of its repo name: `parity.toml` stays the one place that fact is written down.

The scenarios that report documented differences are red: `complex_data_edge_cases_ectestgem`, `ec_analysis_ecyeastgem`, `ec_analysis_ecyeastgem911`, `ec_fva_ecyeastgem`, `ec_fva_light_humangem`, `ec_lookup_edge_cases_ectestgem`, `ec_model_build_yarrowia`, `ec_model_edit_edge_cases_ectestgem`, `ec_model_io_eciyali5`, `ec_model_io_ectestgem`, `ec_proteomics_ecyeastgem`, `enzyme_usage_edge_cases_ectestgem`, `evotune_edge_cases_ectestgem`, `flexibilize_edge_cases_ectestgem`, `flux_data_edge_cases_ectestgem`, `kcat_apply_edge_cases_ectestgem`, `kcat_assignment_pipeline_ectestgem`, `kcat_io_edge_cases_ectestgem`, `kcat_merge_edge_cases_ectestgem`, `kcat_review_stats_ectestgem`, `protein_data_filters_ectestgem`, `sensitivity_tuning_edge_cases_ectestgem`, `upgrade_gecko3_ectestgem`. Red is the correct, documented answer there, the same shape the raven nightly has with `yaml_roundtrip_smallyeast` and `task_checking_smallyeast`.

## Coverage by tier

1. Tier 1, the base of the pipeline: `ec_model_expansion_ectestgem`, `enzyme_annotation_ectestgem`, `protein_pool_ectestgem`, `molecular_weight_sequences`, `ec_model_io_ectestgem`, `ec_model_io_crossload_ectestgem` and `model_edits_ectestgem`.
2. Tier 2, the kcat chain: `kcat_chain_ectestgem`, `okp_kcat_ectestgem` and `fuzzy_kcat_coefficient_ectestgem`.
3. Tier 3, the solver tier: `enzyme_usage_ectestgem`, `ec_fva_ectestgem`, `flexibilize_concs_ectestgem`, `sensitivity_tuning_ectestgem` and `ec_fseof_ectestgem`.
4. Tier 4, real models: `ec_model_full_yeastgem`, `ec_model_build_light_humangem`, `ec_model_build_yarrowia`, `ec_model_load_ecyeastgem`, `ec_model_load_ecyeastgem911`, `ec_model_load_light_humangem`, `ec_model_io_eciyali5`, `ec_analysis_ecyeastgem`, `ec_analysis_ecyeastgem911`, `ec_analysis_light_humangem`, `ec_fva_ecyeastgem`, `ec_fva_light_humangem`, `ec_proteomics_ecyeastgem` and `ec_tuning_ecyeastgem`, with `ec_subset_genes_ectestgem` and `standard_kcat_subsystem_ectestgem` for the shapes that real models have.
5. Edge cases, which each report DIFFER while a case differs: the scenarios named `*_edge_cases_ectestgem`, `kcat_assignment_pipeline_ectestgem`, `kcat_review_stats_ectestgem`, `protein_data_filters_ectestgem` and `upgrade_gecko3_ectestgem`. `sigma_fitter_edge_cases_ectestgem` and `ec_model_build_edge_cases_ectestgem` match.
6. The nightly runs the gecko pair as a second matrix entry, with a four-repo checkout and its own state and report files.

The verdict of every scenario is in the table under "What exists now".

## Open questions

- **Which ecModel YAML content is the target?** The files differ in `metaData`, charge, SBO
  terms, subsystem and gene name (see `ec_model_io_ectestgem` above). The `ec-*` sections are
  identical; whatever is decided for the rest is decided for `writeYAMLmodel` and
  `write_yaml_model`, not for GECKO.
- **How much of MATLAB's multi-output API is a divergence worth asserting?** `applyComplexData`
  returns `foundComplex` / `proposedComplex` tables that geckopy only logs; `applyCustomKcats`
  and `findMetSmiles` are the same shape of difference. None of it changes a model, so no
  scenario can see it --- it belongs on the ledger rows, and today only some of it is there.
