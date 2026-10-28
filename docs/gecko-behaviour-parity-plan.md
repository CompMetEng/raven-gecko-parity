# Plan: behavioural parity for GECKO ↔ geckopy

The gecko ledger carries **78 rows marked `parity`**, 61 of them linked to at least one scenario. A scenario runs both sides on the same input and compares the results; a row without one is a claim about matching names. This plan describes the scenarios, the fixture they share and what each reports. It is the companion to [behaviour-parity-plan.md](behaviour-parity-plan.md), which does the same job for RAVEN.

## What exists now

33 scenarios, 16 that match and 17 that report differences. A scenario reports DIFFER while any case differs; the differences are described in its `scenario.yml` and in the ledger notes of the functions concerned, and a deliberate difference is said to be deliberate there.

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
| `complex_data_edge_cases_ectestgem` | `applyComplexData`, `calculateFfactor`, `fillEnzConcs` | differs: partial stoichiometry `[2, 0]`, `calculateFfactor` with a missing sample, a zero total and an empty enzyme list, a negative abundance in `fillEnzConcs` |
| `ec_lookup_edge_cases_ectestgem` | `getECfromDatabase`, `findECInDB`, `applyECcodes`, `getSubsetEcModel` | differs: conflicting ECs for one gene |
| `ec_model_edit_edge_cases_ectestgem` | `addNewRxnsToEC`, `mapRxnsToConv`, `ecFVA` | differs: `addNewRxnsToEC` when no new enzyme is needed |
| `ec_model_io_crossload_ectestgem` | `loadEcModel`, `saveEcModel` | differs: `subsystem` in MATLAB's YAML reader |
| `ec_model_io_ectestgem` | `saveEcModel`, `loadEcModel` | differs: content of the written files: `metaData` description, charge, SBO terms, `subsystem`, gene names |
| `enzyme_usage_edge_cases_ectestgem` | `getEnzymeUsage`, `reportEnzymeUsage`, `getEnzymeBottlenecks`, `getPfbaEnzymes`, `getConcControlCoeffs` | differs: the pFBA margin, tied bottlenecks, `zero=false` in the report, a missing usage reaction, a negative `topAbsUsage` |
| `evotune_edge_cases_ectestgem` | `screenKcatLeverage`, `selectTunableMask`, `cmaesKcatTuning`, `nChanged`, `parsimonyFrontier`, `bestParsimonious` | differs: order of tied rows, source label, a model without formulas, a flat objective, a NaN RMSE in `bestParsimonious` |
| `flexibilize_edge_cases_ectestgem` | `flexibilizeEnzConcs`, `relaxProteomicsGreedy`, `constrainEnzConcs` | differs: an infeasible start, `constrainEnzConcs` with a prior lower bound or a negative concentration, an empty `enzymeSet`, tied enzymes |
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

`ec_model_io_crossload_ectestgem` reads the two files with both readers. Charge, SBO terms,
gene names, EC codes, the `ec` data and the `metaData` description all load identically on
both sides, from either file. The one difference is `subsystem` in the MATLAB reader: when only
some reactions carry a `subsystem` entry, as in the Python file, `readYAMLmodel` flattens the
cell array with `[model.subSystems{:}]`, which drops the empty cells, so the six
`Protein usage` entries land on the first six reactions of the file and not on the reactions
they were written for.

One difference from `yaml_roundtrip_smallyeast` worth stating for whoever extends this
scenario: ecTestGEM ships no pre-existing `ecModel.yml`, so this scenario's `direct`
checkpoint is a freshly *built* ecModel, not one read off disk the way `yaml_roundtrip`'s is
from `smallYeast.yml`. Reading a real, larger ecModel --- the tutorial's `ecYeastGEM.yml`,
written by MATLAB --- on both sides remains open, and is a materially different question: it
would test today's readers against a file that may carry the legacy quirks
`raven_toolbox.io.read_yaml_model`'s own docstring says it normalises on load (a top-level
`smiles` field, reverse-direction `usage_prot_*` / `prot_pool_exchange`), none of which a
freshly-built model exercises. Worth its own scenario later, not a substitute for this one.

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
   theoretical maximum: `[0, 259.2, 518.4, 777.6, 1036.8, 1296.0]`. The levels agree to within
   1e-6 relative. MATLAB's `solveLP(model,1)` relaxes the objective it fixes to 0.999999 of its
   maximum, so the top level is 1295.998704, and the scenario tolerance is 2e-6 for that reason.
2. The candidate search space. Both sides consider reactions with a gene rule only, drop the
   `usage_prot_*` reactions from the report and keep the enforced production reaction itself
   (`R3`, whose flux follows the enforced level by definition). Uncatalysed and exchange
   reactions (`R1`, `R4`, `S1`, `S2`) are not targets on either side.
3. The selection criterion. Both sides take a reaction as a target when some absolute flux exceeds 1e-10, its spread is nonzero and the
   absolute Pearson correlation with the enforced flux reaches `corrThreshold` (default 0.9),
   and rank by slope times correlation. On this fixture both sides report `R3` and `R5`.

Only the levels, the target set and its ranked order are compared. `ecFSEOF.m` additionally computes a per-gene
essentiality column (one knockout LP per gene), reports the absolute regression slope and splits
results into `rxnTargets` vs `transportTargets`; `ec_fseof` does none of these, pointing callers
at `cobra.flux_analysis.single_gene_deletion` and a subsystem/metabolite check instead.

### Tier 4 --- a real model

ecTestGEM proves the functions agree on every behaviour MATLAB itself pins down. It does not
prove they agree on yeast-GEM: 7 reactions leave a divergence nowhere to hide. The GECKO
tutorial ships `yeast-GEM.yml`, so a slow job --- separate from the nightly, as the RAVEN plan
says for its genome-scale comparison --- compares a full `makeEcModel` and a full BRENDA-based
kcat assignment against it.

`ec_model_full_yeastgem` runs with the other scenarios in the nightly workflow and weekly on its own (`.github/workflows/gecko-genome-scale.yml`); both download the BRENDA release into the GECKO checkout first (`scripts/ci/stage_brenda.py`). DLKcat is out of scope: it is a separate neural-network tool with its own PyTorch environment, and BRENDA fuzzy matching alone is a substantial genome-scale test of the shared pipeline. Two pieces of scenario-side data preparation are needed that no ecTestGEM-scale scenario has needed: `cobra.io.load_yaml_model` (what `load_conventional_gem` calls for a `.yml` file) does not populate `reaction.annotation['ec-code']` from yeast-GEM.yml's own `eccodes` reaction field, so the Python side reads the file a second time with a plain YAML parse to recover it; and `fuzzyKcatMatching`'s organism-closeness escalation needs a real KEGG taxonomy artefact this repo does not vendor or download, so a small two-organism `PhylDist.mat` stand-in (`YeastGEMAdapterParity.m`) exercises the same code path without it.

The scenario matches: both sides build the same 4671 ec reactions from yeast-GEM's 4102 conventional reactions, the same number of reactions with a kcat, the same protein pool bound and the same kcat for every reaction. Both loaders drop a BRENDA row whose kcat is more than four orders of magnitude below the median of its EC when five or more values exist (`filterOutlierKcats` and `filter_outlier_kcats`), so the 0.002 1/s kcat derived from a specific activity for EC 1.1.1.101 is dropped before matching on both sides. The organism-closeness tier of fuzzy matching is exercised only through the two-organism `PhylDist.mat` stand-in (`YeastGEMAdapterParity.m`), so both sides take the same wildcard-organism fallback for almost every query; the EC-code and substrate tiers are exercised at full genome scale. See `scenarios/ec_model_full_yeastgem/scenario.yml`.

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

The scenarios that report documented differences are red: `complex_data_edge_cases_ectestgem`, `ec_lookup_edge_cases_ectestgem`, `ec_model_edit_edge_cases_ectestgem`, `ec_model_io_crossload_ectestgem`, `ec_model_io_ectestgem`, `enzyme_usage_edge_cases_ectestgem`, `evotune_edge_cases_ectestgem`, `flexibilize_edge_cases_ectestgem`, `flux_data_edge_cases_ectestgem`, `kcat_apply_edge_cases_ectestgem`, `kcat_assignment_pipeline_ectestgem`, `kcat_io_edge_cases_ectestgem`, `kcat_merge_edge_cases_ectestgem`, `kcat_review_stats_ectestgem`, `protein_data_filters_ectestgem`, `sensitivity_tuning_edge_cases_ectestgem`, `upgrade_gecko3_ectestgem`. Red is the correct, documented answer there, the same shape the raven nightly has with `yaml_roundtrip_smallyeast` and `task_checking_smallyeast`.

## Coverage by tier

1. Tier 1, the base of the pipeline: `ec_model_expansion_ectestgem`, `enzyme_annotation_ectestgem`, `protein_pool_ectestgem`, `molecular_weight_sequences`, `ec_model_io_ectestgem`, `ec_model_io_crossload_ectestgem` and `model_edits_ectestgem`.
2. Tier 2, the kcat chain: `kcat_chain_ectestgem`, `okp_kcat_ectestgem` and `fuzzy_kcat_coefficient_ectestgem`.
3. Tier 3, the solver tier: `enzyme_usage_ectestgem`, `ec_fva_ectestgem`, `flexibilize_concs_ectestgem`, `sensitivity_tuning_ectestgem` and `ec_fseof_ectestgem`.
4. Tier 4, a real model: `ec_model_full_yeastgem`.
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
