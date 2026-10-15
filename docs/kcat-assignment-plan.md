# Plan: kcat sources, assignment, reference data and change tracking

Status: plan, 2026-10-15. Replaces the 2026-10-05 version of this file.
Background:
- `kcat-assignment-options.md` (designs A to C, and the choices C1 to C9
  that led to this plan);
- `model-change-tracking-options.md` (option 2, decisions E1 to E4);
- the kcat-benchmark results (`docs/results.md` in that repository),
  which set the presets and defaults.

Everything here is implemented identically in geckopy and GECKO
(MATLAB), with parity scenarios in this repository. The project
configuration is `model_adapter.toml`, shared by both tools
(GECKO#486, geckopy#61, issue #132).

## 1. Decisions

| Topic | Decision |
|---|---|
| Stored kcats (C1) | One file per kcat source in the project's `data/kcats/` |
| Key of a stored row (C2) | Conventional reaction (`_REV` for direction) + enzyme set + substrate; never `model.ec.rxns` |
| Granularity (C3) | One value per reaction per source; aggregation (max, median) happens in the gatherer |
| Bookkeeping (C4) | Self-describing files: `#` header lines with source, generator, parameters and data version. Files are tracked in git, so hand edits show in diffs. No registry file |
| Order of sources (C5) | `[kcat_assignment]` in `model_adapter.toml` |
| Choosing a value (C6) | Stages: `select` (first source wins, with the BRENDA quality tiers), `combine` (weighted by uncertainty), `fill` (isozymes, standard kcat). Named presets for common orders |
| Gatherers (C7) | One dispatcher, `gatherKcats(model, source, ...)`, writes the source file; `saveKcatSource` writes any kcat list (hand-made, special cases) |
| Custom kcats (C8) | An ordinary curated source, `data/kcats/custom.tsv`, placed in the order |
| Entry point (C9) | `assignKcats(model)`; the building blocks stay callable |
| Default preset | `gecko_median` in the project template now; `ec_tree_homologs` once the BRENDA measurement files are released |
| Reference data | Release assets in `CompMetEng/toolbox-assets`, which raven-data becomes (transferred and renamed, M1). All BRENDA files form one asset release; the DLKcat files form another |
| One-way rule | Source files feed `model.ec.kcat`; nothing in the model is written back to them |
| Tracking | Change log next to the model and the `ec.source` format (sections 9 and 10, unchanged from the previous plan) |

## 2. Two layers of data

| Layer | What | Where | Who writes it |
|---|---|---|---|
| Reference data | BRENDA kcats, specific activities, molecular weights, per-measurement table, sequences, substrate dictionary, combination parameters; DLKcat model files | `CompMetEng/toolbox-assets` release assets, cached locally | Maintainers, with geckopy's `brenda-refresh` |
| Project kcat files | The kcats gathered for one model, one file per source, plus curated `custom.tsv` | `data/kcats/` in the project, in git | `gatherKcats`, `saveKcatSource`, and people (custom) |

## 3. Reference data in toolbox-assets

### 3.1 The repository

- `SysBioChalmers/raven-data` was transferred to CompMetEng and renamed
  `toolbox-assets` on 2026-10-16. Releases, tags and assets moved
  unchanged, the repository is public, and GitHub redirects the old URLs,
  so released RAVEN and raven-toolbox versions keep working.
- RAVEN and raven-toolbox point to the new repository: RAVEN#764,
  raven-toolbox#172.
- `manifest.json` gets a new version (`manifest-v2`) with the new URLs and
  the GECKO entries; every file keeps its SHA-256 and size.
- The README covers the GECKO assets and their licences.

### 3.2 Downloading

| Tool | Function | Change |
|---|---|---|
| RAVEN | `fetchToolboxAsset(destDir, tag, asset)` | Replaces `fetchRavenDataAsset`, which is removed (all its callers in RAVEN are updated); base URL `https://github.com/CompMetEng/toolbox-assets/releases/download/`; checks the SHA-256 from the manifest, as Python already does |
| raven-toolbox | `raven_toolbox.manifest`, `data`, `binaries` | Baked registries, manifest and publishing scripts point to toolbox-assets; `RAVEN_PYTHON_MANIFEST` still overrides the manifest |
| GECKO | BRENDA and DLKcat loaders | Fetch through `fetchToolboxAsset` into `databases/` on first use |
| geckopy | BRENDA loaders, `run_dlkcat` | Fetch through raven-toolbox's registry |

GECKO and geckopy no longer carry BRENDA files in their repositories. The
first BRENDA use downloads a few MB, as RAVEN already does for KEGG data.

### 3.3 The BRENDA asset

One release per BRENDA release and build, tag `brenda-<release>-<build>`
(for example `brenda-2026.1-1`). When BRENDA is updated, every file is
rebuilt and released together.

| File | Content | Used by |
|---|---|---|
| `kcat.tsv.gz`, `sa.tsv.gz` | Grouped per EC number, substrate and organism: `max`, `median`, `p75`, `n`, references | `brenda` source (`gecko_classic`, `gecko_median`) |
| `mw.tsv.gz` | Molecular weights per EC number and organism | Converting specific activities |
| `brenda_measurements.tsv.gz` | One row per wild-type measurement: kind (kcat or specific activity), EC number, substrate name and InChIKey, organism and taxonomy ID, UniProt accession, value, temperature, pH, native-substrate flag, PMID | `ec_tree`, `homologs` |
| `brenda_sequences.fasta.gz` | One sequence per UniProt accession | `homologs` (MMseqs2) |
| `brenda_substrates.tsv.gz` | Substrate name variants to InChIKey and ChEBI | Matching model metabolites to measured substrates (homologs), in MATLAB without Python |
| `kcatCombineParams.tsv` | Bias and variance per source and evidence bin, for the parameter sets K1 and K3 | `combine` |
| `README.md` | BRENDA citation and statement of changes, credits for UniProt, ChEBI, Rhea, MetaNetX, NCBI Taxonomy, PubChem and CatPred-DB, BRENDA's benefit-sharing note | |

- **One build:** geckopy's `brenda-refresh` reads the BRENDA JSON once and
  applies one set of filters to every file:
  - wild type only;
  - placeholder values dropped;
  - the outlier screen (a value more than 4 log10 units below its EC
    median, widening to the parent EC below 5 values; and a cluster below
    0.01 /s under a gap of at least 3 log10 units);
  - the corrections checked against the papers.
- Every file starts with the same header line: BRENDA release, build,
  date, licence (CC BY 4.0).
- The screen and corrections change `kcat.tsv` compared with today's
  snapshot (medians of groups with implausibly low values; groups made of
  outliers only; three corrected values). This is stated in the release
  notes of the GECKO and geckopy versions that switch.

### 3.4 DLKcat: build files and image

- **Image:** `ghcr.io/compmeteng/dlkcat-gecko:0.1`, labelled with
  `org.opencontainers.image.source=https://github.com/CompMetEng/toolbox-assets`,
  so the package is linked to toolbox-assets and managed with it. Its
  visibility is set to public (a container package has its own
  visibility, separate from the repository's).
- **Build files:** the asset `dlkcat-gecko-0.1.zip` (`DLKcat.py`,
  `Dockerfile`, `requirements.txt`, the model weights and the
  `.pickle`/`.npy` dictionaries, which GECKO keeps in git LFS today).
- **Building:** a workflow in toolbox-assets unpacks the asset, builds the
  image and pushes it. The same image tag is never pushed twice; a change
  gets a new version, as with the other assets.
- **Running does not need the asset:** `runDLKcat` and `run_dlkcat` pull
  the image and run it; DLKcat reads the copies of the files inside the
  image. Their default image name changes to the CompMetEng one.
- **The old image** `ghcr.io/sysbiochalmers/dlkcat-gecko:0.1` stays
  available for released GECKO and geckopy versions. Container packages
  cannot be transferred between organizations, so the CompMetEng image is
  pushed anew (the same image, retagged, or rebuilt from the asset).
- `src/dlkcat-gecko/` and git LFS are removed from GECKO.
  `DLKcatCurrencyMets.tsv` and `DLKcatIgnoreMets.tsv` stay in the
  repositories, since users extend them.

### 3.5 KEGG phylogenetic distances

`keggPhylDist.mat` is derived from `kegg118_taxonomy.gz`, which is already
an asset. RAVEN ships the `.mat` in its repository and GECKO reads it
from there; geckopy downloads it from RAVEN's repository. geckopy can
instead compute it with `raven_toolbox.reconstruction.kegg.phyl_dist`
from the taxonomy asset, so no separate asset is needed.

## 4. Project kcat files

### 4.1 Format

`data/kcats/<source>.tsv`, one per source. Header lines start with `#`:

```
# source: brenda
# generator: gatherKcats (GECKO 4.x) | gather_kcats (geckopy 4.x)
# parameters: {"aggregation": "median", "max_wildcard_level": 3}
# data: brenda-2026.1-1
# created: 2026-10-15T12:00:00Z
rxn	enzymes	substrate	kcat	log10_sd	source	detail	ec_code	note
```

| Column | Required | Meaning |
|---|---|---|
| `rxn` | no | Conventional reaction, `_REV` for the reverse direction. Empty: every reaction of `enzymes` (curated rows only) |
| `enzymes` | no | UniProt accessions, sorted, joined by ` + ` for a complex. Empty: every enzyme set of `rxn` (curated rows only) |
| `substrate` | no | Metabolite name as in `model.metNames` |
| `kcat` | yes | 1/s |
| `log10_sd` | no | Uncertainty, for `combine` and tuning priors |
| `source` | yes | Origin token (section 9) |
| `detail` | no | `key=value` pairs, e.g. `wc=0; origin=2` |
| `ec_code`, `note` | no | EC number used; free text, reference |

At least one of `rxn` and `enzymes` is filled; gatherers fill both.
Unknown columns are kept when reading and writing.

### 4.2 Matching to a model

`loadKcatSource(model, name)` maps rows to `model.ec` rows, for full and
light models:

1. The base of a `model.ec.rxns` entry is the entry without its trailing
   `_EXP_<n>` (full model) or leading `<nnn>_` (light model).
2. Rows with `rxn` and `enzymes` match the rows with that base and that
   enzyme set.
3. Rows with only `enzymes` match every row with that enzyme set; rows
   with only `rxn` match every row with that base.
4. Partial enzyme-set matches (at least half) are reported, not applied,
   as `applyCustomKcats` does today.
5. Rows that match nothing are reported, never dropped silently.

### 4.3 Writing

- `gatherKcats(model, source, ...)` gathers and writes
  `data/kcats/<source>.tsv` (the `saveAs` option names another file),
  and returns the kcat list.
- `saveKcatSource(model, kcatList, name)` writes any kcat list; the stable
  key is derived from the model (`ec.rxns` → base reaction,
  `ec.rxnEnzMat` → enzyme set).
- Writing never touches a file whose header says `source: custom` or
  `kind: curated`.

## 5. `gatherKcats`

One function per tool with the source as its first argument:

| Source | Gathers | Options |
|---|---|---|
| `brenda` | Fuzzy BRENDA matching (today's `fuzzyKcatMatching`) | `aggregation` (`max`, `median`, `p75`), wildcard and origin limits |
| `dlkcat` | DLKcat through Docker (today's `writeDLKcatInput`, `runDLKcat`, `readDLKcatOutput`) | `overwrite` |
| `okp` | OpenKineticsPredictor (today's submit, fetch and read functions) | `method` (one or several, one file each), `use_stored` |
| `ec_tree` | EC tree estimate from the BRENDA measurement file | `k` (shrinkage, default 5) |
| `homologs` | Identity-weighted homolog mean through MMseqs2 | identity and coverage limits, number of homologs |

- Its help lists every source with its options; an option a source does
  not take is an error.
- Today's public gatherers stay as deprecated wrappers until GECKO 5.0
  (O4).

## 6. Assignment configuration

In `model_adapter.toml`:

```toml
[kcat_assignment]
preset = "gecko_median"     # or leave out and give stages
isozymes = true             # fill from isozymes after the stages
standard = true             # standard kcat for the rest

# Instead of a preset:
[[kcat_assignment.stage]]
action = "select"           # first listed source with a value wins, per reaction
sources = ["custom", "database_top", "dlkcat", "database_bottom"]
top_origin_limit = 6
bottom_origin_limit = 6
wildcard_limit = 3
```

- A stage decides only reactions that are still open.
- `select` uses `mergeKcats` / `merge_kcats`, with the BRENDA tiers
  `database_exact`, `database_top` and `database_bottom` and plain source
  names.
- `combine` (later) combines the listed sources by uncertainty with
  `kcatCombineParams.tsv`; `level = "K1"` or `"K3"`.

Presets:

| Preset | Stages | Status |
|---|---|---|
| `gecko_classic` | custom → BRENDA exact/top → DLKcat → BRENDA bottom, aggregation `max` | Reproduces today's tutorial |
| `gecko_median` | The same with aggregation `median` | Template default now |
| `ec_tree_homologs` | custom → combine(EC tree, homologs) | Default after the BRENDA measurement release |
| `ec_tree` | custom → EC tree | Without MMseqs2 |
| `combined_all` | custom → combine(EC tree, homologs, all OKP methods) | For single kcat values, not recommended for models |

`isozymes` and `standard` apply after every preset. A level correction
(fitting σ) is a separate step, not part of assignment.

## 7. `assignKcats`

`model = assignKcats(model)` / `assign_kcats(model)`:

1. reads `[kcat_assignment]` from the adapter;
2. loads the source files the stages name and maps them to `model.ec`
   (section 4.2); a missing file is an error that names the
   `gatherKcats` call that makes it;
3. runs the stages, then the isozyme fill and the standard kcat;
4. writes `ec.kcat` and `ec.source`, and appends curated notes to
   `ec.notes`;
5. records one summary entry in the change log (E1);
6. writes `output/<model>_kcatDecisions.tsv`: per `model.ec` row the
   stage, source, value, number of competing candidates and their spread;
7. runs `applyKcatConstraints`.

`mergeKcats`, `assignKcatValues`, `getKcatAcrossIsozymes` and
`getStandardKcat` stay callable on their own.

## 8. Custom kcats

- `data/kcats/custom.tsv` in the format of section 4.1, with
  `# kind: curated`.
- `convertCustomKcats` turns `data/customKcats.tsv` into it. The
  `stoicho` column is dropped (complex stoichiometry comes from
  `applyComplexData`).
- Light models get custom kcats; today `applyCustomKcats` skips them.
- `applyCustomKcats` stays as a deprecated wrapper.

## 9. `ec.source` format

```
<origin>[ (<detail>)][ | tuned=<method>[+<method>] (prior=<value>)]

brenda (wc=0, origin=2)
custom
okp_catapro (similarity=42)
combined (sd=0.61)
brenda (wc=0, origin=2) | tuned=cmaes (prior=6.55)
okp_catapro (similarity=42) | tuned=sensitivity+cmaes (prior=5.48)
```

- **Origins:** `custom`, `brenda`, `sabio_rk`, `uniprot`, `dlkcat`,
  `okp_<method>`, `ec_median`, `ec_tree`, `homologs`, `combined`,
  `isozymes`, `standard`, `manual`.
- **Detail keys per origin:**
  - `brenda`: `wc`, `origin`
  - `okp_<method>`: `similarity`
  - `homologs`: `identity`
  - `combined`: `sd`
  - `isozymes`: `from` (the `model.ec` row copied)
  - others: none
- **Tuning methods:** `sensitivity`, `cmaes`. `prior` is the value before the
  first tuning.
- **Who writes:** assignment stages and `setKcatForReactions` (origin `manual`)
  write a new origin and drop any status. Tuning keeps the origin and
  extends the status.
- **One parser:** `parseKcatSource` / `parse_kcat_source` returns origin,
  detail, tuning methods and prior; `formatKcatSource` /
  `format_kcat_source` is its inverse. Every reader uses them:
  - prior grouping in both tuning implementations (replacing MATLAB's
    exact match and geckopy's token match);
  - `getStandardKcat` and `removeStandardKcat`;
  - `sensitivityTuning`;
  - `reviewKcatAssignment`;
  - the corrections report.
- **Uncertainty for tuning priors:** `sd` in the detail gives the
  combination's uncertainty per kcat. Tuning uses it as the prior width when
  present, and the adapter's per-group width otherwise.

## 10. Change log

### 10.1 Files and in-memory form

- In memory: a top-level model field, MATLAB `model.changeLog` and Python
  `model.change_log`, with two tables, `runs` and `changes`. Copies of the
  model carry it, and functions that rebuild `model.ec` do not drop it.
- On disk: `<model>.runs.tsv` and `<model>.changes.tsv` next to the model.
  `saveEcModel` writes them, and `loadEcModel` reads them if present.
- **Integrity checks:**
  - Each run stores `state_after`, a hash of `ec.kcat`, `ec.concs`, the
    protein usage bounds and the pool bound.
  - `loadEcModel` warns when the log is missing, or when the last
    `state_after` does not match the loaded model.
- The YAML model format does not change.

### 10.2 Columns

`runs`:

| Column | Meaning |
|---|---|
| `run` | Identifier (time and short hash) |
| `function`, `tool` | E.g. `cmaesKcatTuning`; `GECKO 4.x / MATLAB R2023b` or `geckopy 4.x` |
| `time` | ISO 8601 |
| `parameters`, `seed` | JSON; RNG seed |
| `inputs` | Data files with checksums |
| `summary` | JSON: RMSE before and after; growth before and after; counts |
| `state_before`, `state_after` | Hashes as above |

`changes`:

| Column | Meaning |
|---|---|
| `event`, `run` | Sequence number; run |
| `parameter` | `kcat`, `conc`, `conc_bound`, `pool`, `sigma`, `f`, `flux_bound` |
| `target` | `model.ec.rxns` entry, UniProt accession, reaction, or `*` for a summary |
| `rxn`, `enzymes` | Stable key for kcats |
| `old`, `new` | Values |
| `action` | Fixed vocabulary, see 8.3 |
| `reason` | `key=value` pairs for the action |
| `note` | Free text |

### 10.3 What each function records

| Function | Parameter | Action | Reason |
|---|---|---|---|
| `assignKcats` | kcat (summary) | `assign` | rows per stage and source |
| `setKcatForReactions` | kcat | `manual` | note |
| `sensitivityTuning` | kcat | `sensitivity_tuning` | steps, fold step, growth before and after, control coefficient |
| `cmaesKcatTuning` | kcat | `cmaes_tuning` | leverage share and rank, prior group and width, run RMSE before and after |
| `removeStandardKcat` | kcat | `remove_standard` | |
| `fillEnzConcs` | conc | `proteomics_fill` | data file, column or condition |
| `constrainEnzConcs` | conc_bound | `constrain` | |
| `flexibilizeEnzConcs` | conc_bound | `flexibilize` | growth reached and expected, condition, fold step, steps |
| `relaxProteomicsGreedy` | conc_bound | `relax` | growth before and after |
| `setProtPoolSize`, `updateProtPool` | pool | `pool_update` | Ptot, f, sigma |
| `sigmaFitter` | sigma | `sigma_fit` | data, fit |
| `constrainFluxData` (E2, when `params.logFluxConstraints` is true) | flux_bound | `flux_constraint` | condition, data file |

Granularity: one event per parameter per run. Iterative steps (sensitivity
tuning, flexibilization) are summarised in `reason`. CMA-ES writes one event
per changed kcat, not per generation.

### 10.4 Helpers and tools

- **Recording:** `run = startChangeRun(model, functionName, parameters,
  inputs)` and `model = recordChanges(model, run, parameter, targets, old,
  new, action, reason, note)`, on both sides.
- **Reporting:** `reportChanges(model, ...)` gives filtered views. The
  corrections table of a CMA-ES run and the flexibilization table are views.
- **Re-tuning:** `priorValues(model)` gives the pre-tuning value per kcat,
  from `ec.source` (`prior=`) or the log.
- **Undo:** `revertChanges(model, run)` writes back the `old` values of one
  run and records the revert.
- **Verification:** `checkChangeLog(model, parent)` replays the log on the
  parent model (for example the stage-2 model) and confirms that the current
  values come out.

### 10.5 Lineage

The log records the history of this model. Runs that are not applied (other
seeds, budgets) stay as reports in `output/`. A condition-specific model
saved under a new name starts with its parent's log. Parallel jobs work on
copies; only the applied result's events join the model's log.

## 11. Parity scenarios

| Scenario | Checks |
|---|---|
| `kcat_source_roundtrip` | A source file written in Python reads identically in MATLAB and back; header lines and extra columns kept |
| `kcat_source_match_full_light` | Matching rules 1 to 5 on full and light ecTestGEM, including curated rows with empty `rxn` or `enzymes` |
| `kcat_assignment_presets` | Every preset gives the same `ec.kcat`, `ec.source` and decision report in both tools |
| `kcat_assignment_classic_reproduces` | `gecko_classic` on the pre-screen BRENDA snapshot reproduces today's tutorial chain exactly, on ecTestGEM and the yeast tutorial |
| `kcat_custom_first` | Custom kcats decide first, also on a light model |
| `kcat_one_way` | After assignment, tuning and flexibilization every file in `data/kcats/` is unchanged |
| `brenda_asset_identical` | Both tools load the same BRENDA asset to the same tables |
| `ec_source_format`, `change_log_*` | As in sections 9 and 10 |
| `ec_tree_fixture`, `homologs_fixture`, `combine_fixed` | New sources on a fixture of a few hundred measurement rows |

## 12. Phases

Each phase is one PR in each affected repository, merged together once its
parity scenarios pass.

| Phase | Content | Repositories | Depends on |
|---|---|---|---|
| T1 | Transfer raven-data to `CompMetEng/toolbox-assets` (done); `fetchToolboxAsset` (RAVEN#764) and its SHA-256 check; manifest-v2 | toolbox-assets, RAVEN, raven-toolbox | — |
| T2 | BRENDA asset (screened and corrected, grouped files only), DLKcat asset and the CompMetEng DLKcat image with its build workflow; GECKO and geckopy load from the assets and run the new image; BRENDA files and git LFS removed from the repositories | toolbox-assets, geckopy (`brenda-refresh`), GECKO | T1 |
| K1 | Project kcat files: format, `saveKcatSource`, `loadKcatSource` with matching; `gatherKcats` for `brenda`, `dlkcat`, `okp`; deprecated wrappers | GECKO, geckopy | T2 |
| K2 | `[kcat_assignment]` in the TOML schema and `adapterParameters.tsv`; `assignKcats`; presets `gecko_classic` and `gecko_median`; decision report; tutorials switched | GECKO, geckopy | K1 |
| K3 | Custom kcats as a source, `convertCustomKcats` | GECKO, geckopy | K1 |
| K4 | BRENDA measurement files, sequences, substrate dictionary, combination parameters in the BRENDA asset; `gatherKcats` for `ec_tree` and `homologs` (MMseqs2); the `combine` stage; presets `ec_tree_homologs`, `ec_tree`, `combined_all`; validation on the benchmark's development split and an ecYeastGEM build; template default switched | toolbox-assets, geckopy, GECKO | K2 |
| L1 | `ec.source` format and parser (section 9) | GECKO, geckopy | — |
| L2 | Change log (section 10) | GECKO, geckopy | L1 |

T1 and L1 can start at once; L1 and L2 ship independently of the kcat work.

## 13. Open decisions

- **O2.** Combination parameter sets: ship K1 and K3, chosen by the number
  of measurements for the organism, overridable in the adapter (proposed).
- **O3.** Substrate keys for homologs: settled by shipping
  `brenda_substrates.tsv`; geckopy may also compute InChIKeys from SMILES.
- **O4.** Deprecation wrappers stay until GECKO 5.0 (proposed).
