# Kcat strategies for different users

GECKO and geckopy gather kcats per source (`kcatsFrom<Source>`, one file per
source in `data/kcats`) and assign them in stages (`collectKcats`,
`assignKcats`, `[kcat_assignment]` in `model_adapter.toml`). This plan adds
what different users need to choose their own strategy: where the kcats come
from, and how several values for one reaction become one.

## 1. Who wants what

| User | Trusts | Wants |
|---|---|---|
| Default | Measured kcats, the model's EC numbers | Measurements where they exist, otherwise the closest estimate from measurements; predictions only where there is nothing else |
| Prediction-first | Predictors more than the model's EC annotation | Kcats from OpenKineticsPredictor (OKP) only, as a consensus over the methods |
| Mixed | Measurements and predictors | Each reaction from all evidence at once, weighted by how reliable each source is |
| Single method | One predictor (e.g. CataPro) | That method's predictions, nothing else |
| Reproducing GECKO 3 | The published tutorial | BRENDA fuzzy matching and DLKcat in the old order |
| Curator | Their own values | Their values first, whatever the rest |

## 2. Sources

| Source file | Function | Needs | Notes |
|---|---|---|---|
| `custom.tsv` | by hand, `saveKcatSource` | — | curated; never overwritten |
| `brenda.tsv` | `kcatsFromBrenda` | BRENDA asset (`kcat`, `sa`, `mw`) | fuzzy matching; origin and wildcard tiers |
| `ec_tree.tsv` | `kcatsFromECTree` | BRENDA asset (`measurements`), EC numbers | any reaction with an EC class with measurements |
| `okp_<method>.tsv` | `kcatsFromOKP` | OKP API key; sequences and substrate SMILES | one file per method; keeps the similarity of each enzyme to the method's training sequences |
| `dlkcat.tsv` | `kcatsFromDLKcat` | Docker | kept for `gecko_classic` and `gecko_median`; OKP has DLKcat too (`okp_dlkcat`) |
| `homologs.tsv` | `kcatsFromHomologs` (later) | BRENDA asset (`sequences`), MMseqs2 | |

The BRENDA files are no longer part of GECKO and geckopy: they are one
release of [toolbox-assets](https://github.com/CompMetEng/toolbox-assets),
`brenda-<release>-<build>`, fetched on first use and checked against the
manifest checksums. A project can keep its own BRENDA files in
`data/brenda` (`geckopy brenda-refresh`).

## 3. Ways to decide a reaction's kcat

A stage decides the reactions still open; stages run in order.

- **`select`**: the first listed source with a value wins (today's stage).
- **`combine`**: all listed sources with a value are averaged in log10,
  each corrected for its bias and weighted by 1 / its error variance; both
  depend on the source and its evidence for that reaction:
  - EC tree: how many measurements its EC has and at which EC level;
  - OKP method: the enzyme's similarity to the method's training sequences;
  - homologs: the identity of the nearest homolog.

  The result also has an uncertainty (`sd`), written to `ec.source` as
  `combined (sd=0.61, n=6)` and usable as the prior width for tuning.
  Bias and variance per source and evidence class are fitted once against
  measured kcats hidden from the sources (the kcat benchmark's development
  split) and ship as `kcatCombineParams.tsv` in the BRENDA asset. A source
  without parameters is an error, so a combination is never weighted by
  guesses.
- **Consensus of predictors** is `combine` over OKP methods only: methods
  that agree, and that know similar enzymes, count more.

Source groups in a stage: `okp` stands for all OKP methods set in
`[kcat_sources.okp] methods` (default: all methods OKP offers).

## 4. Presets

| Preset | Stages | For |
|---|---|---|
| `brenda_ec_tree` (default) | select(custom, BRENDA origin 1, EC tree) → combine(okp) for the rest | Default |
| `ec_tree_okp` | select(custom, BRENDA origin 1) → combine(EC tree, okp) → combine(okp) | Mixed |
| `okp` | select(custom) → combine(okp) | Prediction-first |
| `gecko_median` | select(custom, BRENDA top, DLKcat, BRENDA bottom), BRENDA median | Old order, median |
| `gecko_classic` | select(custom, BRENDA top, DLKcat, BRENDA bottom), BRENDA maximum | Reproducing GECKO 3 |

After the stages, every preset fills from isozymes and gives the rest the
standard kcat (both can be switched off). A single method is a stage, not a
preset:

```toml
[[kcat_assignment.stage]]
sources = ["custom", "okp_catapro"]
```

Without an OKP key, the `okp` stages find no files: they are reported and
left out, and those reactions get the standard kcat.

### Evidence behind the presets

From the kcat benchmark (development split; the enzyme's own measurements
and those of sequences at ≥ 95% identity hidden; differences in mean
absolute log10 error):

- The EC tree is closer than the GECKO median at origins 2 and 4 and with EC
  wildcards (−0.13 to −0.18) and as close at origins 1 and 3.
- Where the EC tree has nothing EC-specific (no EC number, or no
  measurements in its classes; 31 pairs), the combination of all predictors
  is closer than the EC tree's fallback (−0.32, −0.63 to −0.06); DLKcat alone
  is not (−0.19, not significant).
- Where there is no own-organism BRENDA match (2,672 pairs), the
  combination with predictors is closer than the EC tree (−0.06, −0.10 to
  −0.01); no single predictor is significantly closer, and DLKcat is
  further off (+0.21).
- In ecModels, the combination with predictors set high-turnover enzymes
  such as ATP synthase too low, so respiration became too expensive. The
  default therefore uses predictors only where there is no measurement-based
  value; `ec_tree_okp` and `okp` are for users who accept that.

## 5. Choosing in practice

```matlab
% model_adapter.toml: [kcat_assignment] preset = "okp"
updateKcatSources(ecModel);        % submits the OKP jobs; call again to collect
ecModel = assignKcats(ecModel);    % output/<model>_kcatDecisions.tsv shows each decision
```

`collectKcats` gives the decisions without changing the model, so presets
can be compared before one is chosen. `kcatAssignmentConfig` shows what a
preset expands to.

## 6. Work

| Step | What | Where |
|---|---|---|
| S1 | BRENDA asset `brenda-2026.1-1` (kcat, sa, mw, measurements, sequences, substrates, combination parameters) in toolbox-assets and the manifest; GECKO and geckopy fetch it and no longer ship BRENDA files | toolbox-assets, GECKO, geckopy |
| S2 | OKP keeps the similarity per prediction (`similarity=` in the detail column); `[kcat_sources.okp] methods`; `okp` source group; `updateKcatSources` submits and collects all methods | GECKO, geckopy |
| S3 | `combine` stage; combination parameters for the EC tree and every OKP method fitted on the benchmark; presets `ec_tree_okp` and `okp`; `brenda_ec_tree` with combine(okp) instead of DLKcat; `sd` in `ec.source` and the decisions | kcat-benchmark, GECKO, geckopy |
| S4 | Tutorial section "choosing a kcat strategy"; parity tests that each preset gives the same kcats in both tools | GECKO, geckopy |
| S5 | Homologs (`kcatsFromHomologs`) as a further combine source | later |
