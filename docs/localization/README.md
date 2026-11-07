# Compartment assignment and confidence tracking — evidence base

The measurement and design record behind raven-toolbox's compartment-assignment
method (`assign_compartments`, `predict_localization`) and its per-reaction
confidence tracking.

These documents live here rather than on
[raven-docs](https://github.com/CompMetEng/raven-docs) because both features are
still in development and neither is at parity with MATLAB RAVEN. The
user-facing site documents current, settled behaviour in both toolboxes; a
benchmark of a method that is still changing describes a moving target, and
would be read as a statement about what the toolbox does today. When either
feature settles, its evidence moves back alongside the guide page that
describes it.

The numbers are as measured when each document was written. Where the current raven-toolbox code differs from a document, the differences are listed under [Differences from the code](#differences-from-the-code).

## Using it today

[usage.md](usage.md) is the user-facing description of `predict_localization` and
`assign_compartments`: loading predictor scores, fusing several sources, and what the
options do. It moved here with the rest when raven-toolbox retired its documentation
tree, and returns to the site when the feature settles.

## Method design

| Document | What it covers |
|---|---|
| [localization_redesign.md](localization_redesign.md) | The flux-free placement MILP plus materialised-FBA certification that `assign_compartments` is built on, and why placement and certification are separated |
| [multi_localization_design.md](multi_localization_design.md) | Why a naive reaction-level multi-localisation is unsound, and the flux-activity coupling that was designed for it. The code uses the option `multi_localize` (a second compartment kept after a loopless FVA), not the coupling |
| [transport_evidence_scoring.md](transport_evidence_scoring.md) | Replacing the blanket inter-compartment transport penalty with a per-transport, evidence-aware cost, planned for both RAVEN and raven-toolbox, which is why it sits here rather than in either. `predict_localization` takes the resulting cost; `assign_compartments` accepts `transport_cost` and does not use it |
| [confidence_tracking.md](confidence_tracking.md) | Per-reaction, multi-facet confidence scoring: the data model, the YAML and SBML round-trip, and what each facet measures |
| [curation_priority_signals.md](curation_priority_signals.md) | A catalogue of signals for ranking which assignments need manual review. Seven of its signals are implemented in `curation_priority` (section 4); the others are a design menu that was not verified against the implementation |

## Predictor benchmarks

How well DeepLoc's predictions match curated models, across four organisms
chosen for their distance from its training data.

| Document | Organism / question |
|---|---|
| [deeploc_yeast_benchmark.md](deeploc_yeast_benchmark.md) | yeast-GEM — the reference case |
| [deeploc_humangem_benchmark.md](deeploc_humangem_benchmark.md) | Human-GEM — a positive control, and the circularity it exposes |
| [deeploc_aracore_benchmark.md](deeploc_aracore_benchmark.md) | AraCore (*Arabidopsis*) — cross-kingdom, fully independent |
| [deeploc_icre1355_benchmark.md](deeploc_icre1355_benchmark.md) | iCre1355 (*Chlamydomonas*) — the most training-distant test |
| [deeploc_normalisation_benchmark.md](deeploc_normalisation_benchmark.md) | Whether raw predictor probabilities beat rescaled ones |
| [localization_finetuning.md](localization_finetuning.md) | Tuning the DeepLoc-loading hyperparameters against curated yeast-GEM |

## Method validation

| Document | What it establishes |
|---|---|
| [yeast_validation.md](yeast_validation.md) | Three yeast checks of certified assignment, including recovery of curated yeast-GEM |
| [multiorganism_validation.md](multiorganism_validation.md) | Whether the method generalises at genome scale across four kingdoms |
| [yeast_localization_benchmark.md](yeast_localization_benchmark.md) | End-to-end run — model, scoring, MILP — with a predictor-noise sweep |
| [predictlocalization_comparison.md](predictlocalization_comparison.md) | Head-to-head against RAVEN's `predictLocalization`, which optimises the same objective by a different method |

## Comparison with CarveFungi

| Document | What it covers |
|---|---|
| [carvefungi_analysis.md](carvefungi_analysis.md) | How CarveFungi's compartment assignment works, and where the two approaches diverge |
| [carvefungi_milp_benchmark.md](carvefungi_milp_benchmark.md) | This transport cost run on CarveFungi's own carve MILP, as a like-for-like test |

## Parameters

| Document | What it covers |
|---|---|
| [localization_parameters.md](localization_parameters.md) | Measured defaults for `predict_localization`, from the 2026-06-20 parameter campaign |

## Differences from the code

- **Transport cost in `assign_compartments`.** The parameter `transport_cost` is accepted and not used: the placement master is flux-free and score-only, and `minimize_transports` prunes transports by flux, not by cost. `predict_localization` uses the cost. [yeast_validation.md](yeast_validation.md), [localization_redesign.md](localization_redesign.md) and [transport_evidence_scoring.md](transport_evidence_scoring.md) carry a note where they describe the cost as an input of `assign_compartments`.
- **Multi-localisation.** `assign_compartments` has the options `multi_localize`, `multi_localize_threshold`, `multi_localize_eps` and `multi_localize_loopless`. A second compartment is proposed where a gene has a DeepLoc score at or above the threshold and is kept only when a loopless FVA on the materialised model shows a flux of at least `multi_localize_eps`. The flux-activity coupling of [multi_localization_design.md](multi_localization_design.md) is not in the current code, and there is no option `multi_localization`.

## Numbers that differ between documents

The documents report different values for the same yeast-GEM replication. The raw results (`.research_tmp/*.json`) are not in the repository, so the documents alone do not show which counting rule or run each value follows.

| Quantity | Values |
|---|---|
| Transports added | 1261 full-connect and 983 FVA-usable ([localization_redesign.md](localization_redesign.md)); 967 ([yeast_validation.md](yeast_validation.md) section 1, [multiorganism_validation.md](multiorganism_validation.md)); 980 (yeast_validation.md section 3) |
| Wall time of `assign_compartments` on yeast-GEM | about 23 s (localization_redesign.md); about 140 s (multiorganism_validation.md) |
| AraCore version | v2.0 ([deeploc_aracore_benchmark.md](deeploc_aracore_benchmark.md)); v2.1 (multiorganism_validation.md and `scripts/benchmark_certified_multiorg.py`), both with 585 reactions |
