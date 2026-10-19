# GECKO version in ecModel files

Plan for [CompMetEng/GECKO#5](https://github.com/CompMetEng/GECKO/issues/5): ecModel files record which
GECKO or geckopy release wrote them and which ecModel format they follow, and both loaders act on it.

## Current state

| | GECKO (MATLAB) | geckopy |
|---|---|---|
| Model YAML | nothing; `version` in `metaData` is the model's own version | `metaData.geckopy_version` (package version), date reset on every save |
| Kcat source files | `# generator: <function> (GECKO <version.txt>)` | `# generator: <function> (geckopy <version>)` |
| Project | `project_format = 4` in `model_adapter.toml` after `upgradeGECKOproject` | the same, after `upgrade_project` |
| Older ecModels on load | `loadEcModel` reads them unchanged | raven-toolbox flips GECKO 3 `usage_prot_*` / `prot_pool_exchange` with a warning |

RAVEN's `writeYAMLmodel` and `readYAMLmodel` already round-trip any extra `model.annotation` field as a
`metaData` key, and raven-toolbox keeps `metaData` on `model.notes`, so neither RAVEN nor raven-toolbox
needs to change. `version.txt` is updated on `main` at release; on `develop` it names the last release.

## Fields

Two `metaData` keys, written by both tools:

```yaml
- metaData: !!omap
  - id: yeastGEM_v9.1.1
  - version: 9.1.1               # the model's own version, unchanged
  - gecko_format: 4              # ecModel format: what a loader branches on
  - generator: GECKO 4.0.1       # or: geckopy 4.0.1; which release wrote the file
```

- `gecko_format` is an integer with the same numbering as `project_format`: 4 for the current layout
  (forward `usage_prot_*` and `prot_pool_exchange`, `ec-rxns` / `ec-enzymes`, `gecko_light` at the top level).
  It changes only when the file layout or the meaning of a field changes, not with every release, so
  GECKO and geckopy releases can differ without the files differing.
- `generator` replaces `geckopy_version`. It is provenance only; no code branches on it.
- A file without `gecko_format` is format 3 if its usage reactions run in the GECKO 3 direction,
  otherwise format 4 (files written by GECKO 4 before this change).

## Work

| Step | What | Where |
|---|---|---|
| V1 | `saveEcModel`: set `annotation.gecko_format = 4` and `annotation.generator = ['GECKO ' version]` before `writeYAMLmodel` (YAML only; SBML is a lossy export). `save_ec_model`: write `gecko_format` and `generator` instead of `geckopy_version`; keep the model's `date` when it has one, as RAVEN does | GECKO, geckopy |
| V2 | `loadEcModel` / `load_ec_model`: read `gecko_format` (or infer it as above); format 3 is upgraded in memory with `upgradeEcModel` / `upgrade_ec_model` and a warning naming `upgradeGECKOproject`; a format above the supported one is an error that names the release that wrote it; geckopy reads `geckopy_version` as `generator` | GECKO, geckopy |
| V3 | `upgradeEcModel` / `upgrade_ec_model` set `gecko_format = 4`; `upgradeGECKOproject` / `upgrade_project` save the upgraded models with it | GECKO, geckopy |
| V4 | Tests: save writes both keys; load of a file without them (format 3 and 4 fixtures); a newer format errors; parity scenario: a file written by one tool loads in the other with the same `gecko_format` and its own `generator` on re-save | GECKO, geckopy, raven-gecko-parity |
| V5 | Docs: the two keys in the ecModel format page and in "GECKO 3 to GECKO 4"; close GECKO#5 | GECKO, geckopy |

V1 to V3 are one PR per tool; V4 goes with them, the parity scenario in raven-gecko-parity.

## Open points

- Whether the in-memory upgrade on load (V2) should be silent apart from the warning, as in geckopy now,
  or ask the user to run `upgradeGECKOproject` first. The plan follows geckopy, so both tools behave alike.
- Kcat source files keep their `# generator:` line; a format number for them is not needed until their
  columns change.
