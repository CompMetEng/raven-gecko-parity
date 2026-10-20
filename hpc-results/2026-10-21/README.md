# HPC batch 2026-10-21

The MATLAB side of the parity scenarios and the RAVEN and GECKO MATLAB tests, run on the Vera cluster with `scripts/hpc` from the `doc/hpc-handover` branch (`docs/hpc-handover.md`). `MANIFEST.txt` has the commit of every repository and the runtime versions; `STATUS.tsv` has one line per job.

| Job set | Result |
| --- | --- |
| Scenarios (64) | 57 `MATCH`, 5 `DIFFER`, 2 `ERROR` |
| RAVEN test classes (23) | no failures; `tReconstruction`, `tSolver` and `tUtils` are `INCOMPLETE` (skipped tests) |
| GECKO core tests | 98 of 98 pass |

- `DIFFER`: `model_edits_ectestgem`, `replace_metabolite_toy`, `ec_fseof_ectestgem`, `fuzzy_kcat_coefficient_ectestgem`, `ec_model_io_ectestgem`. The nightly of 2026-10-20 reports the same differences.
- `ERROR`: `min_nr_fluxes_toy` requests SCIP, which RAVEN ships for Windows only; `ec_model_full_yeastgem` extends `YeastGEMAdapter`, which GECKO `develop` no longer has.
- `kcat_chain_ectestgem` and `add_reactions_from_model_smallyeast` were rerun after a missing `results/` folder and a stalled MATLAB session; the files here are from the reruns (both `MATCH`).
- Slurm output files are not included.
