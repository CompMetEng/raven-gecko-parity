# Running the MATLAB side on a cluster

This page is for whoever sets up the cluster runs. It covers what to install, how to run everything in parallel, and what to send back.

## Why

RAVEN (MATLAB) and raven-toolbox (Python) should give the same output for the functions the parity ledger marks `parity`. The comparison needs MATLAB. The hosted nightly runs the scenarios one after another in one MATLAB session. A cluster runs every scenario and every MATLAB test class at once, and also covers what the nightly cannot:

| Cannot do in the nightly | On the cluster |
| --- | --- |
| Many scenarios in parallel | one job per scenario |
| The RAVEN function tests with a MILP solver (ftINIT, `getMinNrFluxes`, gap-filling) | one job per test class, with Gurobi on the MATLAB path |
| Genome-scale runs | a separate job with more memory and time |
| Check that the Python side gives the same output on every run | each scenario runs the Python side under three `PYTHONHASHSEED` values |

More runs are welcome. Anything that fails or differs comes back as a result file, and the person reading the results (or the assistant working in the repository) turns it into a fix, a ledger note or a new scenario. Run everything every time the repositories move.

## What you need

| Item | Detail |
| --- | --- |
| MATLAB | R2024b. The jobs call `matlab -batch`, so no display is needed. Each running job holds one MATLAB licence seat, so the number of concurrent jobs is limited by the licence (`MAX_PARALLEL`, default 8, per job array). |
| Gurobi | A licence that works on compute nodes (site, cluster or WLS), plus the Gurobi MATLAB interface (`$GUROBI_HOME/matlab`). A named-user academic licence is tied to one machine and will not authenticate on compute nodes. 10 scenarios need Gurobi on both the MATLAB and the Python side. |
| Python | 3.12 or later, in a virtual environment (see Setup). |
| git | Full clones (not shallow) of the five repositories below. |
| Network | For installing and cloning. Some scenarios and tests download a file or a binary on first use (see Before the first batch). |

Repositories, side by side under one shared directory (`$WORK`):

| Directory | Repository | Branch |
| --- | --- | --- |
| `RAVEN` | `CompMetEng/RAVEN` | `develop` |
| `GECKO` | `CompMetEng/GECKO` | `develop` |
| `raven-toolbox` | `CompMetEng/raven-toolbox` | `develop` |
| `geckopy` | `CompMetEng/geckopy` | `main` |
| `raven-gecko-parity` | `CompMetEng/raven-gecko-parity` | `doc/hpc-handover` (these scripts live there until they are merged into `develop`) |

`parity.toml` names these branches and expects the siblings next to each other. Use a fast shared filesystem for `$WORK`. `$WORK` must be an absolute path without spaces or quotes; `env.sh` refuses to run otherwise. For the GECKO tests it must also not contain the text `doc`: GECKO's installer removes its own folders from the MATLAB path when the checkout path contains it, and `run_matlab_tests.sh gecko` refuses.

## Setup

```bash
export WORK=/shared/fast/parity-work && mkdir -p $WORK && cd $WORK
for r in RAVEN GECKO raven-toolbox; do git clone -b develop https://github.com/CompMetEng/$r; done
git clone -b doc/hpc-handover https://github.com/CompMetEng/raven-gecko-parity   # the scripts are on this branch for now
git clone -b main https://github.com/CompMetEng/geckopy

python3 -m venv venv && source venv/bin/activate
pip install -e raven-gecko-parity "raven-toolbox[excel]" gurobipy
pip install --no-deps -e geckopy      # its metadata pins raven-toolbox to a git ref
python - <<'PY'                       # geckopy's other dependencies, without raven-toolbox
import tomllib
deps = tomllib.load(open("geckopy/pyproject.toml", "rb"))["project"]["dependencies"]
print("\n".join(d for d in deps if not d.strip().startswith("raven-toolbox")))
PY
# pip install what it printed

cd raven-gecko-parity/scripts/hpc
cp env.local.example.sh env.local.sh      # edit: WORK, MATLAB, GRB_LICENSE_FILE, GUROBI_HOME
export PARITY_OUT=$WORK/out/2026-10-21    # names the batch; no default, see below
./make_prefs_template.sh                  # once: a MATLAB preference file with RAVEN's solver set to glpk
./check_licences.sh                       # MATLAB starts; Gurobi solves a tiny LP from MATLAB and from Python
```

Do not copy `env.sh`; site values go in `env.local.sh`, which `env.sh` reads first (a value given on the command line wins over the example file).

`PARITY_OUT` is the folder of one batch and has no default. Give each batch a name and use the same value for `submit_slurm.sh`, the GECKO tests, reruns and `collect.sh`; a default derived from the date would send a rerun on the next day to another folder. It must be an absolute path without spaces, quotes, `$`, backslashes or any of `, % : ; =`.

`make_prefs_template.sh` matters. RAVEN stores its solver choice in the MATLAB preferences. Each job uses its own private preference folder (shared ones get corrupted by concurrent sessions), and that folder starts empty, so without the template several scenarios stop with "RAVEN solver not defined". Each job copies the template into its private folder. A scenario that needs another solver sets it itself and restores glpk afterwards.

The scenario jobs do not reset the MATLAB path (the test jobs do). A `startup.m`, personal path entries or a `pathdef.m` saved earlier by `checkRaven` or the GECKO installer on the account that runs the jobs are therefore on the path of every scenario job and can shadow RAVEN functions. `check_licences.sh` prints the user path and any `startup.m` it finds; run the jobs from an account without them.

Check the setup before submitting anything:

```bash
python -c "import raven_toolbox; print(raven_toolbox.__file__)"   # must point into $WORK/raven-toolbox
parity check                                                       # must print PASS for both pairs
```

If the first command points into `site-packages`, reinstall raven-toolbox from its checkout.

## Before the first batch

1. `./check_licences.sh` (see Setup). It must end with `all checks passed`.
2. Smoke test one scenario:

   ```bash
   ./run_scenario.sh elemental_balance_smallyeast
   ```

   It prints `elemental_balance_smallyeast MATCH`, `DIFFER` or `ERROR`. `DIFFER` is a valid answer; `ERROR` means a run failed, and `$PARITY_OUT/logs/elemental_balance_smallyeast.log` says which side.
3. Smoke test a scenario that needs Gurobi (`task_checking_smallyeast`) and one that needs BLAST (`homology_chain`).
4. Fetch downloads once, by hand, so concurrent jobs do not race: RAVEN downloads BLAST and DIAMOND into its own checkout on first use through one fixed temporary file, with no lock. Running `homology_chain` (step 3) and `./run_matlab_tests.sh raven tReconstruction` once puts the binaries in place. Jobs that need a file from the network (UniProt, KEGG) fail on nodes without internet; fetch those from a login node first or accept `ERROR` for them.

## Full run

```bash
./submit_slurm.sh all              # scenarios, then the RAVEN tests, then the GECKO tests; each starts when the one before ended

# or one part at a time:
./submit_slurm.sh scenarios [raven|gecko]   # 64 scenarios (49 RAVEN, 15 GECKO), one array task each
./submit_slurm.sh tests                      # every RAVEN test class t*.m except tBinaries and tInstallation, one task each
./submit_slurm.sh gecko-tests                # all *Tests.m files of GECKO's test/unit_tests (CI runs only geckoCoreFunctionTests)
```

The GECKO tests run after the scenarios because GECKO's installer saves the MATLAB path and the tests write into shared GECKO fixture folders. `tInstallation` is left out because it needs the network and `checkRaven` saves the MATLAB path; run `./run_matlab_tests.sh raven tInstallation` by hand on a quiet system if you want it. Defaults: 2 CPUs, 16 GB and 2 h per scenario; `ec_model_full_yeastgem` gets 4 CPUs, 64 GB and 12 h. `MAX_PARALLEL` (default 8) caps the concurrent tasks of one array, and every task holds a MATLAB licence seat. `all` keeps the scenario, test and GECKO jobs apart in time (`--dependency=afterany`), so at most `MAX_PARALLEL` plus the one genome-scale scenario run together; separate `scenarios` and `tests` submissions can overlap and double that. The scripts are scheduler-neutral: the unit of work is one call to `run_scenario.sh <id>`, `run_matlab_tests.sh raven <class>` or `run_matlab_tests.sh gecko`. Translate `submit_slurm.sh` for PBS or LSF if needed; keep `PARITY_OUT` identical for every job of a batch (`submit_slurm.sh` exports it).

Do not start a second batch of the same scenario while the first is running. `run_scenario.sh` refuses (`ERROR (another run of this scenario is active)`): the context and result files of a scenario live in its folder of the shared checkout.

## What each job does

`run_scenario.sh <id>` runs the Python side three times (`PYTHONHASHSEED` 0, 1, 2), runs the MATLAB side once through `parity run --impl matlab`, then compares:

1. each Python result with the MATLAB result;
2. the Python results with each other (a difference here means the Python output depends on the interpreter's hash seed).

| Status in `status/<id>.txt` | Meaning |
| --- | --- |
| `RUNNING` | the job is alive, or it was killed before it could report (also what a time-limit kill looks like if the signal could not be handled) |
| `MATCH` | all comparisons agree within the scenario's tolerance |
| `DIFFER` | a comparison failed; the log lists the differences |
| `ERROR` | a run or the comparison itself failed (licence, missing file, crash, kill); the log has the error |

The MATLAB step of a scenario is limited to `MATLAB_TIMEOUT` seconds per attempt (default 1800; `MATLAB_TIMEOUT_HEAVY`, default 18000, for the genome-scale scenario) and runs in its own process group, so a session that hangs is killed. A scenario that timed out is tried once more; a second timeout is an `ERROR` (`MATLAB TIMED OUT` in the log). One run in ten of a 40 s scenario stalled in the first batch, so a single retry is expected to be enough.

`run_matlab_tests.sh` writes `status/test-<name>.txt`:

| Status | Meaning |
| --- | --- |
| `PASS` | all tests ran and passed |
| `FAIL` | at least one test failed |
| `INCOMPLETE` | no failure, but some tests were filtered or skipped (counts in the line) |
| `ERROR` | MATLAB failed before the tests finished, or no tests ran |

A trailing `NO-GUROBI` means the job could not solve a tiny LP with `gurobi()` before running the tests: the interface is not on the path, there is no licence, or no session was free (`GUROBI_ERROR=` in `logs/<name>.log` gives the reason). Tests that need a MILP solver switch to Gurobi only when `gurobi()` is on the path and otherwise fall back to GLPK or SCIP without reporting a skip. A `PASS NO-GUROBI` therefore means the MILP tests ran on the fallback, and a `FAIL NO-GUROBI` may be a licence problem and not a RAVEN failure. `INCOMPLETE` does not tell you which solver was used.

Files, under `$PARITY_OUT`:

| Path | Content |
| --- | --- |
| `python/<id>.seed<N>.json`, `matlab/<id>.json` | result documents |
| `logs/<id>.log`, `logs/<name>.log` | every command and the comparison output; the MATLAB output of test jobs |
| `status/` | one status file per job |
| `matlab-tests/<name>.csv` | one row per MATLAB test |
| `lists/` | the id lists of each submission |
| `slurm/` | Slurm output files |

## Sending the results back

```bash
./collect.sh [batch-folder]       # default $PARITY_OUT; writes <batch-folder>.tar.gz
```

It prints the count per status, lists submitted scenarios that never reported (`MISSING.txt`) and warns about jobs still marked `RUNNING`. Collect after the queue is empty; a partial archive says so. The tarball contains the results, logs, test tables, `STATUS.tsv` and `MANIFEST.txt` (the commit of each repository, the MATLAB and Python runtimes, the Gurobi and cobrapy versions). `collect.sh` does not start MATLAB. Give the tarball to whoever works on the repository, or commit the unpacked folder to a branch named `hpc-results/<date>` in this repository. The commit hashes in `MANIFEST.txt` make the results comparable later.

## What happens with the results

| Result | Next step |
| --- | --- |
| `DIFFER` | Read the log. Decide whether it is a bug in RAVEN, a bug in raven-toolbox, a difference in the scenario, or an intended difference. Fix it, or record it in the ledger row's `notes`. |
| `ERROR` | A harness or environment problem. Fix the scenario or the cluster set-up and rerun that id. |
| Python result depends on the hash seed | A reproducibility bug in raven-toolbox. |
| A MATLAB test fails | Report it with the test name; it is a RAVEN or GECKO problem, not a parity problem. |

A written review of about 130 suspected differences between RAVEN and raven-toolbox exists. The scenarios already in this repository cover few of them. New scenarios are written for those cases, and they run through exactly the same scripts, so a rerun after each batch of new scenarios is the expected loop.

## Known limits

- `run_scenario.sh` takes a file lock (`flock`) in the scenario folder so two runs of one scenario cannot share its context and result files. `check_licences.sh` tests that `flock` works on the shared filesystem; on a filesystem mounted without lock support the jobs report `ERROR (flock failed ...)`. Put `$WORK` on a filesystem that supports locks.

- Scenarios that do not declare a solver compare whatever each side picks. On the cluster the MATLAB side uses GLPK (from the preference template) and the Python side uses Gurobi when `gurobipy` is installed, because cobrapy prefers it. For an LP with several optimal solutions the two can legitimately differ. A `DIFFER` in a scenario whose `scenario.yml` has no `matlab_solver` and `python_solver` should be rerun with the Python side on GLPK before it is treated as a finding; the better fix is to declare the solver in the scenario.
- Scenario MATLAB sessions use the MATLAB path of the account (see Setup).
- Test classes that call the network or download tools can fail on compute nodes without internet; those failures say nothing about parity.

## Troubleshooting

| Symptom | Likely cause |
| --- | --- |
| `env.sh:` error at the start of every script | read the message; paths with spaces or quotes, a missing virtualenv or a wrong `PARITY_ROOT` |
| `License Manager Error` from MATLAB | MATLAB licence not reachable from compute nodes. Load the site's MATLAB module before submitting, or set `MATLAB` to the full path. |
| Gurobi fails only in some jobs | the licence limits concurrent sessions; lower `MAX_PARALLEL`, or run the scenario and test arrays one after the other. |
| "RAVEN solver not defined" | `make_prefs_template.sh` was not run, or the template was made with another `WORK`. |
| Many jobs fail in test-class setup, or MATLAB preference errors appear | a MATLAB preference folder shared between jobs; each job must set its own `MATLAB_PREFDIR` (the scripts do). |
| "Cannot find RAVEN Toolbox" in the GECKO job | GECKO's installer found no RAVEN on the path, or `$WORK` contains `doc` |
| `ERROR` with "no MATLAB result" | MATLAB ended without writing the file; read the end of the log. |
| `ERROR` with "Unable to encode sparse objects" | a scenario harness bug; report the id. |
| Every RAVEN scenario `ERROR` with file not found | `$WORK/RAVEN` is not the `develop` layout (`tutorials/data/models/smallYeast.yml` must exist). |
| `ERROR (another run of this scenario is active)` | an earlier job of the same id is still running; the lock is released when its process ends. Other `ERROR (cannot open ... .run.lock)` or `flock failed` messages mean the checkout or filesystem does not allow it. |
| Status stuck at `RUNNING` | the job is still alive, or the node was lost; rerun that id after `scancel`. |
| Task killed at the time limit | raise `--time` for that id and rerun it with `./run_scenario.sh <id>`. |

## Not tested

The shell scripts were tested on a machine without MATLAB or Slurm, with stand-ins for MATLAB and `sbatch`. That covered the `MATCH`, `DIFFER` and `ERROR` paths of a scenario, a kill at the time limit, the per-scenario lock, the hash-seed comparison, the test-suite status logic, the Slurm list and array construction, `collect.sh` on a partial batch and the checks in `env.sh`. Also covered: a fake `sbatch` that runs the array task from a copy of the script outside the scripts folder, the way `slurmd` does, and the `all` dependency chain. Not run: the MATLAB statements themselves (the test job file and the `restoredefaultpath` and `addpath` sequence), whether `writetable` accepts the test-result table once the `Details` column is removed, `make_prefs_template.sh`, whether `MATLAB_PREFDIR` is honoured by your MATLAB installation, the Gurobi set-up, and real `sbatch` behaviour (environment export and signal delivery at the time limit). The first smoke test is the real test of them.

## Optional extras

| Run | Command | Needs |
| --- | --- | --- |
| raven-toolbox parity tests with full Gurobi | `cd $WORK/raven-toolbox && pytest -m parity --timeout=600` | Gurobi without the size limit of the pip licence |
| geckopy tests | `cd $WORK/geckopy && pytest -q` | the venv above |
| Hosted nightly | GitHub Actions, workflow `nightly`, run with `force` | nothing on the cluster |
