# Day-to-day workflow

Where the ledger fits into the work you were doing anyway.

This page is for maintainers and for anyone working on a MATLAB/Python pair on purpose. A
contributor who opens a PR against RAVEN, GECKO, raven-toolbox or geckopy needs nothing from
this page and installs nothing. The repository's own CI workflow runs `parity mirror` on the
PR and posts the result as a PR comment. For RAVEN and GECKO that workflow is
`.github/workflows/parity.yml`, rendered from `docs/shared/workflows/matlab-parity.core.yml`
in this repository.

## A bug is reported against MATLAB RAVEN

1. Reproduce and fix it in RAVEN as usual.
2. Push and open the PR. RAVEN's `parity.yml` workflow comments on the PR with what the
   change implies for raven-toolbox. To see the same answer before pushing, run this inside
   your RAVEN checkout (it needs `pip install -e ".[dev]"` from this repository once):

   ```bash
   parity mirror --since develop..HEAD
   ```

   Running `python hooks/install.py` once installs a `pre-push` git hook in each of the four
   repositories that runs this command on every `git push`, before anything is sent to
   GitHub. The hook prints information and never blocks the push, unless `PARITY_STRICT=1`
   is set in your shell.
3. If the touched function is `parity`, check whether raven-toolbox has the same bug. It
   usually does; the port followed the original closely. Open the sibling issue with the
   RAVEN issue linked, and record it on the ledger row while it is open:

   ```yaml
   - matlab: fillGaps
     python: raven_toolbox.gapfilling.connect_blocked_reactions
     status: parity
     notes: RAVEN#456 / raven-toolbox#789: same off-by-one in the weight vector.
   ```
4. If the touched function is `matlab-only` or `via-dependency`, there is nothing to mirror.
   The ledger already answered the question.
5. If the row is `unreviewed`, triage it now. A bug report is the cheapest possible moment to
   decide what a function is.

## A feature is developed in Python first

1. Build it in raven-toolbox or geckopy.
2. `parity sync` adds a row for the new export; set its status in the same PR:

   ```yaml
   - python: raven_toolbox.localization.assign_compartments
     status: matlab-pending
     issue: SysBioChalmers/RAVEN#123
     notes: Port plan: MILP via optimizeProb/getMILPParams, reusing parseScores.
   ```

   `matlab-pending` without an `issue` is a warning, not an error: it is fine to defer the
   decision, not fine to lose it.
3. When the back-port lands, flip the row to `parity` and add the MATLAB name. `parity check`
   enforces this: a `matlab-pending` row that names a MATLAB function is an error, because
   the status is now a lie.
4. Consider a scenario. Anything with a solver in it earns one.
5. Before a row is set to `parity`, or its scenario result changes, run `parity ledger-audit --results <dir>`,
   where `<dir>` holds the last scenario results (`<id>.txt` or `status/<id>.txt`, the first word being `MATCH`,
   `DIFFER` or `ERROR`). It fails when a `parity` row has a scenario recorded as `DIFFER` or `ERROR` and a note that
   does not mention a difference or name the scenario (the way to say that a scenario covering several functions differs elsewhere), and when a scenario that covers the function is not listed in the row's
   `scenarios`. The note has to list every known difference; the audit checks that it says something, not that it
   is complete, so read the scenario output as well. A row can state where it differs with `differs`, a list of
   checkpoint path prefixes that the audit checks against the differences rebuilt from `python/<id>.seed<N>.json`
   and `matlab/<id>.json` in the same directory; a ledger with `strict_differs: true` requires the list on every
   `parity` row linked to a scenario that differs and requires every difference group (the first two path
   segments) to be listed by an entry that is the group or goes deeper (see the README).

## A deliberate divergence

Record it with a reason, and the question stops coming back:

```yaml
- matlab: runDynamicFBA
  status: matlab-only
  reason: Maintained Python packages already cover dynamic FBA (dfba, reframed, mewpy).
```

`parity report` renders these as the divergence register, so the rationale lives with the
decision instead of in a pull request comment from two years ago.

## Releases

Do not force version lockstep; it produces empty releases. State the baseline instead, one
line in each CHANGELOG:

```
Parity baseline: RAVEN 2.11. Exceptions: parity/raven.yml in raven-gecko-parity.
```

Before tagging:

```bash
parity check --strict
parity ledger-audit --results <last nightly results>
parity report -o build/
```

and run the scenarios for anything that changed since the last tag.

## Working with a coding agent

The `AGENTS.md` of each code repository, rendered from `docs/shared/AGENTS.core.md`, gives a
coding agent the ledger duties above, and the `parity-ledger` skill of the `raven-gecko`
plugin in [CompMetEng/skills](https://github.com/CompMetEng/skills) has the statuses and
commands. The purpose is the same as that of the git hooks: the cross-check happens without
anyone having to remember it.

## Shared documentation and CI workflows

`docs/shared/` in this repository is the single source for each code repository's
`CONTRIBUTING.md`, `AGENTS.md`, issue and PR templates, and CI workflow files.
Every generated `AGENTS.md` has a "Cross-implementation parity" section that tells an agent
to check this ledger before closing a bug or finishing a feature.

To change any of it:

1. Edit the file under `docs/shared/`: a `*.core.*` file for content every repository
   shares, or an `appendices/<repo>.*` file for one repository's own section.
2. Render it into the sibling checkouts that `parity.toml` configures:

   ```bash
   parity docs render
   ```

   This writes files only on your machine. It does not commit or push.
3. Review the diff in each sibling checkout and open a PR there.

`parity docs check` lists which sibling files are missing or differ from the current
rendered output, without writing anything.

## Decision records

`scripts/decision_records.py` validates the records in a `decisions/` folder (statuses, names, links, and what a frozen record may change) and the `decision-records` job of `ci.yml` runs it on every pull request and on every push to `develop`. The rules of a record are in the docstring of the script and in the `decision-records` skill of the `craft` plugin of CompMetEng/skills. The script has its own tests in `tests/test_decision_records.py`. The four code repositories do not receive the script.
## Writing documentation

The prose rules for every markdown file, docstring, comment and commit or pull request text
in these repositories are in the `writing` skill of the `craft` plugin in
[CompMetEng/skills](https://github.com/CompMetEng/skills). They include the ban on em dashes
and the per-sentence checks.
