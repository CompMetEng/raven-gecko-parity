# Contributing

Thank you for contributing! Anybody is welcome to contribute, whether by
opening issues or by sending pull requests with corrections or additions.

## Reporting issues

Open an issue for a missing feature, a bug or unexpected result, missing
documentation, or any other feedback. Before opening one:

- Search existing issues for the same report.
- Reproduce it against the current default branch.
- Include what's needed to understand it: inputs, versions, error output.

If you're unsure whether something is worth an issue, ask first in
Discussions.

## Contributing changes

Consider opening an issue and assigning it to yourself before starting
non-trivial work, to avoid duplicated effort and surface any in-progress
work you could build on. Issues labelled `help wanted` are open to anyone.

See the repository's README for local development setup.

### Branching model

- The default development branch (`develop`/`devel`) is what every pull
  request is based on and targets.
- `main` carries only tested, released code; it is not committed to
  directly.
- Name your branch after the change, prefixed with its type:
  `feat/`, `fix/`, `docs/`, `chore/`, `refactor/`, `style/`, `test/`.
  Examples: `fix/gap-filling-bounds`, `feat/parallel-solve`.

### Commit messages

Use [Conventional Commits](https://www.conventionalcommits.org/):

```
type: brief description
```

or, with an optional scope:

```
type(scope): brief description
```

| type       | when                                             |
|------------|---------------------------------------------------|
| `feat`     | a new feature or capability                        |
| `fix`      | a bugfix                                           |
| `docs`     | documentation only                                 |
| `style`    | formatting/whitespace, no behavior change          |
| `refactor` | restructuring code without changing behavior       |
| `test`     | adding or adjusting tests                          |
| `build`    | build system, packaging, or dependency changes     |
| `ci`       | CI configuration changes                           |
| `chore`    | anything else: data files, housekeeping, binaries  |

Examples:

```
feat: add YAML export for enzyme-constrained models
fix: optimizeProb parsing of Gurobi status codes
chore: update KEGG reference data to version 83.0
```

A longer explanation is welcome in the commit description. Keep individual
commits reviewable; if a change needs many commits to explain itself,
consider splitting it into separate pull requests instead.

### Pull requests

- Keep pull requests small and focused: one feature or one bug fix.
  A few small, related fixes may be combined; avoid bundling unrelated work.
- The PR author and the approving/merging reviewer cannot be the same
  person.
- At least one approval is required, and CI must pass, before merging.
- PR titles follow the same `type: description` format as commit messages
  above. CI checks this and blocks the merge on a mismatch.
- Wait at least a day before merging, to give other developers a chance to
  look.

### Reviewing pull requests

- **Functionality**: CI passes; new or changed behavior has a test.
- **Documentation**: the reasoning for the change is stated, in the PR
  description or as comments.
- **Reproducibility**: any added or changed scripts do what they claim when
  run.
- **Style**: consistent with the surrounding code; new data files are easy
  to understand.
- Be explicit rather than vague: "please include X here" rather than "X
  could be included here."
- Keep the discussion in the pull request itself rather than in private
  messages.

### Issues and labels

- Type: `bug`, `enhancement`, `docs`.
- Status: `triage`, `confirmed`, `blocked`.
- Entry point: `good first issue`, `help wanted`.
- Reference `Closes #N` / `Fixes #N` in a pull request description to close
  the issue automatically on merge.

## Repository specifics

{{APPENDIX}}
