# Agent instructions

These apply to any AI coding agent working in this repository. Add repository-specific
instructions here (or in the appendix below), not in a tool-specific file that other
contributors' agents do not read.

## No AI attribution

Never mention an AI assistant, model, or vendor anywhere this repository is
visible to others: commit messages, pull request titles and bodies, code,
comments, docstrings, docs, branch names, tags, file or folder names. Write
plain, content-only text. A human reviewed the change and is responsible
for it.

## Commits and pull requests

Follow `CONTRIBUTING.md`: [Conventional Commits](https://www.conventionalcommits.org/)
(`type: description`) for commit messages and pull request titles, and
`feat/`, `fix/`, `docs/`, `chore/`, `refactor/`, `style/`, `test/` branch
prefixes. PR titles are lint-checked in CI, so set the type correctly the first
time instead of relying on a reviewer to catch it.

## Cross-implementation parity

This repo is one half of a MATLAB/Python pair tracked in
[raven-gecko-parity](https://github.com/CompMetEng/raven-gecko-parity).
Before closing a bug or finishing a feature, check
`raven-gecko-parity/ledgers/<pair>.yml` for the functions you touched:

- `parity` -> the same bug or feature probably applies to the other
  implementation. Say so explicitly in your summary, and open or reference a
  sibling issue.
- `matlab-only` / `python-only` / `via-dependency` -> deliberate; do not
  mirror.
- `unreviewed` -> triage it as part of this change.

New public functions need a ledger row in the same pull request;
`parity check` fails without one. Run `parity mirror --since <base>..HEAD`
from the raven-gecko-parity checkout to see what a change implies.

## Decision records

If the repository has a `decisions/` folder, read the records that concern the code
you change before you change it. A record with status `accepted` binds. `proposed`
is context, `rejected` states what was considered and not chosen, and `superseded`
is history that names its successor. When a task conflicts with an accepted record,
stop, say so in your summary, and propose changing the task or adding a record that
supersedes it. An accepted record is edited to mark it superseded and for no other
reason.
`python3 .github/scripts/decision_records.py check` validates the folder.

## Repository specifics

{{APPENDIX}}
