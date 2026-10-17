### What changed and why

<!-- One bullet per change, in the same style as the commit messages, e.g.
- fix: `functionName` failed to do X (fixes #123)
- feat: added option to do Y
- docs: updated the Z section
-->

### Checklist

- [ ] The PR title follows [Conventional Commits](https://www.conventionalcommits.org/)
      (`type: description`); it becomes the squash-merge commit message.
- [ ] Tests cover the new or changed behavior, and CI passes.
- [ ] I checked [raven-gecko-parity](https://github.com/CompMetEng/raven-gecko-parity)
      for parity implications of this change, and updated or referenced the ledger if needed.
- [ ] This PR targets `develop`/`devel` (squash-merge), or `main`, following the release
      process in `CONTRIBUTING.md`.
