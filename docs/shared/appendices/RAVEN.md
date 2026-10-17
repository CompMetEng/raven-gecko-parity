**RAVEN** (MATLAB toolbox)

- Every function starts with a commented header describing it and its
  parameters; check an existing function for the expected style.
- Aim for backward compatibility: avoid MATLAB syntax newer than necessary,
  to accommodate users on older MATLAB versions.
- Run `updateDocumentation()` before a pull request from `devel` to `main`.
- External binary updates (BLAST+, CD-HIT, DIAMOND, HMMER, MAFFT, WoLFPSORT):
  branch from `devel` as `chore/update-binaries`, commit each program/OS
  separately (e.g. `chore: update HMMER (Win) to 3.2.1`), update the
  matching license file and `software/version.txt`.
- Follows [semantic versioning](https://semver.org/); a pull request from
  `devel` to `main` triggers a release, with `version.txt` updated in the
  same PR.
