**GECKO** (MATLAB toolbox)

- Every new script starts with a commented section describing it and its
  inputs/outputs; check an existing function for the expected style.
- New functions need a matching test in
  `tests/unit_tests/geckoCoreFunctionTests.m`, and the CI run against it
  must pass.
- Follows [semantic versioning](https://semver.org/); a pull request from
  `develop` to `main` triggers a release, with `version.txt` updated as part
  of it.
