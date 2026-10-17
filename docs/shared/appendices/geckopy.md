**geckopy** (Python)

- Lint and format with `ruff check` / `ruff format` (configured in
  `pyproject.toml`).
- Tests run with `pytest`; real-scale smoke tests are opt-in via
  `pytest -m smoke` and are skipped by default.
- Versioning follows [PEP 440](https://peps.python.org/pep-0440/); document
  changes in `CHANGELOG.md` using the
  [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) format.
