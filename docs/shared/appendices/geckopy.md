**geckopy** (Python)

Integration branch `main`; pull requests target it. This is the Python implementation of
GECKO, built on raven-toolbox, and the ledger is `ledgers/gecko.yml` in raven-gecko-parity.

### Layout

- Source is `src/geckopy/` (`adapter`, `databases`, `ec_model`, `gather_kcats`,
  `get_enzyme_data`, `kcat_tuning`, `limit_proteins`, `utilities`), with `cli.py` providing
  the `geckopy` command.
- Tests are `tests/test_*.py`, with data in `tests/data/`.
- ecModel YAML reading and writing is in `raven_toolbox.io.ec_data`, not in geckopy.

### Code

- Function names are snake_case ports of GECKO functions (`applyKcatConstraints` becomes
  `apply_kcat_constraints`).
- The public API is flat: every public name is re-exported from `geckopy/__init__.py` and
  listed in its `__all__`. A new name needs a ledger row.
- Docstrings are NumPy style.
- Versions follow GECKO's major version (4.x) and [PEP 440](https://peps.python.org/pep-0440/).
  `CHANGELOG.md` uses the [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) format,
  and entries for added functions start with "New:".

### Checks

```bash
pip install -e ".[dev]"
ruff check src tests
pytest -q
```

- ruff uses line length 100 and rules `E,F,W`. `ruff format` is not enforced.
- pytest excludes the `smoke` marker by default. `pytest -m smoke` runs tests that need
  external models such as Human-GEM.
- A change to geckopy and raven-toolbox together needs both checkouts installed in a
  specific order; the `python-tests` skill of the `raven-gecko` plugin describes it.
