"""`parity seeds`: the Python side of a scenario must not depend on the interpreter's hash seed."""

from __future__ import annotations

from pathlib import Path

import pytest

from parity.cli import main

CONFIG = """
[pairs.demo]
ledger = "ledgers/demo.yml"

[pairs.demo.matlab]
repo = "org/DEMO"
ref = "develop3"
path = "DEMO"

[pairs.demo.python]
repo = "org/demopy"
package = "demopy"
path = "demopy"
"""

NAMES = "alpha bravo charlie delta echo foxtrot golf hotel india juliet kilo lima mike".split()

STABLE = f"def run(ctx):\n    return {{'order': sorted({NAMES!r})}}\n"
UNSTABLE = f"def run(ctx):\n    return {{'order': list({set(NAMES)!r})}}\n".replace(
    repr(set(NAMES)), "set(" + repr(NAMES) + ")"
)


@pytest.fixture
def workspace(tmp_path: Path, monkeypatch) -> Path:
    (tmp_path / "parity.toml").write_text(CONFIG, encoding="utf-8")
    (tmp_path / "DEMO").mkdir()
    (tmp_path / "demopy").mkdir()
    for name, body in (("stable_demo", STABLE), ("unstable_demo", UNSTABLE)):
        folder = tmp_path / "scenarios" / name
        folder.mkdir(parents=True)
        (folder / "scenario.yml").write_text(f"id: {name}\npair: demo\ndescription: d\n", encoding="utf-8")
        (folder / "run.py").write_text(body, encoding="utf-8")
    monkeypatch.delenv("PARITY_CONFIG", raising=False)
    monkeypatch.chdir(tmp_path)
    return tmp_path


def test_a_seed_independent_scenario_passes(workspace: Path, capsys):
    assert main(["seeds", "stable_demo", "--seeds", "0", "1", "2"]) == 0
    assert "same result under seeds 0, 1, 2" in capsys.readouterr().out


def test_a_seed_dependent_scenario_is_reported(workspace: Path, capsys):
    assert main(["seeds", "unstable_demo", "--seeds", "0", "1", "2", "3"]) == 1
    out = capsys.readouterr().out
    assert "unstable_demo: seed" in out and "differs from seed 0" in out
