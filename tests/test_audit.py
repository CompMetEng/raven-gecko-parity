"""`parity ledger-audit`: a parity row has to admit the differences its scenarios show."""

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


def ledger(notes: str = "", scenarios: str = "      - toy_demo\n") -> str:
    return (
        "version: 1\npair: demo\nmatlab: org/DEMO\npython: org/demopy\nentries:\n"
        "  - matlab: doThing\n    python: demopy.do_thing\n    status: parity\n"
        f"    scenarios:\n{scenarios}" + (f"    notes: '{notes}'\n" if notes else "")
    )


@pytest.fixture
def workspace(tmp_path: Path, monkeypatch) -> Path:
    (tmp_path / "parity.toml").write_text(CONFIG, encoding="utf-8")
    (tmp_path / "ledgers").mkdir()
    folder = tmp_path / "scenarios" / "toy_demo"
    folder.mkdir(parents=True)
    (folder / "scenario.yml").write_text(
        "id: toy_demo\npair: demo\ndescription: d\nentries: [doThing]\n", encoding="utf-8"
    )
    (tmp_path / "results" / "status").mkdir(parents=True)
    monkeypatch.delenv("PARITY_CONFIG", raising=False)
    monkeypatch.chdir(tmp_path)
    return tmp_path


def write(workspace: Path, text: str, verdict: str | None = None) -> None:
    (workspace / "ledgers" / "demo.yml").write_text(text, encoding="utf-8")
    if verdict:
        (workspace / "results" / "status" / "toy_demo.txt").write_text(verdict + "\n", encoding="utf-8")


def test_a_match_needs_no_note(workspace: Path, capsys):
    write(workspace, ledger(), "MATCH")
    assert main(["ledger-audit", "--results", "results"]) == 0
    assert "ledger audit clean" in capsys.readouterr().out


def test_a_difference_without_a_note_fails(workspace: Path, capsys):
    write(workspace, ledger(), "DIFFER")
    assert main(["ledger-audit", "--results", "results"]) == 1
    assert "recorded as DIFFER" in capsys.readouterr().out


def test_a_difference_named_in_the_note_passes(workspace: Path):
    write(workspace, ledger("DIVERGENCE: the defaults are not the same."), "DIFFER")
    assert main(["ledger-audit", "--results", "results"]) == 0


def test_a_note_that_names_the_scenario_passes(workspace: Path):
    write(workspace, ledger("toy_demo: no difference found for this function."), "DIFFER")
    assert main(["ledger-audit", "--results", "results"]) == 0


def test_an_error_needs_a_note_too(workspace: Path):
    write(workspace, ledger("Both sides agree."), "ERROR")
    assert main(["ledger-audit", "--results", "results"]) == 1


def test_a_scenario_that_covers_the_row_must_be_listed(workspace: Path, capsys):
    write(workspace, ledger(scenarios="      - other_demo\n"), "MATCH")
    assert main(["ledger-audit", "--results", "results"]) == 1
    assert "is not listed" in capsys.readouterr().out


def test_without_results_only_the_listing_is_checked(workspace: Path):
    write(workspace, ledger())
    assert main(["ledger-audit"]) == 0


def test_rows_that_are_not_parity_are_ignored(workspace: Path):
    write(workspace, ledger().replace("status: parity", "status: python-pending"), "DIFFER")
    assert main(["ledger-audit", "--results", "results"]) == 0
