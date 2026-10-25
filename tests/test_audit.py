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


# ---- structured `differs` lists -------------------------------------------------------


def write_runs(workspace: Path, python: dict, matlab: dict) -> None:
    """Result files of a run: one Python seed and the MATLAB side, as ``collect.sh`` lays them out."""
    import json

    for folder, name, results in (("python", "toy_demo.seed0.json", python), ("matlab", "toy_demo.json", matlab)):
        target = workspace / "results" / folder
        target.mkdir(parents=True, exist_ok=True)
        (target / name).write_text(json.dumps({"scenario": "toy_demo", "results": results}), encoding="utf-8")


def keyed_ledger(differs: str | None, notes: str = "", strict: bool = False, extra_row: str = "") -> str:
    head = "version: 1\npair: demo\nmatlab: org/DEMO\npython: org/demopy\n" + ("strict_differs: true\n" if strict else "")
    row = "  - matlab: doThing\n    python: demopy.do_thing\n    status: parity\n    scenarios:\n      - toy_demo\n"
    if differs is not None:
        row += f"    differs: {differs}\n"
    if notes:
        row += f"    notes: '{notes}'\n"
    return head + "entries:\n" + row + extra_row


DIFFERENT = ({"written": {"lines": [1, 2, 3]}, "kept": 1}, {"written": {"lines": [1, 2]}, "kept": 1})


def test_a_listed_prefix_that_matches_a_difference_passes(workspace: Path):
    write_runs(workspace, *DIFFERENT)
    write(workspace, keyed_ledger("[written.lines]"), "DIFFER")
    assert main(["ledger-audit", "--results", "results"]) == 0


def test_a_shorter_prefix_matches_deeper_differences_and_indices_are_ignored(workspace: Path):
    write_runs(workspace, {"a": {"rows": [{"x": 1}, {"x": 2}]}}, {"a": {"rows": [{"x": 1}, {"x": 3}]}})
    write(workspace, keyed_ledger("[a.rows.x]"), "DIFFER")
    assert main(["ledger-audit", "--results", "results"]) == 0
    write(workspace, keyed_ledger("[a]"), "DIFFER")
    assert main(["ledger-audit", "--results", "results"]) == 0


def test_a_stale_prefix_fails(workspace: Path, capsys):
    write_runs(workspace, *DIFFERENT)
    write(workspace, keyed_ledger("[written.lines, kept]"), "DIFFER")
    assert main(["ledger-audit", "--results", "results"]) == 1
    assert "'kept' matches no difference" in capsys.readouterr().out


def test_a_prefix_for_a_scenario_that_now_matches_fails(workspace: Path, capsys):
    write_runs(workspace, {"kept": 1}, {"kept": 1})
    write(workspace, keyed_ledger("[kept]"), "MATCH")
    assert main(["ledger-audit", "--results", "results"]) == 1
    assert "matches no difference" in capsys.readouterr().out


def test_a_row_with_a_differs_list_is_not_held_to_the_note_rule(workspace: Path):
    write_runs(workspace, *DIFFERENT)
    write(workspace, keyed_ledger("[written]"), "DIFFER")
    assert main(["ledger-audit", "--results", "results"]) == 0


def test_a_ledger_without_the_keys_is_audited_as_before(workspace: Path, capsys):
    write_runs(workspace, *DIFFERENT)
    write(workspace, keyed_ledger(None), "DIFFER")
    assert main(["ledger-audit", "--results", "results"]) == 1
    assert "neither names a difference" in capsys.readouterr().out


def test_strict_needs_a_differs_list_on_a_row_linked_to_a_difference(workspace: Path, capsys):
    write_runs(workspace, *DIFFERENT)
    write(workspace, keyed_ledger(None, notes="Both sides differ somewhere.", strict=True), "DIFFER")
    assert main(["ledger-audit", "--results", "results"]) == 1
    assert "has no differs list" in capsys.readouterr().out


def test_strict_does_not_accept_a_note_in_place_of_the_list(workspace: Path):
    write_runs(workspace, *DIFFERENT)
    write(workspace, keyed_ledger(None, notes="toy_demo: not a divergence.", strict=True), "DIFFER")
    assert main(["ledger-audit", "--results", "results"]) == 1


def test_strict_passes_when_every_group_is_listed(workspace: Path):
    write_runs(workspace, *DIFFERENT)
    write(workspace, keyed_ledger("[written.lines]", strict=True), "DIFFER")
    assert main(["ledger-audit", "--results", "results"]) == 0


def test_strict_reports_a_difference_group_that_no_row_lists(workspace: Path, capsys):
    write_runs(workspace, {"written": {"lines": [1, 2]}, "other": {"x": 1}}, {"written": {"lines": [1]}, "other": {"x": 2}})
    write(workspace, keyed_ledger("[written.lines]", strict=True), "DIFFER")
    assert main(["ledger-audit", "--results", "results"]) == 1
    assert "difference group other.x is listed in the differs of no row" in capsys.readouterr().out


def test_strict_lets_another_row_cover_a_group(workspace: Path):
    write_runs(workspace, {"written": {"lines": [1, 2]}, "other": {"x": 1}}, {"written": {"lines": [1]}, "other": {"x": 2}})
    second = (
        "  - matlab: other\n    python: demopy.other\n    status: parity\n    scenarios:\n      - toy_demo\n"
        "    differs: [other]\n"
    )
    write(workspace, keyed_ledger("[written]", strict=True, extra_row=second), "DIFFER")
    assert main(["ledger-audit", "--results", "results"]) == 0


def test_strict_empty_list_needs_a_note(workspace: Path, capsys):
    write_runs(workspace, *DIFFERENT)
    second = "  - matlab: other\n    python: demopy.other\n    status: parity\n    scenarios:\n      - toy_demo\n    differs: []\n"
    write(workspace, keyed_ledger("[written.lines]", strict=True, extra_row=second), "DIFFER")
    assert main(["ledger-audit", "--results", "results"]) == 1
    assert "empty differs list needs a note" in capsys.readouterr().out
    second = second + "    notes: 'No difference of its own: the scenario differs elsewhere.'\n"
    write(workspace, keyed_ledger("[written.lines]", strict=True, extra_row=second), "DIFFER")
    assert main(["ledger-audit", "--results", "results"]) == 0


def test_strict_flags_a_difference_scenario_without_result_files(workspace: Path, capsys):
    write(workspace, keyed_ledger("[written]", strict=True), "DIFFER")
    assert main(["ledger-audit", "--results", "results"]) == 1
    assert "result files are missing" in capsys.readouterr().out


def test_the_keys_survive_a_round_trip_of_the_ledger(workspace: Path):
    from parity.ledger import dump_ledger, load_ledger

    write(workspace, keyed_ledger("[written, kept]", strict=True), "DIFFER")
    ledger_path = workspace / "ledgers" / "demo.yml"
    loaded = load_ledger(ledger_path)
    assert loaded.strict_differs
    assert loaded.entries[0].differs == ["written", "kept"]
    ledger_path.write_text(dump_ledger(loaded), encoding="utf-8")
    again = load_ledger(ledger_path)
    assert again.strict_differs and again.entries[0].differs == ["written", "kept"]
