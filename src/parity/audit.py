"""Ledger audit: a row declared at ``parity`` must say what is known to differ.

Rules about rows with status ``parity``:

* a scenario recorded as ``DIFFER`` or ``ERROR`` needs a note that mentions the
  difference (the word "differ" or "diverge" in any form) or names the scenario, which is how a
  row says that a scenario covering several functions differs elsewhere;
* a scenario whose ``entries`` name the row's function must be listed in the row's
  ``scenarios``, so the row and the evidence stay together.

Rows may also carry ``differs``, a list of checkpoint path prefixes (written without the
leading ``results.`` and without list indices, for example ``written.lines``). The audit
rebuilds the differences of each linked scenario from the result files of the run
(``python/<id>.seed<N>.json`` against ``matlab/<id>.json``) and checks the list against them:

* every listed prefix has to be the start of the path of a difference of one of the row's
  scenarios, so a list cannot go stale; a scenario without a verdict in the results folder (a
  partial run) is not checked;
* a row with ``differs`` is not held to the note rule above, and ``differs: []`` needs a note
  that says the row has no difference of its own.

A ledger that sets ``strict_differs: true`` is audited more strictly. Every parity row linked
to a scenario that differs must have a ``differs`` list, and every group of differences of such
a scenario (its first two path segments) has to be covered by the list of at least one row that
links the scenario: an entry covers a group when it is the group or goes deeper. A difference at
the root of the results (one side returned no results) belongs to no group and is always
reported. Rows and ledgers without these keys are audited by the first two rules only.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path

from parity.ledger import Ledger, path_segments
from parity.scenarios import Scenario, compare_results

#: Words that count as the note admitting a difference.
_ADMITS = re.compile(r"differ|diverg", re.IGNORECASE)

#: Recorded results that mean the two sides disagree or could not be compared.
BAD_RESULTS = frozenset({"DIFFER", "ERROR"})


@dataclass(frozen=True)
class AuditFinding:
    row: str
    message: str

    def format(self) -> str:
        return f"{self.row}: {self.message}"


def read_results(directory: Path) -> dict[str, str]:
    """Last recorded result per scenario id from ``<directory>/status/<id>.txt`` or ``<directory>/<id>.txt``.

    The file holds the verdict as its first word (``MATCH``, ``DIFFER`` or ``ERROR``).
    """
    folder = directory / "status" if (directory / "status").is_dir() else directory
    results: dict[str, str] = {}
    for path in sorted(folder.glob("*.txt")):
        words = path.read_text(encoding="utf-8").split()
        if words:
            results[path.stem] = words[0].upper()
    return results


def _segments(path: str) -> tuple[str, ...]:
    return path_segments(path)


def _group(path: str) -> tuple[str, ...]:
    return _segments(path)[:2]


def _is_prefix(entry: tuple[str, ...], path: tuple[str, ...]) -> bool:
    """``entry`` is the start of ``path``, at a segment boundary."""
    return bool(entry) and path[: len(entry)] == entry


def _covers(entry: tuple[str, ...], group: tuple[str, ...]) -> bool:
    """``entry`` is the group or goes deeper, so it names no more than the group holds."""
    return bool(group) and len(entry) >= len(group) and entry[: len(group)] == group


def _result_root(directory: Path) -> Path:
    """Folder that holds ``python/`` and ``matlab/``: the given one, or its parent when the
    ``status/`` folder itself was given."""
    for candidate in (directory, directory.parent):
        if (candidate / "matlab").is_dir() or (candidate / "python").is_dir():
            return candidate
    return directory


def read_differences(
    directory: Path, scenarios: list[Scenario], errors: dict[str, str] | None = None
) -> dict[str, list[str]]:
    """Differences per scenario id, rebuilt from the result files of a run.

    Python (every seed) is compared with MATLAB within the scenario's tolerance. A scenario
    without both result files is left out. A result file that cannot be read leaves the scenario
    out too and, when ``errors`` is given, records the reason there.
    """
    root = _result_root(directory)
    out: dict[str, list[str]] = {}
    for scenario in scenarios:
        matlab_file = root / "matlab" / f"{scenario.id}.json"
        python_files = sorted((root / "python").glob(f"{scenario.id}.seed*.json"))
        if not matlab_file.is_file() or not python_files:
            continue
        paths: list[str] = []
        current = matlab_file
        try:
            matlab = json.loads(matlab_file.read_text(encoding="utf-8"))
            for python_file in python_files:
                current = python_file
                python = json.loads(python_file.read_text(encoding="utf-8"))
                for difference in compare_results(python, matlab, scenario.tolerance):
                    if difference.path not in paths:
                        paths.append(difference.path)
        except (ValueError, OSError) as error:  # JSONDecodeError and UnicodeDecodeError are ValueErrors
            if errors is not None:
                errors[scenario.id] = f"cannot read {current.name}: {error}"
            continue
        out[scenario.id] = paths
    return out


def _audit_differs(
    ledger: Ledger,
    scenarios: list[Scenario],
    results: dict[str, str],
    differences: dict[str, list[str]] | None,
    read_errors: dict[str, str] | None = None,
) -> list[AuditFinding]:
    findings: list[AuditFinding] = []
    bad = {s.id for s in scenarios if results.get(s.id) in BAD_RESULTS}
    rows_for: dict[str, list] = {}
    for entry in ledger.entries:
        for scenario_id in entry.scenarios:
            rows_for.setdefault(scenario_id, []).append(entry)

    for scenario_id, message in sorted((read_errors or {}).items()):
        findings.append(AuditFinding(scenario_id, f"its differences cannot be checked: {message}"))

    for entry in ledger.entries:
        if entry.differs is not None and not entry.differs and not (entry.notes or "").strip():
            findings.append(AuditFinding(entry.label, "an empty differs list needs a note that says the row has no difference of its own"))

    if differences is not None:
        for entry in ledger.entries:
            # Only scenarios with a verdict and rebuilt differences in this run can show a stale entry.
            checkable = [s for s in entry.scenarios if s in results and s in differences]
            if not checkable:
                continue
            for item in entry.differs or []:
                wanted = _segments(item)
                seen = any(
                    _is_prefix(wanted, _segments(path))
                    for scenario_id in checkable
                    if scenario_id in bad
                    for path in differences[scenario_id]
                )
                if not seen:
                    findings.append(
                        AuditFinding(entry.label, f"differs entry {item!r} is the start of no difference of the row's scenarios in this run")
                    )

    if not ledger.strict_differs:
        return findings

    for entry in ledger.entries:
        if entry.status == "parity" and entry.differs is None:
            linked = [s for s in entry.scenarios if s in bad]
            if linked:
                findings.append(
                    AuditFinding(entry.label, f"scenario {linked[0]} is recorded as {results[linked[0]]} but the row has no differs list")
                )

    if differences is not None:
        for scenario in scenarios:
            if scenario.id not in bad:
                continue
            if scenario.id not in differences:
                if scenario.id not in (read_errors or {}):
                    findings.append(
                        AuditFinding(scenario.id, f"recorded as {results[scenario.id]} but its result files are missing, so its differences cannot be checked")
                    )
                continue
            listed = [_segments(item) for entry in rows_for.get(scenario.id, []) for item in entry.differs or []]
            groups = sorted({_group(path) for path in differences[scenario.id]})
            for group in groups:
                if not group:
                    findings.append(
                        AuditFinding(scenario.id, "a difference at the root of the results (one side returned no results) belongs to no group and cannot be listed")
                    )
                elif not any(_covers(item, group) for item in listed):
                    findings.append(
                        AuditFinding(scenario.id, f"difference group {'.'.join(group)} is listed in the differs of no row (an entry has to be the group or go deeper)")
                    )
    return findings


def audit_ledger(
    ledger: Ledger,
    scenarios: list[Scenario],
    results: dict[str, str] | None = None,
    differences: dict[str, list[str]] | None = None,
    read_errors: dict[str, str] | None = None,
) -> list[AuditFinding]:
    """Findings for the ``parity`` rows of ``ledger``.

    ``results`` maps scenario id to its last verdict; ``differences`` maps scenario id to the
    difference paths rebuilt from the result files and ``read_errors`` to the reason a result
    file could not be read (see :func:`read_differences`).
    """
    findings: list[AuditFinding] = []
    by_function: dict[str, list[str]] = {}
    for scenario in scenarios:
        for name in scenario.entries:
            by_function.setdefault(name, []).append(scenario.id)

    for entry in ledger.entries:
        if entry.status != "parity":
            continue
        for scenario_id in entry.scenarios:
            verdict = (results or {}).get(scenario_id)
            notes = entry.notes or ""
            if ledger.strict_differs or entry.differs is not None:
                continue
            if verdict in BAD_RESULTS and not (_ADMITS.search(notes) or scenario_id in notes):
                findings.append(
                    AuditFinding(entry.label, f"scenario {scenario_id} is recorded as {verdict} but the note neither names a difference nor the scenario")
                )
        for name in (entry.matlab, entry.python):
            for scenario_id in by_function.get(name or "", []):
                if scenario_id not in entry.scenarios:
                    findings.append(
                        AuditFinding(entry.label, f"scenario {scenario_id} covers {name} but is not listed in the row's scenarios")
                    )
    if results is not None:
        findings.extend(_audit_differs(ledger, scenarios, results, differences, read_errors))
    return findings
