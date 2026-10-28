"""Ledger audit: a row declared at ``parity`` must say what is known to differ.

Two rules, both about rows with status ``parity``:

* a scenario recorded as ``DIFFER`` or ``ERROR`` needs a note that mentions the
  difference (the word "differ" or "diverge" in any form) or names the scenario, which is how a
  row says that a scenario covering several functions differs elsewhere;
* a scenario whose ``entries`` name the row's function must be listed in the row's
  ``scenarios``, so the row and the evidence stay together.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from parity.ledger import Ledger
from parity.scenarios import Scenario

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


def audit_ledger(
    ledger: Ledger, scenarios: list[Scenario], results: dict[str, str] | None = None
) -> list[AuditFinding]:
    """Findings for the ``parity`` rows of ``ledger``; ``results`` maps scenario id to its last verdict."""
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
    return findings
