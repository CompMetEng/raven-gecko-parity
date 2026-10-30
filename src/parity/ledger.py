"""The parity ledger: one declared status per public function, on both sides."""

from __future__ import annotations

import difflib
import io
import re
from dataclasses import dataclass, field, replace
from pathlib import Path

import yaml

#: Keys the top level of a ledger file reads. Other keys are left alone (a ledger may hold a
#: YAML anchor there for merge keys), unless one looks like a misspelling of these.
TOP_LEVEL_KEYS = frozenset({"version", "pair", "matlab", "python", "strict_differs", "entries"})


def path_segments(path: str) -> tuple[str, ...]:
    """Dotted checkpoint path without list indices and without the leading ``results``."""
    parts = [part for part in re.sub(r"\[\d+\]", "", str(path)).split(".") if part]
    if parts and parts[0] == "results":
        parts = parts[1:]
    return tuple(parts)


#: Every public function on either side must carry exactly one of these.
STATUSES: dict[str, str] = {
    "parity": "Both sides implement this and are intended to behave the same.",
    "python-pending": "MATLAB has it; a Python port is queued.",
    "matlab-pending": "Python has it; a MATLAB back-port is queued.",
    "matlab-only": "Deliberately MATLAB-only --- not coming to Python.",
    "python-only": "Deliberately Python-only --- not going to MATLAB.",
    "via-dependency": "The other side gets this from a dependency (cobrapy, COBRA Toolbox) instead.",
    "subsumed": "The other side has the capability inside another function, not as its own entry point.",
    "internal": "Not part of the cross-implementation API (glue, installers, helpers).",
    "unreviewed": "Not yet triaged. Drive this count to zero.",
}

#: Statuses that must justify themselves in prose.
REASON_REQUIRED = frozenset({"matlab-only", "python-only", "via-dependency", "subsumed", "internal"})

#: Statuses that describe a settled one-sided decision rather than queued work.
ONE_SIDED = frozenset({"matlab-only", "python-only", "via-dependency", "subsumed"})

#: Statuses that represent queued work rather than a settled decision.
PENDING = frozenset({"python-pending", "matlab-pending"})

FIELD_ORDER = ("matlab", "python", "status", "reason", "issue", "scenarios", "differs", "notes")


class LedgerError(RuntimeError):
    """Raised when a ledger file cannot be parsed into entries."""


@dataclass
class Entry:
    """One row: a function on one or both sides, plus why it is where it is."""

    status: str
    matlab: str | None = None
    python: str | None = None
    reason: str | None = None
    issue: str | None = None
    notes: str | None = None
    scenarios: list[str] = field(default_factory=list)
    #: Checkpoint path prefixes (without the leading ``results.``) where the linked scenarios
    #: still differ. ``None`` means the row does not say; ``[]`` says the row has no
    #: difference of its own.
    differs: list[str] | None = None

    @property
    def label(self) -> str:
        """A stable human-facing name for the row."""
        return self.matlab or self.python or "<empty>"

    @property
    def sort_key(self) -> tuple[str, str]:
        return ((self.matlab or "~").lower(), (self.python or "~").lower())

    def to_dict(self) -> dict:
        out: dict[str, object] = {}
        for name in FIELD_ORDER:
            value = getattr(self, name)
            if name == "differs":
                if value is not None:
                    out[name] = list(value)
                continue
            if value in (None, "", []):
                continue
            out[name] = value
        return out


@dataclass
class Ledger:
    pair: str
    matlab_repo: str
    python_repo: str
    entries: list[Entry]
    version: int = 1
    path: Path | None = None
    #: When true, ``parity ledger-audit`` requires every parity row linked to a scenario that
    #: differs to list where in its ``differs`` field, and every difference to be listed.
    strict_differs: bool = False

    def by_matlab(self) -> dict[str, Entry]:
        return {e.matlab: e for e in self.entries if e.matlab}

    def by_python(self) -> dict[str, Entry]:
        return {e.python: e for e in self.entries if e.python}

    def with_status(self, *statuses: str) -> list[Entry]:
        wanted = set(statuses)
        return [e for e in self.entries if e.status in wanted]

    def counts(self) -> dict[str, int]:
        out = {status: 0 for status in STATUSES}
        for entry in self.entries:
            out[entry.status] = out.get(entry.status, 0) + 1
        return out


class _UniqueKeyLoader(yaml.SafeLoader):
    """A safe loader that refuses a mapping with a repeated key.

    ``yaml.safe_load`` keeps the last value of a repeated key, so a deleted ``- matlab:`` line
    would merge two entries into one without any message. Merge keys (``<<``) are expanded as
    ``yaml.safe_load`` expands them, and a key written next to a merge overrides the merged value
    without being a repetition; only keys written in the mapping itself are compared. A key that
    cannot be hashed is left to PyYAML, which raises its own error for it.
    """

    def construct_mapping(self, node, deep=False):
        seen: set = set()
        for key_node, _ in node.value:
            if key_node.tag == "tag:yaml.org,2002:merge":
                continue
            key = self.construct_object(key_node, deep=True)
            try:
                if key in seen:
                    raise LedgerError(
                        f"line {key_node.start_mark.line + 1}: key {key!r} appears twice in one mapping"
                    )
                seen.add(key)
            except TypeError:
                break
        return super().construct_mapping(node, deep)


def load_ledger(path: Path) -> Ledger:
    """Read a ledger YAML file."""
    if not path.is_file():
        raise LedgerError(f"ledger not found: {path}")

    try:
        data = yaml.load(path.read_text(encoding="utf-8"), Loader=_UniqueKeyLoader) or {}
    except LedgerError as exc:
        raise LedgerError(f"{path}: {exc}") from None
    if not isinstance(data, dict):
        raise LedgerError(f"{path}: expected a mapping at the top level")

    for key in data:
        if key in TOP_LEVEL_KEYS:
            continue
        close = difflib.get_close_matches(str(key), sorted(TOP_LEVEL_KEYS), n=1, cutoff=0.8)
        if close:
            raise LedgerError(f"{path}: unknown top-level key {key!r}; did you mean {close[0]!r}?")
    strict = data.get("strict_differs", False)
    if not isinstance(strict, bool):
        raise LedgerError(f"{path}: 'strict_differs' must be true or false, not {strict!r}")

    raw_entries = data.get("entries") or []
    if not isinstance(raw_entries, list):
        raise LedgerError(f"{path}: 'entries' must be a list")

    entries: list[Entry] = []
    for index, raw in enumerate(raw_entries, start=1):
        if not isinstance(raw, dict):
            raise LedgerError(f"{path}: entry {index} is not a mapping")
        unknown = set(raw) - set(FIELD_ORDER)
        if unknown:
            raise LedgerError(f"{path}: entry {index} has unknown field(s): {', '.join(sorted(unknown))}")
        if "status" not in raw:
            raise LedgerError(f"{path}: entry {index} has no 'status'")
        scenarios = raw.get("scenarios") or []
        if isinstance(scenarios, str):
            scenarios = [scenarios]
        differs = raw.get("differs")
        if isinstance(differs, str):
            differs = [differs]
        if differs is not None and not isinstance(differs, list):
            raise LedgerError(f"{path}: entry {index} has a 'differs' that is not a list")
        for item in differs or []:
            if not path_segments(item):
                raise LedgerError(
                    f"{path}: entry {index} has a differs entry {item!r} that names no checkpoint path"
                )
        entries.append(
            Entry(
                status=str(raw["status"]),
                matlab=raw.get("matlab"),
                python=raw.get("python"),
                reason=raw.get("reason"),
                issue=raw.get("issue"),
                notes=raw.get("notes"),
                scenarios=list(scenarios),
                differs=None if differs is None else [str(item) for item in differs],
            )
        )

    return Ledger(
        pair=data.get("pair", path.stem),
        matlab_repo=data.get("matlab", ""),
        python_repo=data.get("python", ""),
        entries=entries,
        version=int(data.get("version", 1)),
        path=path,
        strict_differs=strict,
    )


HEADER = """\
# Parity ledger --- {pair}
#
# One row per public function on either side. `parity check` fails if a public
# function is missing here, or if a row names something that no longer exists.
#
# status:
{status_help}
#
# Regenerate the human-readable view with: parity report --pair {pair}
"""


def dump_ledger(ledger: Ledger) -> str:
    """Serialise a ledger deterministically, so diffs stay readable."""
    status_help = "\n".join(f"#   {name:<15} {help_}" for name, help_ in STATUSES.items())
    buffer = io.StringIO()
    buffer.write(HEADER.format(pair=ledger.pair, status_help=status_help))
    buffer.write("\n")

    head = {
        "version": ledger.version,
        "pair": ledger.pair,
        "matlab": ledger.matlab_repo,
        "python": ledger.python_repo,
    }
    if ledger.strict_differs:
        head["strict_differs"] = True
    buffer.write(yaml.safe_dump(head, sort_keys=False, allow_unicode=True))
    buffer.write("\nentries:\n")

    for entry in sorted(ledger.entries, key=lambda e: e.sort_key):
        rendered = yaml.safe_dump(
            entry.to_dict(), sort_keys=False, allow_unicode=True, default_flow_style=False, width=100
        )
        lines = rendered.rstrip("\n").split("\n")
        buffer.write(f"  - {lines[0]}\n")
        for line in lines[1:]:
            buffer.write(f"    {line}\n")

    return buffer.getvalue()


def save_ledger(ledger: Ledger, path: Path | None = None) -> Path:
    target = path or ledger.path
    if target is None:
        raise LedgerError("no path to save the ledger to")
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(dump_ledger(ledger), encoding="utf-8")
    return target


def merge_entries(existing: list[Entry], discovered: list[Entry]) -> list[Entry]:
    """Add *discovered* rows for anything *existing* does not already cover.

    Hand-written decisions always win: a discovered row is dropped if either of its
    sides is already claimed by an existing row.
    """
    claimed_matlab = {e.matlab for e in existing if e.matlab}
    claimed_python = {e.python for e in existing if e.python}

    merged = list(existing)
    for entry in discovered:
        if entry.matlab and entry.matlab in claimed_matlab:
            continue
        if entry.python and entry.python in claimed_python:
            continue
        merged.append(replace(entry))
        if entry.matlab:
            claimed_matlab.add(entry.matlab)
        if entry.python:
            claimed_python.add(entry.python)
    return merged
