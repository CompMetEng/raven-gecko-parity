#!/usr/bin/env python3
"""Checks for the decision records of a repository.

A record is a Markdown file `NNNN-title.md` in the records folder (default `decisions`)
that starts with a frontmatter block:

    ---
    status: accepted
    superseded_by: 0007-new-title.md
    ---

Statuses and what they allow:

    proposed    editable; may become accepted or rejected
    accepted    frozen; may become superseded
    rejected    frozen; may become superseded
    superseded  frozen; final; names its successor in `superseded_by`

The changes to a frozen record are limited to its `status` and `superseded_by` lines and
one appended note line that starts with "> ". It is neither deleted nor renamed.
`README.md` in the folder is not a record; any other file name is a finding.

Usage:
    decision_records.py check [--dir DIR]
        the records in the working tree: names, unique numbers, statuses, successors, links
    decision_records.py diff --base REV --head REV [--dir DIR]
        what the commits between REV and REV change in the records
    decision_records.py --version

A folder that does not exist is not an error. Every finding is printed as one line, and the
exit status is 1 when there is one.

shortcut: a push is checked after it landed, change when the integration branch refuses
direct pushes.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path, PurePosixPath

VERSION = "1"
DEFAULT_DIR = "decisions"
STATUSES = ("proposed", "accepted", "rejected", "superseded")
FROZEN = frozenset({"accepted", "rejected", "superseded"})
TRANSITIONS = {
    "proposed": frozenset({"accepted", "rejected"}),
    "accepted": frozenset({"superseded"}),
    "rejected": frozenset({"superseded"}),
    "superseded": frozenset(),
}
RECORD_NAME = re.compile(r"^(\d{4})-[a-z0-9][a-z0-9-]*\.md$")
IGNORED_NAMES = frozenset({"README.md"})
KEY_LINE = re.compile(r"^([A-Za-z_][A-Za-z0-9_-]*):[ \t]*(.*?)[ \t]*$")
STATE_KEYS = frozenset({"status", "superseded_by"})
LINK = re.compile(r"\]\(([^)\s]+)\)")
URL_SCHEME = re.compile(r"^[A-Za-z][A-Za-z0-9+.-]*:")
NOTE = re.compile(r"^> \S")
FENCE = re.compile(r"^\s*(```|~~~)")


def normalize(text: str) -> str:
    return text.replace("\r\n", "\n").replace("\r", "\n")


def frontmatter(text: str) -> dict[str, str] | None:
    """The `key: value` lines of the leading `---` block, or None when there is none."""
    lines = text.split("\n")
    if lines[0].rstrip() != "---":
        return None
    fields: dict[str, str] = {}
    for line in lines[1:]:
        if line.rstrip() == "---":
            return fields
        match = KEY_LINE.match(line)
        if match:
            fields[match.group(1)] = match.group(2).strip("\"'")
    return None


def without_state_lines(text: str) -> list[str]:
    """The lines of a record without its `status` and `superseded_by` lines and trailing blank lines."""
    lines = text.split("\n")
    in_front = lines[0].rstrip() == "---"
    kept = []
    for index, line in enumerate(lines):
        if in_front and index > 0:
            if line.rstrip() == "---":
                in_front = False
            else:
                match = KEY_LINE.match(line)
                if match and match.group(1) in STATE_KEYS:
                    continue
        kept.append(line)
    while kept and not kept[-1].strip():
        kept.pop()
    return kept


def frozen_change_is_allowed(old: str, new: str) -> bool:
    """True when `new` differs from `old` only in the state lines and one appended note."""
    before, after = without_state_lines(old), without_state_lines(new)
    if after[: len(before)] != before:
        return False
    added = [line for line in after[len(before):] if line.strip()]
    return not added or (len(added) == 1 and NOTE.match(added[0]) is not None)


def links_of(text: str) -> list[str]:
    """Relative link targets outside fenced code blocks, without fragments."""
    targets = []
    fenced = False
    for line in text.split("\n"):
        if FENCE.match(line):
            fenced = not fenced
            continue
        if fenced:
            continue
        for match in LINK.finditer(line):
            target = match.group(1).split("#", 1)[0]
            if target and not URL_SCHEME.match(target):
                targets.append(target)
    return targets


def check_tree(directory: str) -> list[str]:
    root = Path(directory)
    if not root.is_dir():
        return []
    findings: list[str] = []
    names_by_number: dict[str, list[str]] = {}
    status_of: dict[str, str | None] = {}
    successor_of: dict[str, str] = {}
    texts: dict[str, str] = {}

    for path in sorted(root.iterdir()):
        if path.is_dir() or path.name in IGNORED_NAMES:
            continue
        label = f"{directory}/{path.name}"
        match = RECORD_NAME.match(path.name)
        if not match:
            findings.append(f"{label}: not a record name (NNNN-title.md with lowercase letters, digits and hyphens)")
            continue
        names_by_number.setdefault(match.group(1), []).append(path.name)
        text = normalize(path.read_text(encoding="utf-8"))
        texts[path.name] = text
        fields = frontmatter(text)
        if fields is None:
            findings.append(f"{label}: no frontmatter block")
            continue
        status = fields.get("status")
        status_of[path.name] = status
        if status not in STATUSES:
            findings.append(f"{label}: status '{status}' is not one of {', '.join(STATUSES)}")
        if status == "superseded":
            if not fields.get("superseded_by"):
                findings.append(f"{label}: a superseded record needs superseded_by")
            else:
                successor_of[path.name] = fields["superseded_by"]
        elif "superseded_by" in fields:
            findings.append(f"{label}: superseded_by on a record that is not superseded")

    for number, names in sorted(names_by_number.items()):
        if len(names) > 1:
            findings.append(f"{directory}: number {number} is used by {', '.join(names)}")

    for name, target in sorted(successor_of.items()):
        label = f"{directory}/{name}"
        target_match = RECORD_NAME.match(target)
        if not target_match or not (root / target).is_file():
            findings.append(f"{label}: superseded_by '{target}' is not a record in {directory}")
            continue
        if target_match.group(1) <= name[:4]:
            findings.append(f"{label}: the successor {target} does not have a higher number")
        if status_of.get(target) != "accepted":
            findings.append(f"{label}: the successor {target} is not accepted")

    for name, text in sorted(texts.items()):
        for target in links_of(text):
            if not (root / target).exists():
                findings.append(f"{directory}/{name}: link '{target}' does not resolve")
    return findings


def git_lines(*args: str) -> list[str]:
    out = subprocess.run(["git", *args], capture_output=True, check=True).stdout
    return [item.decode("utf-8", errors="replace") for item in out.split(b"\0") if item]


def resolves(rev: str) -> bool:
    probe = subprocess.run(["git", "rev-parse", "--verify", "--quiet", f"{rev}^{{commit}}"], capture_output=True)
    return probe.returncode == 0


def show(rev: str, path: str) -> str | None:
    result = subprocess.run(["git", "show", f"{rev}:{path}"], capture_output=True)
    if result.returncode:
        return None
    return normalize(result.stdout.decode("utf-8", errors="replace"))


def records_at(rev: str, directory: str) -> dict[str, str]:
    """Record number -> path of the records directly in `directory` at `rev`."""
    found = {}
    for path in git_lines("ls-tree", "-r", "-z", "--name-only", rev, "--", directory):
        pure = PurePosixPath(path)
        match = RECORD_NAME.match(pure.name)
        if match and pure.parent == PurePosixPath(directory):
            found.setdefault(match.group(1), path)
    return found


def status_in(text: str | None) -> str | None:
    fields = frontmatter(text) if text is not None else None
    return fields.get("status") if fields else None


def check_diff(base: str, head: str, directory: str) -> list[str]:
    if not resolves(base):
        print(f"decision-records: {base} is not in this clone, nothing to compare", file=sys.stderr)
        return []
    merge_base = subprocess.run(["git", "merge-base", base, head], capture_output=True, text=True)
    if merge_base.returncode:
        return [f"{directory}: {base} and {head} have no merge base"]
    fork = merge_base.stdout.strip()

    at_fork = records_at(fork, directory)
    at_base = records_at(base, directory)
    at_head = records_at(head, directory)
    findings: list[str] = []

    for number in sorted(set(at_fork) | set(at_head)):
        old_path, new_path = at_fork.get(number), at_head.get(number)
        old = show(fork, old_path) if old_path else None
        new = show(head, new_path) if new_path else None
        if old_path == new_path and old == new:
            continue
        if old is None and number in at_base:
            continue
        base_text = show(base, at_base[number]) if number in at_base else None
        statuses = {s for s in (status_in(old), status_in(base_text)) if s}
        frozen = sorted(statuses & FROZEN)
        label = new_path or old_path or number

        if new is None:
            if frozen:
                findings.append(f"{label}: the {frozen[0]} record is deleted")
            continue
        new_status = status_in(new)
        if old is None:
            if new_status == "superseded":
                findings.append(f"{label}: a new record cannot start as superseded")
            continue
        if frozen and old_path != new_path:
            findings.append(f"{label}: the {frozen[0]} record is renamed from {old_path}")
        elif frozen and not frozen_change_is_allowed(old, new):
            findings.append(
                f"{label}: the {frozen[0]} record changes beyond its status, superseded_by and one appended note"
            )
        if new_status in STATUSES and new_status != status_in(old):
            for status in sorted(statuses):
                if new_status != status and new_status not in TRANSITIONS[status]:
                    findings.append(f"{label}: status {status} -> {new_status} is not allowed")
    return findings


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(prog="decision_records.py", description="Checks for decision records.")
    parser.add_argument("--version", action="version", version=VERSION)
    commands = parser.add_subparsers(dest="command", required=True)
    check = commands.add_parser("check", help="check the records in the working tree")
    check.add_argument("--dir", default=DEFAULT_DIR)
    diff = commands.add_parser("diff", help="check what the commits between two revisions change")
    diff.add_argument("--base", required=True)
    diff.add_argument("--head", required=True)
    diff.add_argument("--dir", default=DEFAULT_DIR)
    args = parser.parse_args(argv[1:])

    directory = args.dir.strip("/")
    if args.command == "check":
        found = check_tree(directory)
    else:
        found = check_diff(args.base, args.head, directory)
    for line in found:
        print(f"decision-records: {line}")
    return 1 if found else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
