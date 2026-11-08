"""The decision-record lint that `parity docs render` copies into the code repositories.

The script is a stand-alone file under `docs/shared/scripts/`, not part of the `parity`
package, so it is loaded by path. The diff checks run against real git repositories built in
a temporary directory, because what they decide depends on merge bases, renames and line
endings.
"""

from __future__ import annotations

import importlib.util
import subprocess
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "docs" / "shared" / "scripts" / "decision_records.py"
spec = importlib.util.spec_from_file_location("decision_records", SCRIPT)
dr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(dr)


def record(status: str = "proposed", body: str = "The decision.\n", extra: str = "") -> str:
    return f"---\nstatus: {status}\n{extra}---\n\n# Title\n\n{body}"


def git(repo: Path, *args: str) -> str:
    result = subprocess.run(
        ["git", "-c", "user.name=T", "-c", "user.email=t@example.org", "-c", "commit.gpgsign=false", *args],
        cwd=repo,
        capture_output=True,
        text=True,
        check=True,
    )
    return result.stdout


def write(repo: Path, relative: str, text: str) -> None:
    path = repo / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(text.encode("utf-8"))


def commit(repo: Path, message: str = "change") -> str:
    git(repo, "add", "-A")
    git(repo, "commit", "-q", "-m", message)
    return git(repo, "rev-parse", "HEAD").strip()


@pytest.fixture
def repo(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    git(tmp_path, "init", "-q", "-b", "main")
    git(tmp_path, "config", "core.autocrlf", "false")
    monkeypatch.chdir(tmp_path)
    return tmp_path


def diff(base: str, head: str) -> list[str]:
    return dr.check_diff(base, head, "decisions")


# --- check ---------------------------------------------------------------------------------


def test_a_missing_folder_is_not_an_error(repo: Path):
    assert dr.check_tree("decisions") == []


def test_valid_records_pass(repo: Path):
    write(repo, "decisions/README.md", "# Decisions\n")
    write(repo, "decisions/0001-first.md", record("superseded", extra="superseded_by: 0002-second.md\n"))
    write(repo, "decisions/0002-second.md", record("accepted", body="See [the first](0001-first.md#top).\n"))
    write(repo, "decisions/0003-third.md", record("proposed"))
    assert dr.check_tree("decisions") == []


@pytest.mark.parametrize(
    ("name", "text", "expected"),
    [
        ("0001-a.md", "# No frontmatter\n", "no frontmatter block"),
        ("0001-a.md", "---\ntitle: x\n---\n", "status 'None'"),
        ("0001-a.md", record("draft"), "status 'draft'"),
        ("0001-a.md", record("superseded"), "needs superseded_by"),
        ("0001-a.md", record("accepted", extra="superseded_by: 0002-b.md\n"), "not superseded"),
        ("0001-a.md", record(body="[gone](missing.md)\n"), "link 'missing.md' does not resolve"),
        ("notes.md", record(), "not a record name"),
        ("001-short.md", record(), "not a record name"),
        ("0001-Upper.md", record(), "not a record name"),
    ],
)
def test_check_reports_a_bad_record(repo: Path, name: str, text: str, expected: str):
    write(repo, f"decisions/{name}", text)
    assert any(expected in line for line in dr.check_tree("decisions"))


def test_links_in_code_fences_and_urls_are_not_checked(repo: Path):
    body = "[web](https://example.org/x.md) [mail](mailto:a@b.c) [anchor](#top)\n\n```\n[x](nowhere.md)\n```\n"
    write(repo, "decisions/0001-a.md", record(body=body))
    assert dr.check_tree("decisions") == []


def test_a_successor_must_exist_be_newer_and_be_accepted(repo: Path):
    write(repo, "decisions/0002-old.md", record("superseded", extra="superseded_by: 0001-new.md\n"))
    write(repo, "decisions/0001-new.md", record("proposed"))
    write(repo, "decisions/0003-lost.md", record("superseded", extra="superseded_by: 0009-none.md\n"))
    found = dr.check_tree("decisions")
    assert any("0002-old.md: the successor 0001-new.md does not have a higher number" in f for f in found)
    assert any("0002-old.md: the successor 0001-new.md is not accepted" in f for f in found)
    assert any("0003-lost.md: superseded_by '0009-none.md' is not a record" in f for f in found)


def test_a_number_is_used_once(repo: Path):
    write(repo, "decisions/0003-alpha.md", record())
    write(repo, "decisions/0003-beta.md", record())
    assert any("number 0003 is used by 0003-alpha.md, 0003-beta.md" in f for f in dr.check_tree("decisions"))


def test_two_branches_that_add_the_same_number_are_found_after_the_merge(repo: Path):
    write(repo, "decisions/0001-first.md", record("accepted"))
    commit(repo)
    git(repo, "switch", "-q", "-c", "a")
    write(repo, "decisions/0003-alpha.md", record())
    commit(repo)
    git(repo, "switch", "-q", "main")
    git(repo, "switch", "-q", "-c", "b")
    write(repo, "decisions/0003-beta.md", record())
    commit(repo)
    git(repo, "switch", "-q", "a")
    git(repo, "merge", "-q", "--no-edit", "b")
    assert any("number 0003 is used" in f for f in dr.check_tree("decisions"))


# --- diff ----------------------------------------------------------------------------------


def branch_from_main(repo: Path, initial: dict[str, str]) -> str:
    """Commit `initial` on main and return the SHA."""
    for relative, text in initial.items():
        write(repo, relative, text)
    return commit(repo, "base")


def on_branch(repo: Path, name: str = "pr"):
    git(repo, "switch", "-q", "-c", name)


def test_body_edit_of_an_accepted_record_is_found(repo: Path):
    base = branch_from_main(repo, {"decisions/0001-a.md": record("accepted")})
    on_branch(repo)
    write(repo, "decisions/0001-a.md", record("accepted", body="Another decision.\n"))
    head = commit(repo)
    found = diff(base, head)
    assert found == ["decisions/0001-a.md: the accepted record changes beyond its status, superseded_by and one appended note"]


def test_body_edit_of_a_proposed_record_is_allowed(repo: Path):
    base = branch_from_main(repo, {"decisions/0001-a.md": record("proposed")})
    on_branch(repo)
    write(repo, "decisions/0001-a.md", record("proposed", body="Edited.\n"))
    assert diff(base, commit(repo)) == []


def test_superseding_with_a_note_is_allowed(repo: Path):
    base = branch_from_main(
        repo, {"decisions/0001-a.md": record("accepted"), "decisions/0002-b.md": record("proposed")}
    )
    on_branch(repo)
    write(repo, "decisions/0002-b.md", record("accepted"))
    write(
        repo,
        "decisions/0001-a.md",
        record("superseded", extra="superseded_by: 0002-b.md\n") + "\n> Superseded by 0002 on 2026-11-08.\n",
    )
    assert diff(base, commit(repo)) == []


def test_two_appended_lines_are_not_one_note(repo: Path):
    base = branch_from_main(repo, {"decisions/0001-a.md": record("accepted")})
    on_branch(repo)
    write(repo, "decisions/0001-a.md", record("accepted") + "\n> One.\n> Two.\n")
    assert len(diff(base, commit(repo))) == 1


def test_a_note_after_a_file_without_trailing_newline_is_allowed(repo: Path):
    base = branch_from_main(repo, {"decisions/0001-a.md": record("accepted", body="No newline at the end")})
    on_branch(repo)
    write(repo, "decisions/0001-a.md", record("accepted", body="No newline at the end") + "\n\n> Note.\n")
    assert diff(base, commit(repo)) == []


def test_line_ending_changes_alone_are_not_changes(repo: Path):
    text = record("accepted")
    base = branch_from_main(repo, {"decisions/0001-a.md": text})
    on_branch(repo)
    write(repo, "decisions/0001-a.md", text.replace("\n", "\r\n"))
    assert diff(base, commit(repo)) == []


@pytest.mark.parametrize(
    ("before", "after", "allowed"),
    [
        ("proposed", "accepted", True),
        ("proposed", "rejected", True),
        ("accepted", "superseded", True),
        ("rejected", "superseded", True),
        ("accepted", "proposed", False),
        ("accepted", "rejected", False),
        ("rejected", "proposed", False),
        ("rejected", "accepted", False),
        ("superseded", "accepted", False),
        ("superseded", "proposed", False),
    ],
)
def test_status_transitions(repo: Path, before: str, after: str, allowed: bool):
    extra_before = "superseded_by: 0002-b.md\n" if before == "superseded" else ""
    extra_after = "superseded_by: 0002-b.md\n" if after == "superseded" else ""
    base = branch_from_main(
        repo,
        {"decisions/0001-a.md": record(before, extra=extra_before), "decisions/0002-b.md": record("accepted")},
    )
    on_branch(repo)
    write(repo, "decisions/0001-a.md", record(after, extra=extra_after))
    found = diff(base, commit(repo))
    assert (found == []) is allowed
    if not allowed:
        assert any(f"status {before} -> {after} is not allowed" in line for line in found)


def test_a_frozen_record_cannot_be_deleted_but_a_proposed_one_can(repo: Path):
    base = branch_from_main(
        repo,
        {
            "decisions/0001-a.md": record("accepted"),
            "decisions/0002-b.md": record("rejected"),
            "decisions/0003-c.md": record("proposed"),
        },
    )
    on_branch(repo)
    for name in ("0001-a", "0002-b", "0003-c"):
        (repo / f"decisions/{name}.md").unlink()
    found = diff(base, commit(repo))
    assert sorted(found) == [
        "decisions/0001-a.md: the accepted record is deleted",
        "decisions/0002-b.md: the rejected record is deleted",
    ]


def test_renaming_a_frozen_record_is_found_even_when_the_file_is_short(repo: Path):
    short = "---\nstatus: accepted\n---\nx\n"
    base = branch_from_main(repo, {"decisions/0001-a.md": short, "decisions/0002-b.md": record("proposed")})
    on_branch(repo)
    git(repo, "mv", "decisions/0001-a.md", "decisions/0001-renamed.md")
    git(repo, "mv", "decisions/0002-b.md", "decisions/0002-renamed.md")
    found = diff(base, commit(repo))
    assert found == ["decisions/0001-renamed.md: the accepted record is renamed from decisions/0001-a.md"]


def test_new_records_may_start_accepted_or_rejected_but_not_superseded(repo: Path):
    base = branch_from_main(repo, {"decisions/0001-a.md": record("accepted")})
    on_branch(repo)
    write(repo, "decisions/0002-b.md", record("accepted"))
    write(repo, "decisions/0003-c.md", record("rejected"))
    write(repo, "decisions/0004-d.md", record("superseded", extra="superseded_by: 0002-b.md\n"))
    assert diff(base, commit(repo)) == ["decisions/0004-d.md: a new record cannot start as superseded"]


def test_a_base_that_accepted_the_record_after_the_branch_point_still_protects_it(repo: Path):
    branch_from_main(repo, {"decisions/0001-a.md": record("proposed")})
    on_branch(repo)
    write(repo, "decisions/0001-a.md", record("proposed", body="Edited on the branch.\n"))
    head = commit(repo)
    git(repo, "switch", "-q", "main")
    write(repo, "decisions/0001-a.md", record("accepted"))
    base = commit(repo, "accept")
    assert any("the accepted record changes beyond" in line for line in diff(base, head))


def test_a_branch_that_did_not_touch_a_record_is_not_blamed_for_the_base_moving_it(repo: Path):
    branch_from_main(repo, {"decisions/0001-a.md": record("proposed")})
    on_branch(repo)
    write(repo, "decisions/0002-b.md", record("proposed"))
    head = commit(repo)
    git(repo, "switch", "-q", "main")
    write(repo, "decisions/0001-a.md", record("accepted"))
    base = commit(repo, "accept")
    assert diff(base, head) == []


def test_a_merge_result_is_checked_against_its_first_parent(repo: Path):
    branch_from_main(repo, {"decisions/0001-a.md": record("proposed")})
    on_branch(repo)
    write(repo, "decisions/0001-a.md", record("proposed", body="Edited on the branch.\n"))
    commit(repo)
    git(repo, "switch", "-q", "main")
    write(repo, "decisions/0002-b.md", record("accepted"))
    write(repo, "decisions/0001-a.md", record("accepted"))
    commit(repo, "accept")
    git(repo, "switch", "-q", "-c", "merged", "main")
    git(repo, "merge", "-q", "--no-edit", "-X", "theirs", "pr")
    assert any("0001-a.md: the accepted record changes beyond" in line for line in diff("HEAD^1", "HEAD"))


def test_an_unknown_base_is_skipped(repo: Path, capsys: pytest.CaptureFixture[str]):
    head = branch_from_main(repo, {"decisions/0001-a.md": record("accepted")})
    assert diff("0" * 40, head) == []
    assert "not in this clone" in capsys.readouterr().err


# --- the files that `parity docs render` copies ---------------------------------------------

SHARED = SCRIPT.parents[1]


def test_the_script_has_nothing_that_looks_like_a_docs_placeholder():
    from parity import docs

    assert docs._PLACEHOLDER.findall(SCRIPT.read_text(encoding="utf-8")) == []


def test_the_workflow_template_has_only_the_tracked_branch_placeholder():
    import yaml

    from parity import docs

    text = (SHARED / "workflows" / "decision-records.core.yml").read_text(encoding="utf-8")
    assert set(docs._PLACEHOLDER.findall(text)) == {"{{OWN_REF}}"}
    workflow = yaml.safe_load(text.replace("{{OWN_REF}}", "develop"))
    assert workflow[True]["push"]["branches"] == ["develop"]
    commands = [step["run"] for step in workflow["jobs"]["records"]["steps"] if "run" in step]
    assert all("decision_records.py" in command for command in commands)
    assert len(commands) == 3


# --- command line --------------------------------------------------------------------------


def test_main_exit_status_and_output(repo: Path, capsys: pytest.CaptureFixture[str]):
    write(repo, "decisions/0001-a.md", record("draft"))
    assert dr.main(["decision_records.py", "check"]) == 1
    assert "decision-records: decisions/0001-a.md: status 'draft'" in capsys.readouterr().out
    write(repo, "decisions/0001-a.md", record("proposed"))
    assert dr.main(["decision_records.py", "check"]) == 0
    assert capsys.readouterr().out == ""


def test_main_uses_the_folder_option(repo: Path):
    write(repo, "doc/decisions/0001-a.md", record("draft"))
    assert dr.main(["decision_records.py", "check", "--dir", "doc/decisions/"]) == 1
    assert dr.main(["decision_records.py", "check"]) == 0


def test_version(capsys: pytest.CaptureFixture[str]):
    with pytest.raises(SystemExit) as stop:
        dr.main(["decision_records.py", "--version"])
    assert stop.value.code == 0
    assert capsys.readouterr().out.strip() == dr.VERSION
