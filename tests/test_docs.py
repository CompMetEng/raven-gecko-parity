"""Composing the shared GitHub docs once and writing them into each sibling.

Four repos otherwise mean four hand-maintained copies drifting the moment one is edited and
not the others. These tests build a fake shared source and a fake pair (never the real
RAVEN/GECKO/raven-toolbox/geckopy checkouts) so they say something about the composition and
drift-detection logic itself.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from parity import docs
from parity.cli import main
from parity.config import Config, MatlabSide, Pair, PythonSide


def _write_shared_source(root: Path) -> None:
    shared = root / "docs" / "shared"
    (shared / "appendices").mkdir(parents=True)
    (shared / "ISSUE_TEMPLATE").mkdir()

    (shared / "CONTRIBUTING.core.md").write_text("# Contributing\n\n{{APPENDIX}}\n", encoding="utf-8")
    (shared / "AGENTS.core.md").write_text("# Agents\n\n{{APPENDIX}}\n", encoding="utf-8")
    (shared / "PULL_REQUEST_TEMPLATE.core.md").write_text("### Checklist\n", encoding="utf-8")
    (shared / "CODE_OF_CONDUCT.core.md").write_text("# Code of Conduct\n", encoding="utf-8")
    (shared / "SECURITY.core.md").write_text(
        "# Security\n\nReport at https://github.com/{{FULL_REPO}}/security/advisories/new\n", encoding="utf-8"
    )
    (shared / "ISSUE_TEMPLATE" / "bug_report.core.md").write_text(
        "---\nname: Bug report\n---\n\n### Environment\n\n{{APPENDIX}}\n\nSearch: "
        "https://github.com/{{FULL_REPO}}/issues\n",
        encoding="utf-8",
    )
    (shared / "ISSUE_TEMPLATE" / "feature_request.core.md").write_text(
        "---\nname: Feature request\n---\n\nSearch: https://github.com/{{FULL_REPO}}/issues\n", encoding="utf-8"
    )
    (shared / "ISSUE_TEMPLATE" / "config.core.yml").write_text(
        "blank_issues_enabled: false\ncontact_links:\n  - url: https://github.com/{{FULL_REPO}}/discussions\n",
        encoding="utf-8",
    )

    (shared / "workflows").mkdir()
    (shared / "workflows" / "matlab-tests.core.yml").write_text(
        "name: Tests\nsteps:\n{{EXTRA_CHECKOUT}}\ncommand: |\n{{TEST_COMMAND}}\n", encoding="utf-8"
    )
    (shared / "workflows" / "python-tests.core.yml").write_text(
        "name: Tests\ninstall: {{INSTALL_EXTRAS}}\nargs: {{PYTEST_ARGS}}\n", encoding="utf-8"
    )
    (shared / "workflows" / "matlab-parity.core.yml").write_text(
        "pair: {{PAIR}}\nenv_prefix: {{PAIR_UPPER}}\nsibling: {{SIBLING_REPO}}@{{SIBLING_REF}}\n",
        encoding="utf-8",
    )
    (shared / "workflows" / "matlab-lint.core.yml").write_text(
        "name: Lint\nexclude: {{EXCLUDE_DIRS}}\n", encoding="utf-8"
    )
    (shared / "workflows" / "python-lint.core.yml").write_text("name: Lint\nmypy: src\n", encoding="utf-8")
    (shared / "workflows" / "matlab-docs.core.yml").write_text(
        "name: Documentation\nbranches: [{{OWN_REF}}, develop, main]\nsteps:\n{{EXTRA_CHECKOUT}}\ncommand: |\n"
        "{{DOC_COMMAND}}\n",
        encoding="utf-8",
    )
    (shared / "workflows" / "python-install-smoke-test.core.yml").write_text(
        "name: Install smoke test\nimport: {{PACKAGE}}\n", encoding="utf-8"
    )
    (shared / "workflows" / "python-release.core.yml").write_text(
        "name: Release\npackage: {{PACKAGE}}\nrepo: {{REPO}}\n", encoding="utf-8"
    )
    (shared / "workflows" / "matlab-release.core.yml").write_text(
        "name: Release\nbranch: {{OWN_REF}}\n", encoding="utf-8"
    )

    (shared / "appendices" / "FAKE.md").write_text("- fake matlab specifics\n", encoding="utf-8")
    (shared / "appendices" / "fakepy.md").write_text("- fake python specifics\n", encoding="utf-8")
    (shared / "appendices" / "FAKE.issue-env.md").write_text("- MATLAB version: ...\n", encoding="utf-8")
    (shared / "appendices" / "fakepy.issue-env.md").write_text("- Python version: ...\n", encoding="utf-8")
    (shared / "appendices" / "FAKE.extra-checkout.yml").write_text("", encoding="utf-8")
    (shared / "appendices" / "FAKE.test-command.txt").write_text("  fake matlab test command", encoding="utf-8")
    (shared / "appendices" / "fakepy.install-extras.txt").write_text("dev", encoding="utf-8")
    (shared / "appendices" / "fakepy.pytest-args.txt").write_text("--fake-arg", encoding="utf-8")
    (shared / "appendices" / "FAKE.exclude-dirs.txt").write_text("{'vendored'}", encoding="utf-8")
    (shared / "appendices" / "FAKE.doc-command.txt").write_text("  fake_doc_build();", encoding="utf-8")


@pytest.fixture
def config(tmp_path: Path) -> Config:
    _write_shared_source(tmp_path)
    matlab_root = tmp_path / "FAKE"
    python_root = tmp_path / "fakepy"
    matlab_root.mkdir()
    python_root.mkdir()
    pair = Pair(
        name="fake",
        ledger=tmp_path / "ledgers" / "fake.yml",
        matlab=MatlabSide(repo="org/FAKE", path=matlab_root),
        python=PythonSide(repo="org/fakepy", path=python_root, package="fakepy"),
    )
    return Config(root=tmp_path, pairs={"fake": pair})


def _matlab_target(config: Config) -> docs.Target:
    return next(t for t in docs.targets(config, None) if t.side == "matlab")


def _python_target(config: Config) -> docs.Target:
    return next(t for t in docs.targets(config, None) if t.side == "python")


def test_compose_substitutes_the_repo_appendix(config: Config):
    documents = docs.compose(config, _matlab_target(config))
    assert "fake matlab specifics" in documents[".github/CONTRIBUTING.md"]
    assert "fake matlab specifics" in documents["AGENTS.md"]
    assert "{{APPENDIX}}" not in documents[".github/CONTRIBUTING.md"]
    assert "{{APPENDIX}}" not in documents["AGENTS.md"]


def test_compose_substitutes_the_full_repo_slug(config: Config):
    documents = docs.compose(config, _matlab_target(config))
    assert "https://github.com/org/FAKE/discussions" in documents[".github/ISSUE_TEMPLATE/config.yml"]
    assert "https://github.com/org/FAKE/issues" in documents[".github/ISSUE_TEMPLATE/bug_report.md"]


def test_compose_picks_the_matlab_workflow_template_for_a_matlab_side(config: Config):
    documents = docs.compose(config, _matlab_target(config))
    tests_yml = documents[".github/workflows/tests.yml"]
    assert "fake matlab test command" in tests_yml
    assert "install:" not in tests_yml


def test_compose_picks_the_python_workflow_template_for_a_python_side(config: Config):
    documents = docs.compose(config, _python_target(config))
    tests_yml = documents[".github/workflows/tests.yml"]
    assert "install: dev" in tests_yml
    assert "args: --fake-arg" in tests_yml
    assert "fake matlab test command" not in tests_yml


def test_compose_substitutes_pair_and_sibling_placeholders(config: Config):
    documents = docs.compose(config, _matlab_target(config))
    parity_yml = documents[".github/workflows/parity.yml"]
    assert "pair: fake" in parity_yml
    assert "env_prefix: FAKE" in parity_yml
    assert "sibling: org/fakepy@main" in parity_yml


def test_python_side_gets_no_parity_workflow(config: Config):
    assert ".github/workflows/parity.yml" not in docs.compose(config, _python_target(config))


def test_compose_fails_loudly_on_a_missing_appendix(config: Config):
    """A repo without an appendix file must error, not silently ship a document with the
    literal placeholder text still in it."""
    bogus = docs.Target(
        pair="fake",
        side="matlab",
        repo="nonexistent",
        full_repo="org/nonexistent",
        sibling_repo="org/nonexistent-py",
        sibling_ref="main",
        own_ref="main",
        package="",
        path=Path("."),
    )
    with pytest.raises(docs.DocsError):
        docs.compose(config, bogus)


def test_render_writes_every_document_into_both_sides(config: Config):
    written = docs.render(config, None)
    assert len(written) == 2  # matlab side + python side of the one pair

    for target, paths in written:
        expected = set(docs.compose(config, target))
        assert set(paths) == expected
        for output in expected:
            assert paths[output].is_file()

    matlab_path = config.pairs["fake"].matlab.path
    assert "fake matlab specifics" in (matlab_path / ".github" / "CONTRIBUTING.md").read_text(encoding="utf-8")
    assert (matlab_path / ".github" / "ISSUE_TEMPLATE" / "bug_report.md").is_file()
    assert "fake matlab test command" in (matlab_path / ".github" / "workflows" / "tests.yml").read_text(
        encoding="utf-8"
    )

    python_path = config.pairs["fake"].python.path
    assert "install: dev" in (python_path / ".github" / "workflows" / "tests.yml").read_text(encoding="utf-8")


def test_check_reports_everything_missing_before_a_render(config: Config):
    drift = docs.check(config, None)
    expected = sum(len(docs.compose(config, t)) for t in docs.targets(config, None))
    assert len(drift) == expected
    assert all(item.status == "missing" for item in drift)


def test_check_is_clean_immediately_after_a_render(config: Config):
    docs.render(config, None)
    assert docs.check(config, None) == []


def test_check_reports_stale_after_a_hand_edit(config: Config):
    docs.render(config, None)
    agents_md = config.pairs["fake"].matlab.path / "AGENTS.md"
    agents_md.write_text(agents_md.read_text(encoding="utf-8") + "\nhand edit\n", encoding="utf-8")

    drift = docs.check(config, None)
    assert len(drift) == 1
    assert drift[0].filename == "AGENTS.md"
    assert drift[0].status == "stale"


CLI_CONFIG = """
[pairs.fake]
ledger = "ledgers/fake.yml"

[pairs.fake.matlab]
repo = "org/FAKE"
path = "FAKE"

[pairs.fake.python]
repo = "org/fakepy"
package = "fakepy"
path = "fakepy"
"""


def _cli_config(tmp_path: Path) -> Path:
    _write_shared_source(tmp_path)
    (tmp_path / "FAKE").mkdir()
    (tmp_path / "fakepy").mkdir()
    path = tmp_path / "parity.toml"
    path.write_text(CLI_CONFIG, encoding="utf-8")
    return path


def test_cli_docs_render(tmp_path: Path, capsys):
    config_path = _cli_config(tmp_path)
    assert main(["-c", str(config_path), "docs", "render"]) == 0
    assert (tmp_path / "FAKE" / ".github" / "CONTRIBUTING.md").is_file()
    assert "wrote" in capsys.readouterr().out


def test_cli_docs_check_exits_nonzero_before_a_render(tmp_path: Path):
    config_path = _cli_config(tmp_path)
    assert main(["-c", str(config_path), "docs", "check"]) == 1


def test_cli_docs_check_exits_zero_after_a_render(tmp_path: Path):
    config_path = _cli_config(tmp_path)
    main(["-c", str(config_path), "docs", "render"])
    assert main(["-c", str(config_path), "docs", "check"]) == 0
