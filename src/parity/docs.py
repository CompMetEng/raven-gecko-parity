"""Composing the shared GitHub docs and writing them into the siblings.

Four repos otherwise mean four hand-maintained copies of CONTRIBUTING.md, AGENTS.md, issue
and PR templates, CI workflows, and so on, drifting the moment one of them is edited and
not the others. The shared source lives here (`docs/shared/`); each document is a core file
plus zero or more per-repo appendices, each substituted into its own `{{NAME}}` placeholder,
keyed by the repo's short name (`RAVEN`, `GECKO`, `raven-toolbox`, `geckopy`) so a MATLAB and
a Python side never share an appendix by accident. Most documents need only one substitution
point (`{{APPENDIX}}`); a workflow file that differs from its sibling in more than one place
(e.g. an optional extra checkout step, and a separate test-command block) declares one
placeholder per substitution instead of forcing everything through a single spot. Seven
placeholders are always available, filled in from `parity.toml` itself rather than an
appendix file: `{{FULL_REPO}}` (this repo's own `org/repo` slug, for issue search/Discussions/
security-advisory links), `{{REPO}}` (just the short name, e.g. `raven-toolbox`, for display
text), `{{PAIR}}`/`{{PAIR_UPPER}}` (the pair name, e.g. `raven`/`RAVEN` ---
matches the `PARITY_<PAIR>_<SIDE>_PATH` environment-variable convention `parity.toml` uses),
`{{SIBLING_REPO}}`/`{{SIBLING_REF}}` (the *other* side's `org/repo` and tracked branch; for
a workflow that needs to check its sibling out, e.g. GECKO's CI checking out geckopy),
`{{OWN_REF}}` (this side's own tracked branch, e.g. RAVEN's `develop`), and `{{PACKAGE}}`
(the Python import package name; empty on the matlab side).

A `DocSpec` normally applies to every repo on a matching side; set `repos=(...)` on it when
something is genuinely one repo's own concern rather than every repo of that language (e.g. a
stand-in check for geckopy specifically, not "every Python repo").
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from parity.config import Config, Pair

_PLACEHOLDER = re.compile(r"\{\{[A-Z_]+\}\}")


class DocsError(RuntimeError):
    """Raised when a shared doc source or a repo's appendix is missing, or a
    placeholder is left unresolved."""


@dataclass(frozen=True)
class DocSpec:
    #: Path written into the sibling repo, relative to its root.
    output: str
    #: Path read under `docs/shared/`.
    source: str
    #: Placeholder name -> appendix file suffix (read from
    #: `docs/shared/appendices/<repo><suffix>`). Empty for a document with no per-repo
    #: section.
    appendices: dict[str, str] = None  # type: ignore[assignment]
    #: Which sides this document applies to: both by default. A MATLAB-only or
    #: Python-only workflow file (the two languages need structurally different CI, not
    #: just a substituted detail) is scoped with e.g. ``sides=("matlab",)``; two specs may
    #: then share the same `output` path as long as their `sides` never overlap.
    sides: tuple[str, ...] = ("matlab", "python")
    #: Which repos (short name) this document applies to: every repo on a matching side
    #: by default. Set for something only one repo needs (e.g. geckopy has no Sphinx docs
    #: to build yet, so a stand-in Markdown check is geckopy-only, not "every Python repo").
    repos: tuple[str, ...] | None = None

    def __post_init__(self):
        if self.appendices is None:
            object.__setattr__(self, "appendices", {})


DOCS: tuple[DocSpec, ...] = (
    DocSpec(".github/CONTRIBUTING.md", "CONTRIBUTING.core.md", appendices={"APPENDIX": ".md"}),
    DocSpec("AGENTS.md", "AGENTS.core.md", appendices={"APPENDIX": ".md"}),
    DocSpec(".github/pull_request_template.md", "PULL_REQUEST_TEMPLATE.core.md"),
    DocSpec(
        ".github/ISSUE_TEMPLATE/bug_report.md",
        "ISSUE_TEMPLATE/bug_report.core.md",
        appendices={"APPENDIX": ".issue-env.md"},
    ),
    DocSpec(".github/ISSUE_TEMPLATE/feature_request.md", "ISSUE_TEMPLATE/feature_request.core.md"),
    DocSpec(".github/ISSUE_TEMPLATE/config.yml", "ISSUE_TEMPLATE/config.core.yml"),
    DocSpec(".github/CODE_OF_CONDUCT.md", "CODE_OF_CONDUCT.core.md"),
    DocSpec(".github/SECURITY.md", "SECURITY.core.md"),
    DocSpec(
        ".github/workflows/tests.yml",
        "workflows/matlab-tests.core.yml",
        appendices={"EXTRA_CHECKOUT": ".extra-checkout.yml", "TEST_COMMAND": ".test-command.txt"},
        sides=("matlab",),
    ),
    DocSpec(
        ".github/workflows/tests.yml",
        "workflows/python-tests.core.yml",
        appendices={"INSTALL_EXTRAS": ".install-extras.txt", "PYTEST_ARGS": ".pytest-args.txt"},
        sides=("python",),
    ),
    DocSpec(
        ".github/workflows/parity.yml",
        "workflows/matlab-parity.core.yml",
        sides=("matlab",),
    ),
    DocSpec(
        ".github/workflows/lint.yml",
        "workflows/matlab-lint.core.yml",
        appendices={"EXCLUDE_DIRS": ".exclude-dirs.txt"},
        sides=("matlab",),
    ),
    DocSpec(
        ".github/workflows/lint.yml",
        "workflows/python-lint.core.yml",
        sides=("python",),
    ),
    DocSpec(
        ".github/workflows/documentation.yml",
        "workflows/matlab-docs.core.yml",
        appendices={"EXTRA_CHECKOUT": ".extra-checkout.yml", "DOC_COMMAND": ".doc-command.txt"},
        sides=("matlab",),
    ),
    DocSpec(
        ".github/workflows/install-smoke-test.yml",
        "workflows/python-install-smoke-test.core.yml",
        sides=("python",),
    ),
    DocSpec(
        ".github/workflows/release.yml",
        "workflows/python-release.core.yml",
        sides=("python",),
    ),
    DocSpec(
        ".github/workflows/release.yml",
        "workflows/matlab-release.core.yml",
        sides=("matlab",),
    ),
)


@dataclass(frozen=True)
class Target:
    pair: str
    side: str
    repo: str  # short name, e.g. "RAVEN"; keys the appendix files
    full_repo: str  # "org/repo", e.g. "org/RAVEN"; used in links
    sibling_repo: str  # "org/repo" of the *other* side of this pair
    sibling_ref: str  # the branch the other side is tracked at, per parity.toml
    own_ref: str  # this side's own tracked branch, per parity.toml
    package: str  # the Python import package name; "" on the matlab side
    path: Path


@dataclass(frozen=True)
class Drift:
    target: Target
    filename: str
    status: str  # "missing" or "stale"


def _short_name(repo: str) -> str:
    return repo.rsplit("/", 1)[-1]


def _shared_root(config: Config) -> Path:
    return config.root / "docs" / "shared"


def _read(path: Path) -> str:
    if not path.is_file():
        raise DocsError(f"missing {path}")
    return path.read_text(encoding="utf-8")


def _pairs(config: Config, requested: str | None) -> list[Pair]:
    if requested:
        return [config.pair(requested)]
    return [config.pairs[name] for name in sorted(config.pairs)]


def targets(config: Config, requested: str | None) -> list[Target]:
    """One entry per side of each selected pair: what gets rendered, and where."""
    out = []
    for pair in _pairs(config, requested):
        for side in ("matlab", "python"):
            side_config = pair.side(side)
            other_config = pair.side("python" if side == "matlab" else "matlab")
            out.append(
                Target(
                    pair=pair.name,
                    side=side,
                    repo=_short_name(side_config.repo),
                    full_repo=side_config.repo,
                    sibling_repo=other_config.repo,
                    sibling_ref=other_config.ref,
                    own_ref=side_config.ref,
                    package=getattr(side_config, "package", ""),
                    path=side_config.path,
                )
            )
    return out


def compose(config: Config, target: Target) -> dict[str, str]:
    """Render every document in `DOCS` for one target repo."""
    shared = _shared_root(config)

    documents = {}
    for spec in DOCS:
        if target.side not in spec.sides:
            continue
        if spec.repos is not None and target.repo not in spec.repos:
            continue
        text = _read(shared / spec.source)
        for placeholder, suffix in spec.appendices.items():
            appendix_path = shared / "appendices" / f"{target.repo}{suffix}"
            text = text.replace("{{" + placeholder + "}}", _read(appendix_path).rstrip("\n"))
        text = text.replace("{{FULL_REPO}}", target.full_repo)
        text = text.replace("{{REPO}}", target.repo)
        text = text.replace("{{PAIR}}", target.pair)
        text = text.replace("{{PAIR_UPPER}}", target.pair.upper())
        text = text.replace("{{SIBLING_REPO}}", target.sibling_repo)
        text = text.replace("{{SIBLING_REF}}", target.sibling_ref)
        text = text.replace("{{OWN_REF}}", target.own_ref)
        text = text.replace("{{PACKAGE}}", target.package)

        leftover = _PLACEHOLDER.findall(text)
        if leftover:
            raise DocsError(f"{spec.output} for {target.repo} still contains {leftover!r}")
        documents[spec.output] = text
    return documents


def render(config: Config, requested: str | None) -> list[tuple[Target, dict[str, Path]]]:
    """Write the composed docs into each target's sibling checkout.

    Only ever writes under the sibling paths `parity.toml` configures, never this repo's
    own tree, other than reading the shared source under it.
    """
    written = []
    for target in targets(config, requested):
        documents = compose(config, target)
        paths = {}
        for output, text in documents.items():
            out = target.path / output
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_text(text, encoding="utf-8")
            paths[output] = out
        written.append((target, paths))
    return written


def check(config: Config, requested: str | None) -> list[Drift]:
    """Compare what's on disk in each sibling against the current rendered output."""
    drift = []
    for target in targets(config, requested):
        documents = compose(config, target)
        for output, text in documents.items():
            out = target.path / output
            if not out.is_file():
                drift.append(Drift(target, output, "missing"))
            elif out.read_text(encoding="utf-8") != text:
                drift.append(Drift(target, output, "stale"))
    return drift
