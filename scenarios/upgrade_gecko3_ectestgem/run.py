"""Python side of the GECKO 3 upgrade scenario.

Every case builds its own state in a scratch folder, runs the call and records what it left
behind and whether it raised. A raise is a value (`raised`). Nothing reaches the network.
"""

import shutil
import warnings
from pathlib import Path

from geckopy import ModelAdapter, load_conventional_gem, make_ec_model, save_ec_model
from geckopy.upgrade import convert_custom_kcats, upgrade_ec_model, upgrade_project
from geckopy.utilities.load_ec_model import load_ec_model

POOL = "prot_pool"
EXCHANGE = "prot_pool_exchange"


def _copy_project(template: str, target: Path) -> Path:
    if target.exists():
        shutil.rmtree(target)
    shutil.copytree(template, target, ignore=shutil.ignore_patterns("__pycache__"))
    return target


def _model(project: Path):
    adapter = ModelAdapter.from_folder(project)
    conv = load_conventional_gem(adapter)
    return make_ec_model(conv, adapter), adapter


def _protein_rxns(model):
    return [r for r in model.reactions if r.id.startswith("usage_prot_") or r.id == EXCHANGE]


def _summary(model) -> list:
    """Pool coefficient, bounds and objective coefficient of the protein usage reactions."""
    pool = model.metabolites.get_by_id(POOL)
    rows = []
    for rxn in sorted(_protein_rxns(model), key=lambda r: r.id):
        rows.append({"rxn": rxn.id, "pool_coef": float(rxn.metabolites.get(pool, 0.0)),
                     "lb": float(rxn.lower_bound), "ub": float(rxn.upper_bound),
                     "c": float(rxn.objective_coefficient)})
    return rows


def _make_reverse(model) -> None:
    """The GECKO 3 shape: usage and pool exchange reactions reverse, with objective weights."""
    for rxn in _protein_rxns(model):
        rxn.add_metabolites({m: -2 * c for m, c in rxn.metabolites.items()})
        lb, ub = rxn.bounds
        rxn.bounds = (-ub, -lb)
    model.reactions.get_by_id(EXCHANGE).objective_coefficient = 2.0
    model.reactions.get_by_id("usage_prot_P1").objective_coefficient = -1.0


def _make_forward_with_negative_lb(model) -> None:
    """A GECKO 4 model whose usage reaction of P1 has a negative lower bound."""
    model.reactions.get_by_id("usage_prot_P1").bounds = (-5.0, 1000.0)
    model.reactions.get_by_id(EXCHANGE).objective_coefficient = 2.0


def _attempt(fn):
    try:
        return False, fn()
    except Exception:  # noqa: BLE001 --- whether the call raises is the finding, not the type
        return True, None


# --------------------------------------------------------------------------- #
# in-memory and file round trip
# --------------------------------------------------------------------------- #

def _case_in_memory(inputs, root: Path) -> dict:
    out = {}
    project = _copy_project(inputs["project_template"], root / "memory")
    for name, make in (("reverse", _make_reverse), ("forward_negative_lb", _make_forward_with_negative_lb)):
        model, _ = _model(project)
        make(model)
        before = _summary(model)
        changed = upgrade_ec_model(model)
        after = _summary(model)
        again = upgrade_ec_model(model)
        out[name] = {"before": before, "changed": sorted(changed), "after": after,
                     "changed_again": sorted(again), "after_again": _summary(model)}
    return out


def _case_yaml(inputs, root: Path) -> dict:
    out = {}
    project = _copy_project(inputs["project_template"], root / "yaml")
    for name, make in (("reverse", _make_reverse), ("forward_negative_lb", _make_forward_with_negative_lb)):
        model, _ = _model(project)
        make(model)
        path = project / "models" / f"{name}.yml"
        save_ec_model(model, path)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            loaded = load_ec_model(path)
        # What the reader returns and which reactions the upgrade turns are not recorded: the
        # two readers differ in whether they turn a reaction with a negative lower bound.
        upgrade_ec_model(loaded)
        out[name] = {"after": _summary(loaded)}
    return out


# --------------------------------------------------------------------------- #
# projects
# --------------------------------------------------------------------------- #

def _listing(project: Path) -> list:
    today = __import__("datetime").date.today().isoformat()
    names = []
    for p in project.rglob("*"):
        if p.is_file():
            names.append(p.relative_to(project).as_posix().replace(today, "<date>"))
    return sorted(names)


def _rows(report) -> list:
    # MATLAB names the argument 'ecModel', geckopy ec_model.
    return [{"step": str(r.step), "status": str(r.status),
             "note": str(r.note).replace("(ec_model)", "(<ecModel argument>)")} for r in report.itertuples()]


def _records(path: Path) -> list:
    """Rows of a kcat source file: rxn, enzymes, kcat, source and note."""
    if not path.is_file():
        return []
    rows = []
    header = None
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("#") or not line:
            continue
        cells = line.split("\t")
        if header is None:
            header = cells
            continue
        rec = dict(zip(header, cells + [""] * (len(header) - len(cells))))
        rows.append({k: rec[k] for k in ("rxn", "enzymes", "kcat", "source", "note")})
    return sorted(rows, key=lambda r: (r["rxn"], r["enzymes"], r["kcat"]))


def _project_case(inputs, root: Path, name: str, *, models=(), no_toml=False, dry_run=True, with_model=True,
                  twice=False, dlkcat="valid") -> dict:
    project = _copy_project(inputs["project_template"], root / name)
    if dlkcat == "valid":
        shutil.copy(inputs["dlkcat_valid"], project / "data" / "DLKcat.tsv")
    model, adapter = _model(project)
    for model_name, make in models:
        m, _ = _model(project)
        make(m)
        save_ec_model(m, project / "models" / f"{model_name}.yml")
    if no_toml:
        (project / "model_adapter.toml").unlink()
        shutil.copy(inputs["class_adapter"], project)
    rec = {}
    raised, report = _attempt(lambda: upgrade_project(project, dry_run=dry_run, ec_model=model if with_model else None))
    rec["raised"] = raised
    rec["rows"] = [] if raised else _rows(report)
    if twice and not raised:
        raised2, report2 = _attempt(lambda: upgrade_project(project, dry_run=dry_run, ec_model=model))
        rec["rows_second_run"] = [] if raised2 else _rows(report2)
    rec["files"] = _listing(project)
    rec["custom"] = _records(project / "data" / "kcats" / "custom.tsv")
    rec["dlkcat"] = _records(project / "data" / "kcats" / "dlkcat.tsv")
    toml = project / "model_adapter.toml"
    rec["toml_first_line"] = toml.read_text(encoding="utf-8").splitlines()[0] if toml.is_file() else ""
    return rec


def _case_projects(inputs, root: Path) -> dict:
    out = {}
    out["dry_run_without_model"] = _project_case(inputs, root, "dry_nomodel", with_model=False,
                                                 models=[("gecko3", _make_reverse)])
    out["dry_run_with_model"] = _project_case(inputs, root, "dry_model", models=[("gecko3", _make_reverse)])
    out["apply_twice"] = _project_case(inputs, root, "apply", models=[("gecko3", _make_reverse)],
                                       dry_run=False, twice=True)
    out["forward_model_negative_lb"] = _project_case(inputs, root, "forward_lb",
                                                     models=[("forward", _make_forward_with_negative_lb)],
                                                     dry_run=False)
    out["stored_dlkcat_foreign_substrates"] = _project_case(inputs, root, "foreign", dry_run=False,
                                                            dlkcat="foreign")
    out["no_toml"] = _project_case(inputs, root, "no_toml", no_toml=True)
    return out


# --------------------------------------------------------------------------- #
# convertCustomKcats
# --------------------------------------------------------------------------- #

def _case_convert(inputs, root: Path) -> dict:
    out = {}
    for name, key in (("stock", "custom_stock"), ("empty_note_and_duplicate", "custom_edge")):
        project = _copy_project(inputs["project_template"], root / f"convert_{name}")
        model, adapter = _model(project)
        target = project / "data" / "kcats" / "custom.tsv"
        raised, result = _attempt(lambda: convert_custom_kcats(
            model, custom_kcats=inputs[key], out_file=target, adapter=adapter))
        out[name] = {"raised": raised, "rows": 0 if raised else int(result[1]["rows"]),
                     "records": _records(target)}
    return out


def run(ctx):
    inputs = ctx["inputs"]
    root = Path(inputs["scratch"])
    return {
        "in_memory": _case_in_memory(inputs, root),
        "yaml_round_trip": _case_yaml(inputs, root),
        "projects": _case_projects(inputs, root),
        "convert_custom": _case_convert(inputs, root),
    }
