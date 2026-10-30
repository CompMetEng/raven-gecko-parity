"""Python side of the kcat assignment pipeline scenario.

Every case builds its own project folder in a scratch folder, runs the call and records what
it returned and whether it raised. A raise is a value (`raised`). Nothing reaches the network.
"""

import shutil
from pathlib import Path

import numpy as np
import pandas as pd

from geckopy import ModelAdapter, load_conventional_gem, make_ec_model
from geckopy.databases import load_brenda_measurements
from geckopy.gather_kcats import assign_kcats, collect_kcats, kcats_from_ec_tree
from geckopy.gather_kcats.ec_tree import ECTree
from geckopy.gather_kcats.kcat_assignment import kcat_assignment_config, okp_methods
from geckopy.get_enzyme_data.ec_from_gem import fill_eccodes_from_gem

SELECT = '''
[kcat_assignment]
isozymes = false
standard = false

[[kcat_assignment.stage]]
action = "select"
sources = ["custom", "database_top", "dlkcat", "database_bottom"]
top_origin_limit = 4
'''
SELECT_FILL = SELECT.replace("isozymes = false", "isozymes = true").replace("standard = false", "standard = true")
COMBINE = '''
[kcat_assignment]
preset = "ec_tree_homologs"
isozymes = false
standard = false
'''
CONFIGS = {
    "default": "",
    "preset_gecko_classic": '[kcat_assignment]\npreset = "gecko_classic"\n',
    "two_stages": ('[[kcat_assignment.stage]]\naction = "select"\nsources = ["custom", "database_top"]\n'
                   'top_origin_limit = 1\n\n[[kcat_assignment.stage]]\naction = "combine"\n'
                   'sources = ["ec_tree"]\ncriteria = "median"\n'),
    "unknown_key": "[kcat_assignment]\nbogus = 1\n",
    "stage_sources_string": '[[kcat_assignment.stage]]\nsources = "custom"\n',
    "preset_empty": '[kcat_assignment]\npreset = ""\n',
    "criteria_capitalised": '[[kcat_assignment.stage]]\nsources = ["custom"]\ncriteria = "Median"\n',
    "okp_methods_unknown": '[kcat_sources.okp]\nmethods = ["Nope"]\n',
    "okp_methods_two": '[kcat_sources.okp]\nmethods = ["CataPro", "UniKP"]\n',
    "preset_and_stage": ('[kcat_assignment]\npreset = "okp"\n\n[[kcat_assignment.stage]]\n'
                         'sources = ["custom"]\n'),
    "brenda_aggregation_capitalised": '[kcat_sources.brenda]\naggregation = "Median"\n',
    "preset_unknown": '[kcat_assignment]\npreset = "nonsense"\n',
}


def _num(value) -> dict:
    value = float(value)
    if np.isnan(value):
        return {"nan": True, "value": 0.0}
    return {"nan": False, "value": value}


def _plain(value):
    """A config value as both sides can state it: float, str or list of str."""
    if isinstance(value, bool):
        return float(value)
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, (list, tuple)):
        return [str(v) for v in value]
    return str(value)


def _attempt(fn):
    try:
        return False, fn()
    except Exception:  # noqa: BLE001 --- whether the call raises is the finding, not the type
        return True, None


def _copy_project(template: str, target: Path) -> Path:
    if target.exists():
        shutil.rmtree(target)
    shutil.copytree(template, target, ignore=shutil.ignore_patterns("__pycache__"))
    return target


def _project(inputs: dict, root: Path, name: str, files: dict | None = None, toml: str = ""):
    """A project copy with the BRENDA fixtures, the named kcat source files and extra toml."""
    project = _copy_project(inputs["project_template"], root / name)
    brenda = project / "data" / "brenda"
    brenda.mkdir(exist_ok=True)
    for f in Path(inputs["brenda_dir"]).glob("*.tsv"):
        shutil.copy(f, brenda / f.name)
    kcats = project / "data" / "kcats"
    for target, source in (files or {}).items():
        kcats.mkdir(parents=True, exist_ok=True)
        shutil.copy(Path(inputs["kcats_dir"]) / f"{source}.tsv", kcats / f"{target}.tsv")
    if toml:
        path = project / "model_adapter.toml"
        path.write_text(path.read_text(encoding="utf-8") + "\n" + toml, encoding="utf-8")
    return project


def _model(project: Path):
    adapter = ModelAdapter.from_folder(project)
    model = make_ec_model(load_conventional_gem(adapter), adapter)
    fill_eccodes_from_gem(model)
    return model, adapter


# --------------------------------------------------------------------------- #
# kcat_assignment_config
# --------------------------------------------------------------------------- #

def _config_case(inputs: dict, root: Path, name: str) -> dict:
    project = _project(inputs, root, f"config_{name}", toml=CONFIGS[name])

    def read():
        adapter = ModelAdapter.from_folder(project)
        config = kcat_assignment_config(adapter)
        methods = okp_methods(config)
        return config, methods

    raised, result = _attempt(read)
    if raised:
        return {"raised": True, "preset": "", "isozymes": False, "standard": False, "stages": [],
                "sources": [], "okp_methods": []}
    config, methods = result
    stages = [{"action": s.action, "sources": [str(x) for x in s.sources], "criteria": s.criteria,
               "top_origin_limit": float(s.top_origin_limit), "bottom_origin_limit": float(s.bottom_origin_limit),
               "wildcard_limit": float(s.wildcard_limit)} for s in config.stages]
    sources = [{"source": src, "option": key, "value": _plain(config.sources[src][key])}
               for src in sorted(config.sources) for key in sorted(config.sources[src])]
    return {"raised": False, "preset": str(config.preset), "isozymes": bool(config.isozymes),
            "standard": bool(config.standard), "stages": stages, "sources": sources,
            "okp_methods": [str(m) for m in methods]}


# --------------------------------------------------------------------------- #
# collect_kcats / assign_kcats
# --------------------------------------------------------------------------- #

def _decisions(table: pd.DataFrame) -> list:
    rows = []
    for r in table.sort_values("ec_rxn").itertuples():
        rows.append({"ec_rxn": str(r.ec_rxn), "stage": str(r.stage), "source": str(r.source),
                     "kcat": _num(r.kcat), "candidates": int(r.candidates), "spread": _num(r.spread_log10),
                     "sd": _num(r.sd_log10)})
    return rows


def _collect_case(inputs: dict, root: Path, name: str, files: dict, toml: str, assign: bool = False) -> dict:
    project = _project(inputs, root, name, files=files, toml=toml)
    model, adapter = _model(project)
    raised, result = _attempt(lambda: collect_kcats(model, adapter=adapter, write_report=False))
    if raised:
        return {"raised": True, "decisions": [], "missing": [], "assigned": []}
    decisions, report, _ = result
    out = {"raised": False, "decisions": _decisions(decisions), "missing": sorted(report["missing"]),
           "assigned": []}
    if assign:
        raised, _ = _attempt(lambda: assign_kcats(model, adapter=adapter, write_report=False))
        out["assign_raised"] = raised
        out["assigned"] = [{"ec_rxn": str(r), "kcat": _num(k), "source": str(s)}
                           for r, k, s in sorted(zip(model.ec.rxns, model.ec.kcat, model.ec.source))]
    return out


# --------------------------------------------------------------------------- #
# EC tree
# --------------------------------------------------------------------------- #

QUERIES = ["1.1.1.1", "1.1.1.1;2.7.1.2", "", "9.9.9.9", "1.1.2.1;1.1.2.1", "3.5.1.n3", "x.y",
           "2.7.1.1;2.7.1.2", "1.1.1.1;1.1.1.2", "1.1.1.2;1.1.1.9", "2.7.1.2;4.1.1.1", " 1.1.1.1 ; 2.7.1.1 "]


def _estimates(tree_or_error, queries) -> list:
    out = []
    for q in queries:
        value, n, depth, ec = tree_or_error.estimate(q)
        out.append({"query": q, "log10": _num(value), "n": int(n), "depth": int(depth), "ec": str(ec)})
    return out


def _tree_case(inputs: dict) -> dict:
    table = load_brenda_measurements(Path(inputs["brenda_dir"]))
    out = {}
    for label, k in (("k5", 5.0), ("k1", 1.0), ("k0", 0.0), ("k_negative", -1.0), ("k_nan", float("nan"))):
        raised, tree = _attempt(lambda: ECTree(table, k=k))
        out[label] = {"raised": raised, "estimates": [] if raised else _estimates(tree, QUERIES)}
    nan_table = pd.DataFrame({"ec": ["1.1.1.1", "1.1.1.1", "2.7.1.1", "2.7.1.1", "2.7.1.1"],
                              "log10_kcat": [1.0, np.nan, 2.0, 3.0, 4.0]})
    raised, tree = _attempt(lambda: ECTree(nan_table, k=5.0))
    out["nan_in_table"] = {"raised": raised,
                           "estimates": [] if raised else _estimates(tree, ["1.1.1.1", "2.7.1.1", "9.9.9.9"])}
    return out


def _rows(frame: pd.DataFrame) -> list:
    rows = [{"rxn": str(r.rxn_id), "eccode": str(r.eccode), "kcat": _num(r.kcat), "detail": str(r.detail)}
            for r in frame.itertuples()]
    return sorted(rows, key=lambda r: (r["rxn"], r["kcat"]["value"]))


def _ec_tree_model_case(inputs: dict, root: Path) -> dict:
    project = _project(inputs, root, "ec_tree_model")
    model, adapter = _model(project)
    ids = sorted(model.ec.rxns)
    variants = {
        "all": {},
        "k1": {"k": 1.0},
        "subset": {"ec_rxns": ids[:3]},
        "subset_with_repeat": {"ec_rxns": [ids[0], ids[1], ids[1]]},
        "empty_subset": {"ec_rxns": []},
    }
    out = {"n_ec_rxns": len(ids)}
    for name, kwargs in variants.items():
        raised, frame = _attempt(lambda: kcats_from_ec_tree(model, save=False, adapter=adapter, **kwargs))
        out[name] = {"raised": raised, "rows": [] if raised else _rows(frame)}
    return out


def run(ctx):
    inputs = ctx["inputs"]
    root = Path(inputs["scratch"])
    select_files = {"custom": "custom", "brenda": "brenda", "dlkcat": "dlkcat"}
    combine_files = {"custom": "custom", "brenda": "brenda", "ec_tree": "ec_tree", "homologs": "homologs"}
    return {
        "config": {name: _config_case(inputs, root, name) for name in CONFIGS},
        "collect": {
            "select_complete": _collect_case(inputs, root, "select_complete", select_files, SELECT),
            "select_dlkcat_without_rows": _collect_case(
                inputs, root, "select_empty", {**select_files, "dlkcat": "empty_dlkcat"}, SELECT),
            "select_no_custom_file": _collect_case(
                inputs, root, "select_nocustom", {"brenda": "brenda", "dlkcat": "dlkcat"}, SELECT),
            "select_isozymes_and_standard": _collect_case(
                inputs, root, "select_fill", select_files, SELECT_FILL),
            "select_assigned": _collect_case(inputs, root, "select_assign", select_files, SELECT, assign=True),
            "combine_complete": _collect_case(inputs, root, "combine_complete", combine_files, COMBINE),
            "combine_homologs_without_rows": _collect_case(
                inputs, root, "combine_empty", {**combine_files, "homologs": "empty_homologs"}, COMBINE),
        },
        "ec_tree": _tree_case(inputs),
        "ec_tree_model": _ec_tree_model_case(inputs, root),
    }
