"""Python side of the kcat I/O edge-case scenario.

Every case builds its own state, runs one call and records what it left behind and whether
it raised. A raise is a value (`raised`). No case reaches PubChem: every metabolite name
used is in the vendored smilesDB.tsv of its case.
"""

import importlib.util
import sys
from pathlib import Path

import cobra
import numpy as np
import pandas as pd

from geckopy import (
    DLKcatIgnoreLists,
    ModelAdapter,
    assign_standard_kcat,
    find_met_smiles,
    load_conventional_gem,
    load_kcat_source,
    load_uniprot_tsv,
    make_ec_model,
    save_kcat_source,
)
from geckopy.gather_kcats.open_kinetics_predictor import build_okp_input_csv


def _adapter(inputs: dict) -> ModelAdapter:
    base = ModelAdapter.from_folder(inputs["adapter_python"])
    fixture = Path(inputs["fixture_dir"])
    base.params.path = fixture
    base.params.conv_gem = fixture / "models" / "testModel.xml"
    module_path = Path(inputs["adapter_python"]) / "adapter.py"
    spec = importlib.util.spec_from_file_location("ectestgem_adapter", module_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module.TestGEMAdapter(base.params)


def _model(inputs: dict, light: bool = False, extra: bool = False):
    adapter = _adapter(inputs)
    conv = load_conventional_gem(adapter)
    if extra:
        _add_extra_reactions(conv)
    return make_ec_model(conv, adapter, gecko_light=light)


def _add_extra_reactions(conv: cobra.Model) -> None:
    """Three reactions without a GPR: an ordinary one, a sink with two substrates and one
    named like a biomass pseudoreaction."""
    m1 = conv.metabolites.get_by_id("m1c")
    m2 = conv.metabolites.get_by_id("m2c")
    specs = [("R6", "R6", {m1: -1.0, m2: 1.0}),
             ("R7", "R7", {m1: -1.0, m2: -1.0}),
             ("R8", "Biomass Pseudoreaction", {m1: -1.0, m2: 1.0})]
    for rxn_id, name, stoich in specs:
        rxn = cobra.Reaction(rxn_id, name=name)
        rxn.lower_bound = 0.0
        rxn.upper_bound = 1000.0
        rxn.add_metabolites(stoich)
        conv.add_reactions([rxn])


def _attempt(fn):
    try:
        return False, fn()
    except Exception:  # noqa: BLE001 --- whether the call raises is the finding, not the type
        return True, None


def _num(value) -> dict:
    value = float(value)
    if np.isnan(value):
        return {"nan": True, "value": 0.0}
    return {"nan": False, "value": value}


# --------------------------------------------------------------------------- #
# findMetSmiles
# --------------------------------------------------------------------------- #

def _smiles_snapshot(model) -> dict:
    return {m.id: str(m.annotation.get("smiles", "")) for m in sorted(model.metabolites, key=lambda m: m.id)}


def _case_smiles(inputs) -> dict:
    out = {}
    model = _model(inputs)
    raised, _ = _attempt(lambda: find_met_smiles(
        model, cache_path=Path(inputs["dup_root"]) / "data" / "smilesDB.tsv"))
    out["duplicate_cache_name"] = {"raised": raised, "smiles": _smiles_snapshot(model)}

    model = _model(inputs)
    model.metabolites.get_by_id("m1c").name = ""
    raised, _ = _attempt(lambda: find_met_smiles(
        model, cache_path=Path(inputs["empty_root"]) / "data" / "smilesDB.tsv"))
    out["empty_metabolite_name"] = {"raised": raised, "smiles": _smiles_snapshot(model)}
    return out


# --------------------------------------------------------------------------- #
# OpenKineticsPredictor input
# --------------------------------------------------------------------------- #

def _case_okp(inputs) -> dict:
    model = _model(inputs)
    find_met_smiles(model, cache_path=Path(inputs["okp_root"]) / "data" / "smilesDB.tsv")
    ignore = DLKcatIgnoreLists(ignore_names=[], ignore_smiles=[], currency_pairs=[])
    raised, text = _attempt(lambda: build_okp_input_csv(
        model, ignore, ec_rxns=list(model.ec.rxns), only_with_smiles=True))
    lines = [] if raised else sorted(line for line in text.strip("\n").split("\n")[1:] if line)
    return {"all_reactions": {"raised": raised, "lines": lines}}


# --------------------------------------------------------------------------- #
# kcat source files
# --------------------------------------------------------------------------- #

def _case_save_source(inputs) -> dict:
    model = _model(inputs)
    kcat_list = pd.DataFrame({
        "rxn_id": ["R3", "R5"],
        "kcat": [25.0, 30.0],
        "substrates": ["m1", "m2"],
        "eccode": ["1.1.2.1", ""],
        "source": ["x", "x"],
        "detail": ["a\tb", "; x"],
    })
    folder = Path(inputs["save_dir"])
    raised, path = _attempt(lambda: save_kcat_source(model, kcat_list, "edge", folder=folder))
    lines = []
    if not raised:
        for line in Path(path).read_text(encoding="utf-8").splitlines():
            if not line.startswith(("# generator:", "# created:")):
                lines.append(line)
    return {"raised": raised, "lines": lines}


def _case_load_source(inputs) -> dict:
    model = _model(inputs)
    unclosed_raised, _ = _attempt(lambda: load_kcat_source(model, Path(inputs["odd_quote_source"])))
    raised, result = _attempt(lambda: load_kcat_source(model, Path(inputs["odd_source"])))
    if raised:
        return {"raised": True, "rows": [], "n_invalid": 0, "n_unmatched": 0, "n_partial": 0,
                "unclosed_quote_raised": unclosed_raised}
    df, report = result
    rows = sorted(
        ({"reaction": str(r.rxn_id), "kcat": float(r.kcat)} for r in df.itertuples()),
        key=lambda r: r["reaction"],
    )
    return {
        "raised": False,
        "rows": rows,
        "n_invalid": len(report["invalid"]),
        "n_unmatched": len(report["unmatched"]),
        "n_partial": len(report["partial"]),
        "unclosed_quote_raised": unclosed_raised,
    }


# --------------------------------------------------------------------------- #
# getStandardKcat
# --------------------------------------------------------------------------- #

def _rows(model, ids) -> list[dict]:
    ec = model.ec
    index = {r: i for i, r in enumerate(ec.rxns)}
    records = []
    for r in ids:
        if r not in index:
            continue
        i = index[r]
        records.append({"reaction": r, "kcat": _num(ec.kcat[i]),
                        "source": str(ec.source[i]) if ec.source[i] else ""})
    return sorted(records, key=lambda e: e["reaction"])


def _case_standard(inputs) -> dict:
    out = {}
    uniprot = load_uniprot_tsv(inputs["uniprot_path"], id_type="taxonomy_id", auto_download=False)

    model = _model(inputs, extra=True)
    model.ec.kcat[:] = 10.0
    before = set(model.ec.rxns)
    raised, _ = _attempt(lambda: assign_standard_kcat(model, uniprot, fill_zero_kcat=True))
    after = [str(r) for r in model.ec.rxns]
    out["extra_reactions"] = {
        "raised": raised,
        "new_ec_rxns": sorted(set(after) - before),
        "rows": _rows(model, ["R6", "R7", "R8"]),
    }

    model = _model(inputs, light=True)
    model.ec.kcat[:] = 10.0
    steps = {}
    for label in ("first_call", "second_call"):
        raised, _ = _attempt(lambda: assign_standard_kcat(model, uniprot, fill_zero_kcat=True))
        steps[label] = {
            "raised": raised,
            "enzymes": sorted(str(e) for e in model.ec.enzymes),
            "ec_rxns": sorted(str(r) for r in model.ec.rxns),
        }
    out["light_called_twice"] = steps
    return out


def run(ctx):
    inputs = ctx["inputs"]
    return {
        "smiles": _case_smiles(inputs),
        "okp_input": _case_okp(inputs),
        "save_source": _case_save_source(inputs),
        "load_source": _case_load_source(inputs),
        "standard": _case_standard(inputs),
    }
