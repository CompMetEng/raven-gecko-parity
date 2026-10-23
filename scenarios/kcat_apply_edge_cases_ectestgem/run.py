"""Python side of the kcat-apply edge-case scenario.

Every case builds its own state (a fresh full or light ecTestGEM model, or a small
hand-built light model), runs one call and records what the call left behind: the kcat
rows, the enzyme coefficients of the S matrix, and whether the call raised. A raise is a
value (`raised`), not a failure, so a case where one side raises and the other returns
shows up as a difference in that value.
"""

import importlib.util
import sys
from pathlib import Path

import numpy as np
from scipy import sparse

from geckopy import (
    EcData,
    ModelAdapter,
    apply_custom_kcats,
    apply_kcat_constraints,
    fill_kcats_from_isozymes,
    load_conventional_gem,
    make_ec_model,
    set_kcat_for_reactions,
)


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


def _model(inputs: dict, light: bool):
    adapter = _adapter(inputs)
    conv = load_conventional_gem(adapter)
    return make_ec_model(conv, adapter, gecko_light=light)


def _prot_coefficients(model) -> list[dict]:
    # Coefficients of the metabolic reactions only: the usage reactions and the pool
    # exchange carry fixed -1 / +1 entries that do not depend on any kcat.
    records = []
    for rxn in model.reactions:
        if rxn.id.startswith("usage_prot_") or rxn.id == "prot_pool_exchange":
            continue
        for met, coeff in rxn.metabolites.items():
            if met.id.startswith("prot_") and coeff != 0:
                records.append({"reaction": rxn.id, "metabolite": met.id, "coefficient": float(coeff)})
    return sorted(records, key=lambda e: (e["reaction"], e["metabolite"]))


def _num(value: float) -> dict:
    value = float(value)
    if np.isnan(value):
        return {"nan": True, "value": 0.0}
    return {"nan": False, "value": value}


def _rows(model, ids=None) -> list[dict]:
    ec = model.ec
    records = []
    for i, rxn in enumerate(ec.rxns):
        if ids is not None and rxn not in ids:
            continue
        records.append({
            "reaction": str(rxn),
            "kcat": _num(ec.kcat[i]),
            "source": str(ec.source[i]) if ec.source[i] else "",
            "notes": str(ec.notes[i]) if ec.notes[i] else "",
        })
    return sorted(records, key=lambda r: r["reaction"])


def _set_kcats(model, values: dict) -> None:
    index = {r: i for i, r in enumerate(model.ec.rxns)}
    for rxn, value in values.items():
        model.ec.kcat[index[rxn]] = value


def _attempt(fn) -> bool:
    try:
        fn()
        return False
    except Exception:  # noqa: BLE001 --- whether the call raises is the finding, not the type
        return True


# --------------------------------------------------------------------------- #
# hand-built light models
# --------------------------------------------------------------------------- #

def _light_model(adapter_inputs, kcat, enz_rows, mw):
    import cobra

    model = cobra.Model("light_edge")
    model.add_metabolites([cobra.Metabolite("prot_pool", compartment="c")])
    rxn = cobra.Reaction("R1")
    rxn.lower_bound = 0.0
    rxn.upper_bound = 1000.0
    model.add_reactions([rxn])
    n = len(kcat)
    model.ec = EcData(
        gecko_light=True,
        rxns=[f"{i + 1:03d}_R1" for i in range(n)],
        kcat=np.array(kcat, dtype=float),
        source=[""] * n,
        notes=[""] * n,
        eccodes=[""] * n,
        rxn_enz_mat=sparse.csr_matrix(np.array(enz_rows, dtype=float)),
        genes=[f"G{i + 1}" for i in range(len(mw))],
        enzymes=[f"P{i + 1}" for i in range(len(mw))],
        mw=np.array(mw, dtype=float),
        sequence=[""] * len(mw),
    )
    return model


def _case_light_subset(inputs) -> dict:
    out = {}
    # R1 has two isozymes; P1 costs 10000/200 = 50, P2 costs 5000/10 = 500. After the
    # second isozyme is tuned to a worse kcat, the caller re-applies the reaction it
    # changed. MATLAB names the cobra reaction (all isozymes compete); Python names the
    # ec rows it changed.
    for label, py_ids in (("changed_isozyme_only", ["002_R1"]), ("both_isozymes", ["001_R1", "002_R1"])):
        model = _light_model(inputs, [200.0, 10.0], [[1, 0], [0, 1]], [10000.0, 5000.0])
        apply_kcat_constraints(model)
        model.ec.kcat[1] = 20.0
        apply_kcat_constraints(model, update_rxns=py_ids)
        out[label] = _prot_coefficients(model)
    return out


def _case_light_empty_row(inputs) -> dict:
    # 002_R1 has a kcat but no enzyme in the coupling matrix.
    model = _light_model(inputs, [50.0, 30.0], [[1, 0], [0, 0]], [10000.0, 5000.0])
    apply_kcat_constraints(model)
    return {"coefficients": _prot_coefficients(model)}


def _case_light_ids(inputs) -> dict:
    out = {}
    for label, ids in (("cobra_id", ["R1"]), ("ec_id", ["001_R1"])):
        model = _light_model(inputs, [200.0, 10.0], [[1, 0], [0, 1]], [10000.0, 5000.0])
        out[label] = {"raised": _attempt(lambda: apply_kcat_constraints(model, update_rxns=ids))}
    return out


# --------------------------------------------------------------------------- #
# applyKcatConstraints on the full model
# --------------------------------------------------------------------------- #

def _all_kcat(model, value):
    model.ec.kcat[:] = value


def _case_full_zero_reapply(inputs) -> dict:
    model = _model(inputs, light=False)
    _all_kcat(model, 10.0)
    apply_kcat_constraints(model)
    first = _prot_coefficients(model)
    _all_kcat(model, 0.0)
    apply_kcat_constraints(model)
    return {"after_first_apply": first, "after_zero_reapply": _prot_coefficients(model)}


def _case_full_update_list(inputs) -> dict:
    out = {}
    model = _model(inputs, light=False)
    _all_kcat(model, 10.0)
    apply_kcat_constraints(model)
    _all_kcat(model, 20.0)
    out["empty_list"] = {
        "raised": _attempt(lambda: apply_kcat_constraints(model, update_rxns=[])),
        "coefficients": _prot_coefficients(model),
    }
    model = _model(inputs, light=False)
    _all_kcat(model, 10.0)
    out["unknown_id"] = {
        "raised": _attempt(lambda: apply_kcat_constraints(model, update_rxns=["NOSUCH"])),
        "coefficients": _prot_coefficients(model),
    }
    return out


# --------------------------------------------------------------------------- #
# setKcatForReactions
# --------------------------------------------------------------------------- #

def _case_set_kcat(inputs) -> dict:
    out = {}
    for label, rxn_ids in (("default_apply", ["R3"]), ("lowercase_id", ["r3"]), ("unknown_id", ["NOSUCH"])):
        model = _model(inputs, light=False)
        raised = _attempt(lambda: set_kcat_for_reactions(model, rxn_ids, 100.0))
        out[label] = {"raised": raised, "rows": _rows(model), "coefficients": _prot_coefficients(model)}
    model = _model(inputs, light=True)
    raised = _attempt(lambda: set_kcat_for_reactions(model, ["R3"], 100.0))
    out["light_base_id"] = {"raised": raised, "rows": _rows(model)}
    return out


# --------------------------------------------------------------------------- #
# getKcatAcrossIsozymes
# --------------------------------------------------------------------------- #

def _case_isozymes(inputs) -> dict:
    out = {}
    model = _model(inputs, light=False)
    nan = float("nan")
    _set_kcats(model, {"R2_EXP_1": nan, "R2_EXP_2": 5.0, "R2_REV_EXP_1": 0.0, "R2_REV_EXP_2": nan, "R3": 8.0})
    raised = _attempt(lambda: fill_kcats_from_isozymes(model, apply=False))
    out["nan_siblings"] = {"raised": raised, "rows": _rows(model)}

    model = _model(inputs, light=False)
    _set_kcats(model, {"R2_EXP_1": 4.0, "R2_EXP_2": 0.0, "R3": 8.0})
    raised = _attempt(lambda: fill_kcats_from_isozymes(model))   # default apply
    out["default_apply"] = {"raised": raised, "rows": _rows(model), "coefficients": _prot_coefficients(model)}
    return out


# --------------------------------------------------------------------------- #
# applyCustomKcats
# --------------------------------------------------------------------------- #

CUSTOM_FILES = ["ck_empty_notes", "ck_blank_kcat", "ck_lowercase", "ck_protein_only",
                "ck_quoted_list", "ck_partial_complex"]


def _case_custom(inputs) -> dict:
    out = {}
    folder = Path(inputs["custom_dir"])
    for name in CUSTOM_FILES:
        model = _model(inputs, light=False)
        raised = _attempt(lambda: apply_custom_kcats(model, folder / f"{name}.tsv", apply=False))
        out[name] = {"raised": raised, "rows": _rows(model)}
    for label, name in (("light_model", "ck_empty_notes"), ("light_protein_only", "ck_protein_only")):
        model = _model(inputs, light=True)
        raised = _attempt(lambda: apply_custom_kcats(model, folder / f"{name}.tsv", apply=False))
        out[label] = {"raised": raised, "rows": _rows(model)}
    return out


def run(ctx):
    inputs = ctx["inputs"]
    return {
        "light_subset": _case_light_subset(inputs),
        "light_empty_isozyme_row": _case_light_empty_row(inputs),
        "light_update_ids": _case_light_ids(inputs),
        "full_zero_reapply": _case_full_zero_reapply(inputs),
        "full_update_list": _case_full_update_list(inputs),
        "set_kcat": _case_set_kcat(inputs),
        "isozymes": _case_isozymes(inputs),
        "custom": _case_custom(inputs),
    }
