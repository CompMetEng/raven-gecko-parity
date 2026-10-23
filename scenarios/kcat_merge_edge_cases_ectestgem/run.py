"""Python side of the kcat merge edge-case scenario.

Every case builds its own state, runs one call and records what came back and whether the
call raised. A raise is a value (`raised`), so a case where one side raises and the other
returns shows up as a difference in that value. Lists whose order is the result (the
written DLKcat input, the fuzzy matches of a subset) are emitted in the order returned.
"""

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from geckopy import (
    ModelAdapter,
    apply_kcat_list,
    fill_eccodes_from_gem,
    fuzzy_kcat_matching,
    load_brenda_data,
    load_conventional_gem,
    load_dlkcat_ignore_lists,
    load_phyl_dist,
    make_ec_model,
    merge_dlkcat_and_fuzzy_kcats,
    merge_kcats,
    read_dlkcat_output,
    write_dlkcat_input,
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


def _model(inputs: dict):
    adapter = _adapter(inputs)
    model = make_ec_model(load_conventional_gem(adapter), adapter, gecko_light=False)
    fill_eccodes_from_gem(model)
    return model


def _num(value) -> dict:
    value = float(value)
    if np.isnan(value):
        return {"nan": True, "value": 0.0}
    return {"nan": False, "value": value}


def _attempt(fn):
    try:
        return False, fn()
    except Exception:  # noqa: BLE001 --- whether the call raises is the finding, not the type
        return True, None


def _index(model, rxn_id: str) -> int:
    return list(model.ec.rxns).index(rxn_id)


# --------------------------------------------------------------------------- #
# merging kcat lists
# --------------------------------------------------------------------------- #

def _kcat_list(rows: list[dict]) -> pd.DataFrame:
    n = len(rows)
    return pd.DataFrame({
        "rxn_id": [r["rxn"] for r in rows],
        "source": [r["source"] for r in rows],
        "eccode": [""] * n,
        "substrates": [[""] for _ in rows],
        "genes": [[""] for _ in rows],
        "kcat": [float(r["kcat"]) for r in rows],
        "wildcard_level": pd.array([r.get("wc", pd.NA) for r in rows], dtype="Int64"),
        "origin": pd.array([r.get("origin", pd.NA) for r in rows], dtype="Int64"),
    })


def _merged(df: pd.DataFrame) -> list[dict]:
    records = [
        {"reaction": str(r.rxn_id), "kcat": float(r.kcat), "source": str(r.source)}
        for r in df.itertuples()
    ]
    return sorted(records, key=lambda r: (r["reaction"], r["source"]))


def _case_origin_limit() -> dict:
    # A BRENDA match at wildcard level 0 competes with a DLKcat prediction for the same
    # reaction; its origin is 5, 6, 4 and 3. The default limit of origins that count as the
    # best database tier decides which of the two wins.
    fuzzy = _kcat_list([
        {"rxn": "R1", "source": "brenda", "kcat": 5, "wc": 0, "origin": 5},
        {"rxn": "R2", "source": "brenda", "kcat": 6, "wc": 0, "origin": 6},
        {"rxn": "R3", "source": "brenda", "kcat": 4, "wc": 0, "origin": 4},
        {"rxn": "R5", "source": "brenda", "kcat": 3, "wc": 0, "origin": 3},
    ])
    dlkcat = _kcat_list([
        {"rxn": r, "source": "DLKcat", "kcat": k} for r, k in (("R1", 7), ("R2", 8), ("R3", 9), ("R5", 10))
    ])
    priority = ["database_top", "dlkcat", "database_bottom"]
    return {
        "merge_kcats": _merged(merge_kcats(fuzzy, dlkcat, source_priority=priority)),
        "merge_dlkcat_and_fuzzy": _merged(merge_dlkcat_and_fuzzy_kcats(dlkcat, fuzzy)),
    }


def _case_okp_prefix() -> dict:
    okp = _kcat_list([
        {"rxn": "R1", "source": "OKP-BRENDA", "kcat": 3},
        {"rxn": "R2", "source": "OKP-CataPro", "kcat": 4},
    ])
    dlkcat = _kcat_list([
        {"rxn": "R1", "source": "DLKcat", "kcat": 9},
        {"rxn": "R2", "source": "DLKcat", "kcat": 9},
    ])
    priority = ["database_exact", "catapro", "dlkcat"]
    return {"merged": _merged(merge_kcats(okp, dlkcat, source_priority=priority))}


# --------------------------------------------------------------------------- #
# selecting one kcat per reaction
# --------------------------------------------------------------------------- #

def _case_assign_zero_row(inputs) -> dict:
    model = _model(inputs)
    kcat_list = pd.DataFrame({"rxn_id": ["R3", "R5"], "kcat": [0.0, 5.0], "source": ["a", "b"]})
    raised, _ = _attempt(lambda: apply_kcat_list(model, kcat_list, criteria="max"))
    rows = []
    for rxn in ("R3", "R5"):
        i = _index(model, rxn)
        rows.append({"reaction": rxn, "kcat": _num(model.ec.kcat[i]),
                     "source": str(model.ec.source[i]) if model.ec.source[i] else ""})
    return {"raised": raised, "rows": rows}


# --------------------------------------------------------------------------- #
# DLKcat input and output
# --------------------------------------------------------------------------- #

DLKCAT_FILES = ["dl_unknown_substrate", "dl_unknown_rxn_na", "dl_inf", "dl_dotted", "dl_all_none"]


def _dlkcat_rows(df: pd.DataFrame) -> list[dict]:
    records = [
        {
            "reaction": str(r.rxn_id),
            "gene": str(r.genes[0]) if len(r.genes) else "",
            "substrate": str(r.substrates[0]) if len(r.substrates) else "",
            "kcat": _num(r.kcat),
        }
        for r in df.itertuples()
    ]
    return sorted(records, key=lambda r: (r["reaction"], r["gene"]))


def _case_read_dlkcat(inputs) -> dict:
    out = {}
    folder = Path(inputs["data_dir"])
    for name in DLKCAT_FILES:
        model = _model(inputs)
        raised, df = _attempt(lambda: read_dlkcat_output(model, folder / f"{name}.tsv"))
        out[name] = {"raised": raised, "rows": [] if raised else _dlkcat_rows(df)}
    return out


def _case_write_dlkcat(inputs) -> dict:
    out = {}
    ignore_lists = load_dlkcat_ignore_lists()
    target = Path(inputs["dlkcat_write_target"])
    for label, ids in (("empty_selection", []), ("selection_order", ["R5", "R3"])):
        model = _model(inputs)
        raised, df = _attempt(lambda: write_dlkcat_input(
            model, target, ignore_lists, ec_rxns=ids, overwrite=True, only_with_smiles=False))
        rows = [] if raised else [
            {"reaction": str(r.rxn_id), "gene": str(r.gene), "substrate": str(r.substrate)}
            for r in df.itertuples()
        ]
        out[label] = {"raised": raised, "rows": rows}
    return out


# --------------------------------------------------------------------------- #
# BRENDA loading and fuzzy matching
# --------------------------------------------------------------------------- #

def _case_outlier(inputs) -> dict:
    brenda = load_brenda_data(inputs["outlier_brenda_dir"])
    rows = [
        {"eccode": str(r.ec_code), "substrate": str(r.substrate), "organism": str(r.organism),
         "kcat": float(r.kcat)}
        for r in brenda.kcat_max.itertuples()
    ]
    rows.sort(key=lambda r: (r["eccode"], r["substrate"], r["organism"]))
    return {"kcat_max": rows}


def _match_rows(df: pd.DataFrame) -> list[dict]:
    return [
        {
            "reaction": str(r.rxn_id),
            "kcat": float(r.kcat),
            "eccode": str(r.eccode) if r.eccode else "",
            "origin": -1 if pd.isna(r.origin) else int(r.origin),
            "wildcard_level": -1 if pd.isna(r.wildcard_level) else int(r.wildcard_level),
        }
        for r in df.itertuples()
    ]


def _case_fuzzy(inputs) -> dict:
    out = {}
    brenda = load_brenda_data(inputs["brenda_dir"])
    phyl = load_phyl_dist(inputs["phyl_dist_path"])

    def run(model, ids):
        return fuzzy_kcat_matching(model, brenda, phyl, ec_rxns=ids)

    model = _model(inputs)
    model.ec.eccodes[_index(model, "R3")] = "9.9.9.9"
    raised, df = _attempt(lambda: run(model, ["R3"]))
    out["full_wildcard"] = {"raised": raised, "matches": [] if raised else _match_rows(df)}

    model = _model(inputs)
    model.ec.eccodes[_index(model, "R2_EXP_1")] = "1.1.1.5;1.1.1.1"
    raised, df = _attempt(lambda: run(model, ["R2_EXP_1"]))
    out["tie_between_ec_codes"] = {"raised": raised, "matches": [] if raised else _match_rows(df)}

    model = _model(inputs)
    raised, df = _attempt(lambda: run(model, ["R3", "R2_EXP_1", "R2_EXP_1"]))
    out["subset_order_and_repeat"] = {"raised": raised, "matches": [] if raised else _match_rows(df)}

    model = _model(inputs)
    raised, df = _attempt(lambda: run(model, []))
    out["empty_subset"] = {"raised": raised, "n_matches": 0 if raised else len(df)}
    return out


def run(ctx):
    inputs = ctx["inputs"]
    return {
        "merge_origin_limit": _case_origin_limit(),
        "merge_okp_prefix": _case_okp_prefix(),
        "assign_zero_row": _case_assign_zero_row(inputs),
        "read_dlkcat": _case_read_dlkcat(inputs),
        "write_dlkcat": _case_write_dlkcat(inputs),
        "brenda_outlier": _case_outlier(inputs),
        "fuzzy": _case_fuzzy(inputs),
    }
