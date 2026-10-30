"""Python side of the ecModel edit edge-case scenario (see scenario.yml).

A call that raises is reported as `raised: true` with an empty result; which exception it was
is not the finding.
"""

import importlib.util
import sys
from pathlib import Path

import cobra
import pandas as pd

from geckopy import (
    ModelAdapter,
    add_new_rxns_to_ec,
    apply_kcat_constraints,
    ec_fva,
    fill_eccodes_from_gem,
    load_conventional_gem,
    make_ec_model,
    map_rxns_to_conv,
    set_prot_pool_size,
)
from geckopy.utilities.add_new_rxns_to_ec import NewEnzyme


def _bound(value: float) -> tuple[float, str]:
    value = float(value)
    if value == float("inf"):
        return 0.0, "+inf"
    if value == float("-inf"):
        return 0.0, "-inf"
    return value, "finite"


def _adapter(inputs: dict) -> ModelAdapter:
    """geckopy's ecTestGEM adapter, repointed at the GECKO copy of the fixture, as an adapter
    object with the TestGEMAdapter methods (ec_fva_ectestgem builds it the same way)."""
    base = ModelAdapter.from_folder(inputs["adapter_python"])
    fixture = Path(inputs["fixture_dir"])
    base.params.path = fixture
    base.params.conv_gem = fixture / "models" / "testModel.xml"
    spec = importlib.util.spec_from_file_location("ectestgem_adapter", Path(inputs["adapter_python"]) / "adapter.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module.TestGEMAdapter(base.params)


def _summary(model) -> dict:
    reactions = []
    for rxn in sorted(model.reactions, key=lambda r: r.id):
        lower, lower_kind = _bound(rxn.lower_bound)
        upper, upper_kind = _bound(rxn.upper_bound)
        reactions.append({
            "id": rxn.id, "lower_bound": lower, "lower_kind": lower_kind,
            "upper_bound": upper, "upper_kind": upper_kind,
        })
    associations = sorted(
        ({"reaction": rxn.id, "gene": gene.id} for rxn in model.reactions for gene in rxn.genes),
        key=lambda e: (e["reaction"], e["gene"]),
    )
    ec = model.ec
    matrix = ec.rxn_enz_mat
    if hasattr(matrix, "toarray"):
        matrix = matrix.toarray()
    coupling = sorted(
        ({"reaction": ec.rxns[i], "enzyme": ec.enzymes[j], "coefficient": float(matrix[i][j])}
         for i in range(len(ec.rxns)) for j in range(len(ec.enzymes)) if matrix[i][j] != 0),
        key=lambda e: (e["reaction"], e["enzyme"]),
    )
    enzymes = sorted(
        ({"enzyme": str(e), "gene": str(g), "mw": float(m)} for e, g, m in zip(ec.enzymes, ec.genes, ec.mw)),
        key=lambda e: e["enzyme"],
    )
    return {
        "reactions": reactions,
        "metabolites": sorted(met.id for met in model.metabolites),
        "gene_associations": associations,
        "ec_rxns": sorted(str(r) for r in ec.rxns),
        "ec_enzymes": enzymes,
        "ec_coupling": coupling,
    }


def _parse_equation(equation: str) -> tuple[list[str], list[str], bool]:
    reversible = "<=>" in equation
    left, right = (side.strip() for side in equation.split("<=>" if reversible else "=>"))
    return [m.strip() for m in left.split("+")], [m.strip() for m in right.split("+")], reversible


def _cobra_id(token: str) -> str:
    name, compartment = token.rstrip("]").split("[")
    return f"{name}{compartment}"


def _add_case(inputs: dict, case: dict) -> dict:
    adapter = _adapter(inputs)
    model = make_ec_model(load_conventional_gem(adapter), adapter, gecko_light=False)
    spec = inputs["new_reaction"]
    reactants, products, reversible = _parse_equation(spec["equation"])
    new_rxn = cobra.Reaction(id=spec["id"], name=spec["name"])
    new_rxn.lower_bound = -1000.0 if reversible else 0.0
    new_rxn.upper_bound = 1000.0
    for token in reactants:
        new_rxn.add_metabolites({model.metabolites.get_by_id(_cobra_id(token)): -1.0})
    for token in products:
        new_rxn.add_metabolites({model.metabolites.get_by_id(_cobra_id(token)): 1.0})
    new_rxn.gene_reaction_rule = str(case["gr_rule"])
    enzymes = [NewEnzyme(enzyme=e["enzyme"], gene=e["gene"], mw=float(e["mw"])) for e in case["enzymes"]]
    try:
        result = add_new_rxns_to_ec(model, [new_rxn], enzymes)
    except Exception:  # noqa: BLE001 --- whether it raised is the finding, not the exception type
        return {"raised": True, "rxns_added": [], "enz_added": [], "model": {}}
    return {
        "raised": False,
        "rxns_added": sorted(result.rxns_added),
        "enz_added": sorted(result.enz_added),
        "model": _summary(model),
    }


def _renamed(inputs: dict, case: dict):
    adapter = _adapter(inputs)
    conv = load_conventional_gem(adapter).copy()
    conv.reactions.get_by_id(str(case["from"])).id = str(case["to"])
    return adapter, conv


def _map_rxns(adapter, conv) -> dict:
    model = make_ec_model(conv, adapter, gecko_light=False)
    ranked = {rid: float(rank) for rank, rid in enumerate(sorted(r.id for r in model.reactions), start=1)}
    try:
        result = map_rxns_to_conv(model, conv, pd.Series(ranked))
    except Exception:  # noqa: BLE001
        return {"raised": True, "mapped": [], "usage": []}
    mapped = sorted(
        ({"reaction": rxn.id, "flux": float(v)} for rxn, v in zip(conv.reactions, result.mapped_flux)),
        key=lambda e: e["reaction"],
    )
    usage = sorted(
        ({"enzyme": label, "flux": float(v)} for label, v in zip(result.usage_enz, result.enz_usage_flux)),
        key=lambda e: e["enzyme"],
    )
    return {"raised": False, "mapped": mapped, "usage": usage}


def _ec_fva(inputs: dict, adapter, conv) -> dict:
    model = make_ec_model(conv, adapter, gecko_light=False)
    model.solver = str(inputs["python_solver"])
    fill_eccodes_from_gem(model)
    index = {r: i for i, r in enumerate(model.ec.rxns)}
    for rxn_id, kcat in inputs["kcats"].items():
        model.ec.kcat[index[str(rxn_id)]] = float(kcat)
        model.ec.source[index[str(rxn_id)]] = "manual"
    apply_kcat_constraints(model)
    set_prot_pool_size(model)
    try:
        df = ec_fva(model, conv, progress=False, n_proc=int(inputs["python_n_proc"]))
    except Exception:  # noqa: BLE001
        return {"raised": True, "fva": []}
    rows = sorted(
        ({"reaction": str(rid), "min_flux": float(row.min_flux), "max_flux": float(row.max_flux)}
         for rid, row in df.iterrows()),
        key=lambda r: r["reaction"],
    )
    return {"raised": False, "fva": rows}


def run(ctx):
    inputs = ctx["inputs"]
    cobra.Configuration().solver = str(inputs["python_solver"])
    ids = {}
    for case in inputs["id_cases"]:
        adapter, conv = _renamed(inputs, case)
        ids[str(case["name"])] = {
            "map_rxns": _map_rxns(adapter, conv),
            "ec_fva": _ec_fva(inputs, adapter, conv),
        }
    return {
        "add_new_rxns": {str(c["name"]): _add_case(inputs, c) for c in inputs["add_cases"]},
        "ids_with_suffix_text": ids,
    }
