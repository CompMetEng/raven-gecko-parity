"""Python side of the relax_enz_concs scenario.

relax_enz_concs(per_complex=...) on ecTestGEM with a complex written into ec.rxn_enz_mat. A case that
raises is a result (`raised`), so the two sides can disagree about raising.
"""

import importlib.util
import math
import sys
from pathlib import Path

import cobra
import numpy as np

from geckopy import (
    ModelAdapter,
    apply_kcat_constraints,
    constrain_enz_concs,
    fill_enz_concs,
    fill_eccodes_from_gem,
    level_complex_concs,
    relax_enz_concs,
    load_conventional_gem,
    load_prot_data,
    make_ec_model,
    set_prot_pool_size,
)

PROTEINS = ["P1", "P2", "P3", "P4", "P5"]


def _adapter(inputs: dict) -> ModelAdapter:
    """geckopy's ecTestGEM adapter, repointed at the GECKO copy of the fixture."""
    base = ModelAdapter.from_folder(inputs["adapter_python"])
    fixture = Path(inputs["fixture_dir"])
    base.params.path = fixture
    base.params.conv_gem = fixture / "models" / "testModel.xml"

    adapter_module_path = Path(inputs["adapter_python"]) / "adapter.py"
    spec = importlib.util.spec_from_file_location("ectestgem_adapter", adapter_module_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module.TestGEMAdapter(base.params)


def num(x) -> dict:
    """A number as a record, so that NaN is a value both sides can state."""
    x = float(x)
    return {"is_nan": bool(math.isnan(x)), "value": 0.0 if math.isnan(x) else x}


def run_case(fn) -> dict:
    """Result of one case; a raise is a result on both sides, not a failure."""
    try:
        value = fn()
    except Exception:
        return {"raised": True, "value": {}}
    return {"raised": False, "value": value}


def build(inputs: dict):
    """The R2/R4-blocked single-route ecTestGEM of flexibilize_edge_cases_ectestgem, concentrations
    filled and not yet constrained."""
    adapter = _adapter(inputs)
    conv = load_conventional_gem(adapter)
    for rxn_id in ("R2", "R4"):
        rxn = conv.reactions.get_by_id(rxn_id)
        rxn.lower_bound = 0.0
        rxn.upper_bound = 0.0

    model = make_ec_model(conv, adapter, gecko_light=False)
    model.solver = str(inputs["python_solver"])
    fill_eccodes_from_gem(model)
    model.adapter = adapter

    index = {r: i for i, r in enumerate(model.ec.rxns)}
    for rxn_id, kcat in inputs["kcats"].items():
        i = index[str(rxn_id)]
        model.ec.kcat[i] = float(kcat)
        model.ec.source[i] = "manual"

    prot_data = load_prot_data(Path(adapter.params.path) / "data" / "proteomics.tsv", [1])
    fill_enz_concs(model, prot_data)
    return model


def finish(model) -> None:
    constrain_enz_concs(model)
    apply_kcat_constraints(model)
    set_prot_pool_size(model)


def set_layout(model, layout: dict) -> None:
    """Replace the rows of ec.rxn_enz_mat by the copies of the layout."""
    rxns = list(model.ec.rxns)
    enzymes = list(model.ec.enzymes)
    matrix = np.zeros(model.ec.rxn_enz_mat.shape)
    for rxn, subunits in layout.items():
        for enzyme, copies in (subunits or {}).items():
            matrix[rxns.index(str(rxn)), enzymes.index(str(enzyme))] = float(copies)
    model.ec.rxn_enz_mat = type(model.ec.rxn_enz_mat)(matrix)


def set_concs(model, **concs) -> None:
    enzymes = list(model.ec.enzymes)
    for protein, value in concs.items():
        model.ec.concs[enzymes.index(protein)] = float(value)


def concs_list(model) -> list:
    enzymes = list(model.ec.enzymes)
    return [{"protein": p, "conc": num(model.ec.concs[enzymes.index(p)])} for p in PROTEINS]


def usage_ub_list(model) -> list:
    return [
        {"protein": p, "ub": num(model.reactions.get_by_id("usage_prot_" + p).upper_bound)}
        for p in PROTEINS
    ]


def relax_case(inputs: dict, *, statistic, exp_growth: float, per_complex: bool) -> dict:
    """statistic None leaves the measured concentrations of the complex as they are."""
    model = build(inputs)
    set_layout(model, inputs["layouts"]["binding"])
    if statistic is not None:
        level_complex_concs(model, statistic=statistic)
    finish(model)
    start = float(model.optimize().objective_value)
    result = relax_enz_concs(model, exp_growth=exp_growth, bio_rxn="R5", per_complex=per_complex)
    ids = [str(p) for p in result.uniprot_ids]
    relaxed = sorted(
        (
            {
                "protein": ids[i],
                "old_conc": num(result.old_concs[i]),
                "flex_conc": num(result.flex_concs[i]),
                "ratio_incr": num(result.ratio_incr[i]),
            }
            for i in range(len(ids))
        ),
        key=lambda r: r["protein"],
    )
    final = model.optimize()
    optimal = final.status == "optimal"
    return {
        "concs": concs_list(model),
        "start_growth": num(start),
        "relaxed": relaxed,
        "relaxed_set": sorted(ids),
        "final_is_optimal": bool(optimal),
        "final_growth": num(final.objective_value if optimal else 0.0),
        "usage_ub": usage_ub_list(model),
        "pool_ub": num(model.reactions.get_by_id("prot_pool_exchange").upper_bound),
    }


def run(ctx):
    inputs = ctx["inputs"]
    cobra.Configuration().solver = str(inputs["python_solver"])
    cases = {}
    targets = {
        "measured": (None, float(inputs["exp_growth"])),
        "median": ("median", float(inputs["exp_growth"])),
        "high": (None, float(inputs["exp_growth_high"])),
    }
    for label, (statistic, target) in targets.items():
        for per_complex in (False, True):
            name = f"relax_{label}_{'per_complex' if per_complex else 'default'}"
            cases[name] = run_case(
                lambda s=statistic, t=target, p=per_complex: relax_case(
                    inputs, statistic=s, exp_growth=t, per_complex=p))
    return cases
