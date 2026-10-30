"""Python side of the flux-data edge-case scenario.

loadFluxData and constrainFluxData on files and arguments the happy-path scenario does not
reach: a ragged row, a file with only the three fixed columns, a header with nested parentheses,
a missing file, a condition given by name and one that is not in the file, an exchange reaction
that is not in the model and an invalid maxMinGrowth. A case that raises is a result
(`raised`), so the two sides can disagree about raising; load_control and constrain_control
return real content on both sides.
"""

from pathlib import Path

import cobra

from geckopy import apply_flux_data_constraints
from geckopy.databases.flux_data import load_flux_data
import importlib.util
import math
import sys
from pathlib import Path

import cobra

from geckopy import (
    ModelAdapter,
    apply_kcat_constraints,
    constrain_enz_concs,
    fill_enz_concs,
    fill_eccodes_from_gem,
    load_conventional_gem,
    load_prot_data,
    make_ec_model,
    set_prot_pool_size,
)


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


def numlist(values) -> list:
    return [num(v) for v in values]


def run_case(fn) -> dict:
    """Result of one case; a raise is a result on both sides, not a failure."""
    try:
        value = fn()
    except Exception:
        return {"raised": True, "value": {}}
    return {"raised": False, "value": value}


def build(inputs: dict, *, constrain: bool = True, kcats: dict | None = None):
    """The R2/R4-blocked single-route ecTestGEM of enzyme_usage_ectestgem."""
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
    for rxn_id, kcat in (kcats or inputs["kcats"]).items():
        i = index[str(rxn_id)]
        model.ec.kcat[i] = float(kcat)
        model.ec.source[i] = "manual"

    prot_data = load_prot_data(Path(adapter.params.path) / "data" / "proteomics.tsv", [1])
    fill_enz_concs(model, prot_data)
    if constrain:
        constrain_enz_concs(model)
    apply_kcat_constraints(model)
    set_prot_pool_size(model)
    return model


def enzyme_index(model, protein: str) -> int:
    return list(model.ec.enzymes).index(protein)


def usage_flux(model, solution, protein: str) -> float:
    return abs(float(solution.fluxes["usage_prot_" + protein]))


def tie_p4(model) -> None:
    """Make P4's cap equal to the flux it carries at the optimum, so that P4 and P5 both bind."""
    sol = model.optimize()
    model.ec.concs[enzyme_index(model, "P4")] = usage_flux(model, sol, "P4")
    constrain_enz_concs(model)


def bounds(model, rxn_id: str) -> dict:
    rxn = model.reactions.get_by_id(rxn_id)
    return {"lb": num(rxn.lower_bound), "ub": num(rxn.upper_bound)}


def _load(inputs, name: str) -> dict:
    flux_data = load_flux_data(Path(inputs["data_dir"]) / name)
    return {
        "loaded": True,
        "conds": [str(c) for c in flux_data.conds],
        "p_tot": numlist(flux_data.p_tot),
        "gr_rate": numlist(flux_data.gr_rate),
        "exch_fluxes": [numlist(row) for row in flux_data.exch_fluxes],
        "exch_mets": [str(m) for m in flux_data.exch_mets],
        "exch_rxn_ids": [str(r) for r in flux_data.exch_rxn_ids],
    }


def _constrain(inputs, name: str, condition, max_min_growth: str = "max") -> dict:
    flux_data = load_flux_data(Path(inputs["data_dir"]) / name)
    model = build(inputs, constrain=False)
    apply_flux_data_constraints(
        model, flux_data,
        condition=condition, max_min_growth=max_min_growth, loose_strict_flux="loose",
        bio_rxn="R5", c_source="S1",
    )
    return {"S1": bounds(model, "S1"), "S2": bounds(model, "S2"), "R5": bounds(model, "R5")}


def run(ctx):
    inputs = ctx["inputs"]
    cobra.Configuration().solver = str(inputs["python_solver"])
    cases = {}

    cases["load_control"] = run_case(lambda: _load(inputs, "control_two_cond.tsv"))
    cases["load_ragged_row"] = run_case(lambda: _load(inputs, "ragged.tsv"))
    cases["load_fixed_columns_only"] = run_case(lambda: _load(inputs, "fixed_only.tsv"))
    cases["load_nested_parentheses_header"] = run_case(lambda: _load(inputs, "paren_header.tsv"))
    cases["load_missing_file"] = run_case(lambda: _load(inputs, "does_not_exist.tsv"))

    cases["constrain_control"] = run_case(lambda: _constrain(inputs, "control_two_cond.tsv", 0))
    cases["constrain_condition_by_name"] = run_case(lambda: _constrain(inputs, "control_two_cond.tsv", "cond2"))
    cases["constrain_unknown_condition"] = run_case(lambda: _constrain(inputs, "control_two_cond.tsv", "nope"))
    cases["constrain_missing_exchange"] = run_case(lambda: _constrain(inputs, "missing_exchange.tsv", 0))
    cases["constrain_invalid_max_min_growth"] = run_case(
        lambda: _constrain(inputs, "control_two_cond.tsv", 0, max_min_growth="median"))

    return cases
