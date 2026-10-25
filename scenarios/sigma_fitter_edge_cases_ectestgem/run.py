"""Python side of the sigma-fit edge cases.

Every case builds a fresh single-route ecModel (R2 and R4 blocked), optionally edits the conventional
GEM first, and records what fit_sigma returns and the state it leaves the model in: the fitted sigma,
the protein pool bound and whether the biomass reaction still has equal bounds. The fixed_biomass cases
also record the bounds themselves; elsewhere sigmaFitter narrows an infinite upper bound to 1000, which is
not the point of those cases.
"""

import importlib.util
import sys
from pathlib import Path

import cobra

from geckopy import (
    ModelAdapter,
    apply_kcat_constraints,
    fill_eccodes_from_gem,
    fit_sigma,
    load_conventional_gem,
    make_ec_model,
    set_prot_pool_size,
)


def _adapter(inputs: dict) -> ModelAdapter:
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


def _model(adapter, conv, inputs: dict, bounds: dict):
    """Single-route ecModel; ``bounds`` maps a reaction id to (lower, upper) on the conventional GEM."""
    conv = conv.copy()
    for rxn_id in ("R2", "R4"):
        conv.reactions.get_by_id(rxn_id).bounds = (0.0, 0.0)
    for rxn_id, (lower, upper) in bounds.items():
        conv.reactions.get_by_id(rxn_id).bounds = (lower, upper)

    model = make_ec_model(conv, adapter, gecko_light=False)
    fill_eccodes_from_gem(model)
    model.adapter = adapter

    index = {r: i for i, r in enumerate(model.ec.rxns)}
    for rxn_id, kcat in inputs["kcats"].items():
        i = index[str(rxn_id)]
        model.ec.kcat[i] = float(kcat)
        model.ec.source[i] = "manual"

    model.solver = str(inputs["python_solver"])
    model.objective = str(inputs["bio_rxn"])
    apply_kcat_constraints(model)
    set_prot_pool_size(model)
    return model


def _case(adapter, conv, inputs: dict, bounds: dict, growth: float, method, report_bounds: bool) -> dict:
    model = _model(adapter, conv, inputs, bounds)
    kwargs = {} if method is None else {"method": method}
    result = fit_sigma(
        model,
        growth_rate=float(growth),
        p_tot=float(inputs["p_tot"]),
        f=float(inputs["f"]),
        **kwargs,
    )
    pool = model.reactions.get_by_id("prot_pool_exchange")
    bio = model.reactions.get_by_id(str(inputs["bio_rxn"]))
    out = {
        "sigma": float(result.sigma),
        "model_pool_ub": float(pool.upper_bound),
        "objective_fixed": bool(bio.lower_bound == bio.upper_bound),
    }
    if report_bounds:
        out["objective_lb"] = float(bio.lower_bound)
        out["objective_ub"] = float(bio.upper_bound)
    return out


def run(ctx):
    inputs = ctx["inputs"]
    cobra.Configuration().solver = str(inputs["python_solver"])

    adapter = _adapter(inputs)
    conv = load_conventional_gem(adapter)

    uptake = {"S1": (float(inputs["uptake_lower_bound"]), 1000.0)}
    fixed = {str(inputs["bio_rxn"]): (float(inputs["fixed_biomass"]),) * 2}

    cases = {}
    for name, bounds, growth in (
        ("on_grid", {}, inputs["growth_on_grid"]),
        ("between_grid", {}, inputs["growth_between_grid_points"]),
        ("uptake_cap", uptake, inputs["growth_above_uptake_limit"]),
        ("fixed_biomass", fixed, inputs["growth_between_grid_points"]),
    ):
        # method="grid" is the algorithm of sigmaFitter; "default" leaves the choice to fit_sigma.
        report = name == "fixed_biomass"
        cases[f"{name}_default"] = _case(adapter, conv, inputs, bounds, growth, None, report)
        cases[f"{name}_grid"] = _case(adapter, conv, inputs, bounds, growth, "grid", report)
    return cases
