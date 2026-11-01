"""Python side of the flexibilize edge-case scenario.

flexibilizeEnzConcs, relaxProteomicsGreedy and constrainEnzConcs on inputs the happy-path
scenarios do not reach: a growth target that a non-enzyme limit makes unreachable, a usage bound
that is already above the relaxed value, an infeasible start, two enzymes that bind at once, an
empty candidate set, a prior lower bound on a usage reaction and a negative concentration. A case
that raises is a result (`raised`), so the two sides can disagree about raising; `flex_control`
returns real content on both sides.
"""

import cobra

from geckopy import constrain_enz_concs, flexibilize_enz_concs, relax_proteomics_greedy
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


def _flex(model, inputs, exp_growth=None) -> dict:
    result = flexibilize_enz_concs(
        model,
        exp_growth=float(inputs["exp_growth"] if exp_growth is None else exp_growth),
        fold_change=float(inputs["fold_change"]),
        iter_per_enzyme=int(inputs["iter_per_enzyme"]),
        bio_rxn="R5",
        verbose=False,
    )
    flexed = sorted(
        (
            {
                "protein": str(result.uniprot_ids[i]),
                "old_conc": num(result.old_concs[i]),
                "flex_conc": num(result.flex_concs[i]),
                "ratio_incr": num(result.ratio_incr[i]),
                "frequence": num(result.frequence[i]),
            }
            for i in range(len(result.uniprot_ids))
        ),
        key=lambda r: r["protein"],
    )
    final = model.optimize()
    optimal = final.status == "optimal"
    return {
        "flexed": flexed,
        "final_is_optimal": bool(optimal),
        "final_growth": num(final.objective_value if optimal else 0.0),
        "p5_usage_bounds": bounds(model, "usage_prot_P5"),
    }


def _relax(model, ordered=True, **kwargs) -> dict:
    result = relax_proteomics_greedy(model, **kwargs)
    if not ordered:
        # Which enzyme is relaxed first follows the dual solution of the LP, so only the
        # outcome of the loop is recorded.
        return {
            "relaxed": sorted(str(k) for k in result.relaxed),
            "final_growth": num(result.final_growth),
            "converged": bool(result.converged),
            "n_steps": len(result.trace),
        }
    return {
        "relaxed": sorted(str(k) for k in result.relaxed),
        "final_growth": num(result.final_growth),
        "converged": bool(result.converged),
        "trace": [
            {
                "iteration": num(s.iteration),
                "relaxed": str(s.relaxed_uniprot),
                "growth_before": num(s.growth_before),
                "growth_after": num(s.growth_after),
                "shadow_price": num(s.shadow_price),
            }
            for s in result.trace
        ],
    }


def _growth(model) -> float:
    return float(model.optimize().objective_value)


def run(ctx):
    inputs = ctx["inputs"]
    cobra.Configuration().solver = str(inputs["python_solver"])
    cases = {}

    cases["flex_control"] = run_case(lambda: _flex(build(inputs), inputs))

    # A non-enzyme limit (the upper bound of R1, which carries the whole chain flux) caps growth
    # at 5% above the enzyme-limited growth, below exp_growth.
    def capped_growth():
        m = build(inputs)
        m.reactions.get_by_id("R1").upper_bound = 1.05 * _growth(m)
        return _flex(m, inputs)
    cases["flex_growth_capped_by_uptake"] = run_case(capped_growth)

    # The usage bound of P5 is already above conc * (1 + increase) while P5 still limits growth.
    def preset_bound():
        m = build(inputs)
        conc = float(m.ec.concs[enzyme_index(m, "P5")])
        m.reactions.get_by_id("usage_prot_P5").upper_bound = 10.0 * conc
        return _flex(m, inputs, exp_growth=float(inputs["exp_growth_high"]))
    cases["flex_usage_bound_above_target"] = run_case(preset_bound)

    def infeasible_start():
        m = build(inputs)
        m.reactions.get_by_id("R5").lower_bound = 50.0
        return _flex(m, inputs)
    cases["flex_infeasible_start"] = run_case(infeasible_start)

    # A prior lower bound on a usage reaction, and a negative concentration.
    def prior_lb():
        m = build(inputs, constrain=False)
        m.reactions.get_by_id("usage_prot_P5").lower_bound = 0.05
        m.ec.concs[enzyme_index(m, "P5")] = 0.01
        constrain_enz_concs(m)
        return {
            "bounds": bounds(m, "usage_prot_P5"),
            "pool_coefficient": num(m.reactions.get_by_id("usage_prot_P5").get_coefficient("prot_pool")),
        }
    cases["constrain_prior_lower_bound"] = run_case(prior_lb)

    def negative_conc():
        m = build(inputs, constrain=False)
        m.ec.concs[enzyme_index(m, "P4")] = -0.5
        constrain_enz_concs(m)
        return {"bounds": bounds(m, "usage_prot_P4")}
    cases["constrain_negative_concentration"] = run_case(negative_conc)

    cases["relax_control"] = run_case(
        lambda: _relax(build(inputs), minimal_growth=float(inputs["relax_target"])))

    def relax_tied():
        m = build(inputs)
        tie_p4(m)
        return _relax(m, ordered=False, minimal_growth=1.5 * _growth(m))
    cases["relax_tied_enzymes"] = run_case(relax_tied)

    cases["relax_empty_enzyme_set"] = run_case(
        lambda: _relax(build(inputs), minimal_growth=float(inputs["relax_target"]), enzyme_set=set()))

    return cases
