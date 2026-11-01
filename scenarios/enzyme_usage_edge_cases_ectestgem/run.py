"""Python side of the enzyme-usage edge-case scenario.

Every case builds its own ecModel (enzyme_usage_ectestgem's R2/R4-blocked single-route
fixture) and records what getEnzymeUsage, reportEnzymeUsage, getEnzymeBottlenecks,
getPfbaEnzymes and getConcControlCoeffs return for an input the happy-path scenario does
not reach. A case that raises is a result (`raised`), so the two sides can disagree about
whether to raise. The control cases return real content on both sides.
"""

import numpy as np

import cobra

from geckopy import (
    constrain_enz_concs,
    enzyme_usage,
    get_conc_control_coeffs,
    get_enzyme_bottlenecks,
    pfba_enzymes,
    report_enzyme_usage,
)
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


def _usage_rows(usage) -> list[dict]:
    return sorted(
        (
            {
                "protein": str(p),
                "ub": num(ub),
                "abs_usage": num(a),
                "cap_usage": num(c),
            }
            for p, ub, a, c in zip(usage.prot_id, usage.ub, usage.abs_usage, usage.cap_usage)
        ),
        key=lambda r: r["protein"],
    )


def _table_rows(df, value_col: str) -> list[dict]:
    return sorted(
        (
            {"protein": str(row.prot_id), "abs_usage": num(row.abs_usage), value_col: num(getattr(row, value_col))}
            for row in df.itertuples()
        ),
        key=lambda r: r["protein"],
    )


def _report(model, usage, inputs, *, high_cap=None, top=None) -> dict:
    report = report_enzyme_usage(
        model, usage,
        high_cap_usage=float(inputs["high_cap_usage"] if high_cap is None else high_cap),
        top_abs_usage=int(inputs["top_abs_usage"] if top is None else top),
    )
    return {
        "high_cap_usage": _table_rows(report.high_cap_usage, "cap_usage"),
        "top_abs_usage": _table_rows(report.top_abs_usage, "perc_usage"),
        "total_usage_flux": num(report.total_usage_flux),
    }


def _bottlenecks(model, top: int, fold_change: float | None = None, price: bool = False) -> dict:
    kwargs = {} if fold_change is None else {"fold_change": fold_change}
    df = get_enzyme_bottlenecks(model, top=top, **kwargs)
    order = [str(u) for u in df.index]
    rows = []
    for u, r in df.iterrows():
        row = {
            "protein": str(u),
            "objective_gain": num(r.objective_gain),
            "flux": num(r.flux),
            "cap_usage": num(r.cap_usage),
            "upper_bound": num(r.upper_bound),
            "multiple_active": bool(r.multiple_active),
        }
        if price:
            row["shadow_price"] = num(r.shadow_price)
        rows.append(row)
    rows.sort(key=lambda row: row["protein"])
    return {
        "order": order,
        "n_resolved": num(df.attrs["n_resolved"]),
        "active_constraints": num(df.attrs["active_constraints"]),
        "rows": rows,
    }


def _pfba(model, **kwargs) -> dict:
    optimum = float(model.optimize().objective_value)
    sol = pfba_enzymes(model, **kwargs)
    usage_sum = sum(abs(float(sol.fluxes[r.id])) for r in model.reactions if r.id.startswith("usage_prot_"))
    growth = float(sol.fluxes["R5"])
    return {
        "growth": num(growth),
        "growth_over_optimum": num(growth / optimum),
        "usage_sum": num(usage_sum),
    }


def _control_coeffs(model) -> list[dict]:
    enz, coeffs = get_conc_control_coeffs(model, fold_change=2.0, limit=0.0)
    return sorted(
        (
            {"protein": str(p), "analysed": bool(e), "coeff": num(c)}
            for p, e, c in zip(model.ec.enzymes, enz, coeffs)
        ),
        key=lambda r: r["protein"],
    )


def run(ctx):
    inputs = ctx["inputs"]
    cobra.Configuration().solver = str(inputs["python_solver"])
    cases = {}

    model = build(inputs)
    sol = model.optimize()
    usage = enzyme_usage(model, sol.fluxes)
    cases["control_usage"] = run_case(lambda: {"objective": num(sol.objective_value), "proteins": _usage_rows(usage)})
    cases["control_report"] = run_case(lambda: _report(model, usage, inputs))
    cases["control_bottlenecks"] = run_case(lambda: _bottlenecks(model, int(inputs["top_abs_usage"]), price=True))
    cases["bottlenecks_fold_change_3"] = run_case(
        lambda: _bottlenecks(model, int(inputs["top_abs_usage"]), fold_change=float(inputs["fold_change_3"]))
    )
    cases["bottlenecks_fold_change_1"] = run_case(
        lambda: _bottlenecks(model, int(inputs["top_abs_usage"]), fold_change=1.0)
    )
    cases["control_pfba"] = run_case(lambda: _pfba(model))

    # A measured concentration of 0 gives a usage upper bound of 0.
    def zero_ub():
        m = build(inputs)
        m.ec.concs[enzyme_index(m, "P1")] = 0.0
        constrain_enz_concs(m)
        s = m.optimize()
        return {"proteins": _usage_rows(enzyme_usage(m, s.fluxes))}
    cases["zero_ub_enzyme"] = run_case(zero_ub)

    # Enzymes with zero usage left out of the usage table, then summarised.
    def exclude_zero():
        u = enzyme_usage(model, sol.fluxes, include_zero=False)
        return {"proteins": _usage_rows(u), "report": _report(model, u, inputs)}
    cases["usage_exclude_zero"] = run_case(exclude_zero)

    cases["report_negative_top"] = run_case(lambda: _report(model, usage, inputs, top=-1))
    cases["report_top_one"] = run_case(lambda: _report(model, usage, inputs, top=1))

    # An enzyme without a usage reaction.
    def missing_usage_reaction():
        m = build(inputs)
        m.reactions.get_by_id("usage_prot_P1").remove_from_model()
        s = m.optimize()
        return {"proteins": _usage_rows(enzyme_usage(m, s.fluxes))}
    cases["usage_missing_reaction"] = run_case(missing_usage_reaction)

    cases["pfba_fraction_09"] = run_case(lambda: _pfba(model, fraction_of_optimum=0.9))
    cases["pfba_unknown_objective"] = run_case(lambda: _pfba(model, objective="NOSUCH"))

    def pfba_no_objective():
        m = build(inputs)
        for rxn in m.reactions:
            rxn.objective_coefficient = 0.0
        return _pfba(m)
    cases["pfba_no_objective"] = run_case(pfba_no_objective)

    # P4 and P5 both bind: the shadow prices of the two enzymes are not unique.
    def tied_bottlenecks():
        m = build(inputs)
        tie_p4(m)
        return _bottlenecks(m, int(inputs["top_abs_usage"]))
    cases["tied_bottlenecks"] = run_case(tied_bottlenecks)

    # Both caps bind, so each enzyme at its cap has an objective_gain of 0 (the other cap still
    # holds) and the LP gives a nonzero dual to at least one of them.
    def dual_without_gain():
        m = build(inputs)
        tie_p4(m)
        df = get_enzyme_bottlenecks(m, top=100)
        at_cap = df[(df.upper_bound > 0) & (df.cap_usage >= 1 - 1e-6)]
        return {
            "n_at_cap": num(len(at_cap)),
            "max_objective_gain_at_cap": num(at_cap.objective_gain.max()),
            "dual_nonzero_at_cap": bool((at_cap.shadow_price.abs() > 1e-9).any()),
        }
    cases["dual_without_gain"] = run_case(dual_without_gain)

    # Protein pool binding, no enzyme concentration constraints.
    def pool_binding_coeffs():
        m = build(inputs, constrain=False)
        pool = m.reactions.get_by_id("prot_pool_exchange")
        s = m.optimize()
        pool.upper_bound = 0.5 * abs(float(s.fluxes["prot_pool_exchange"]))
        return {"coeffs": _control_coeffs(m)}
    cases["control_coeffs_pool_binding"] = run_case(pool_binding_coeffs)

    return cases
