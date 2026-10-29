"""Python side of the published-ecYeastGEM analysis scenario.

The post-solve analysis functions run on ecYeastGEM 9.1.1 of the ecModels repository. Every
case records the result of the call and a case that raises is a result. See scenario.yml for
what is compared and why the flux vector is compared separately from the quantities that are
the same at every optimum.
"""

import math
from pathlib import Path

import cobra
import numpy as np
import pandas as pd

from geckopy import (
    ModelAdapter,
    apply_kcat_constraints,
    ec_fseof,
    enzyme_usage,
    get_conc_control_coeffs,
    get_enzyme_bottlenecks,
    load_ec_model,
    pfba_enzymes,
    report_enzyme_usage,
)


def num(x) -> dict:
    """A number as a record, so that NaN is a value both sides can state."""
    x = float(x)
    return {"is_nan": bool(math.isnan(x)), "value": 0.0 if math.isnan(x) else x}


def run_case(fn) -> dict:
    """Result of one case; a raise is a result on both sides, not a failure."""
    try:
        value = fn()
    except Exception as error:
        print(f"case raised: {type(error).__name__}: {error}")
        return {"raised": True, "value": {}}
    return {"raised": False, "value": value}


def _adapter(inputs: dict) -> ModelAdapter:
    adapter = ModelAdapter.from_folder(Path(inputs["adapter_python"]))
    folder = Path(inputs["adapter_python"])
    adapter.params.path = folder
    return adapter


def _set_objective(model, rxn_id: str, coefficient: float = 1.0) -> None:
    model.objective = {model.reactions.get_by_id(rxn_id): coefficient}
    model.objective_direction = "max"


def _set_bounds(model, rxn_id: str, lower=None, upper=None) -> None:
    rxn = model.reactions.get_by_id(rxn_id)
    if lower is not None and lower > rxn.upper_bound:
        rxn.upper_bound = float(upper if upper is not None else lower)
    if upper is not None and upper < rxn.lower_bound:
        rxn.lower_bound = float(lower if lower is not None else upper)
    if lower is not None:
        rxn.lower_bound = float(lower)
    if upper is not None:
        rxn.upper_bound = float(upper)


def fresh(inputs: dict):
    """The tutorial ecModel with glucose uptake open and growth as the objective."""
    cobra.Configuration().solver = str(inputs["python_solver"])
    adapter = _adapter(inputs)
    model = load_ec_model(str(inputs["model_file"]), adapter)
    model.solver = str(inputs["python_solver"])
    model.tolerance = float(inputs["python_tolerance"])
    model.adapter = adapter
    _set_bounds(model, str(inputs["c_source"]), lower=-1000.0)
    _set_objective(model, str(inputs["bio_rxn"]), 1.0)
    return model


def _usage_rows(usage) -> list[dict]:
    return sorted(
        (
            {"protein": str(p), "ub": num(ub), "abs_usage": num(a), "cap_usage": num(c)}
            for p, ub, a, c in zip(usage.prot_id, usage.ub, usage.abs_usage, usage.cap_usage)
        ),
        key=lambda r: r["protein"],
    )


def _usage_summary(usage) -> dict:
    """What every optimal vertex shares: the total usage and the largest usages as numbers."""
    values = sorted((float(a) for a in usage.abs_usage), reverse=True)
    return {"total": num(sum(values)), "largest": [num(v) for v in values[:5]]}


def _table_rows(df, value_col: str) -> list[dict]:
    return sorted(
        (
            {"protein": str(row.prot_id), "abs_usage": num(row.abs_usage), value_col: num(getattr(row, value_col))}
            for row in df.itertuples()
        ),
        key=lambda r: r["protein"],
    )


def _report(model, usage, inputs) -> dict:
    report = report_enzyme_usage(
        model,
        usage,
        high_cap_usage=float(inputs["high_cap_usage"]),
        top_abs_usage=int(inputs["top_abs_usage"]),
    )
    return {
        "high_cap_usage": _table_rows(report.high_cap_usage, "cap_usage"),
        "top_abs_usage": _table_rows(report.top_abs_usage, "perc_usage"),
        "top_abs_usage_order": [str(p) for p in report.top_abs_usage.prot_id],
        "total_usage_flux": num(report.total_usage_flux),
    }


def _bottlenecks(model, top: int) -> dict:
    df = get_enzyme_bottlenecks(model, top=top)
    order = [str(u) for u in df.index]
    rows = sorted(
        (
            {
                "protein": str(u),
                "shadow_price": num(r.shadow_price),
                "flux": num(r.flux),
                "cap_usage": num(r.cap_usage),
                "upper_bound": num(r.upper_bound),
            }
            for u, r in df.iterrows()
        ),
        key=lambda r: r["protein"],
    )
    return {"order": order, "rows": rows}


def _control_coeffs(model) -> list[dict]:
    enz, coeffs = get_conc_control_coeffs(model, fold_change=2.0, limit=0.001)
    return sorted(
        (
            {"protein": str(p), "analysed": bool(e), "coeff": num(c)}
            for p, e, c in zip(model.ec.enzymes, enz, coeffs)
        ),
        key=lambda r: r["protein"],
    )


def _exchange_fluxes(sol, inputs) -> list[dict]:
    return [{"reaction": str(r), "flux": num(sol.fluxes[str(r)])} for r in inputs["exchange_rxns"]]


def _growth(model, inputs):
    sol = model.optimize()
    return sol, {
        "objective": num(sol.objective_value),
        "exchange": _exchange_fluxes(sol, inputs),
    }


def _pfba(model) -> dict:
    optimum = float(model.optimize().objective_value)
    sol = pfba_enzymes(model)
    usage_sum = sum(abs(float(sol.fluxes[r.id])) for r in model.reactions if r.id.startswith("usage_prot_"))
    objective = float(sum(r.objective_coefficient * sol.fluxes[r.id] for r in model.reactions))
    return {
        "objective": num(objective),
        "objective_over_optimum": num(objective / optimum),
        "usage_sum": num(usage_sum),
    }


def _fseof(model, inputs) -> dict:
    result = ec_fseof(
        model,
        str(inputs["fseof_target"]),
        str(inputs["c_source"]),
        n_steps=int(inputs["fseof_steps"]),
        bio_rxn=str(inputs["bio_rxn"]),
    )
    ranked = [str(row.reaction) for row in result.targets.itertuples()]
    return {
        "enforced_levels": [num(v) for v in result.enforced],
        "targets": sorted(ranked),
        "targets_ranked": ranked,
    }


def _pick_edit_rows(ec) -> dict:
    """One rule, shared with the MATLAB side, from the order of ec.rxns."""
    rxns = [str(r) for r in ec.rxns]
    matrix = ec.rxn_enz_mat.tocsr()
    per_row = np.diff(matrix.indptr)
    kcat = np.asarray(ec.kcat, dtype=float)

    def first(mask):
        hits = np.flatnonzero(mask)
        return int(hits[0]) if len(hits) else -1

    iso = np.array(["_EXP_" in r for r in rxns])
    rev = np.array([r.endswith("_REV") for r in rxns])
    return {
        "isozyme": first(iso & (kcat > 0)),
        "reversible": first(rev & ~iso & (kcat > 0)),
        "complex": first((per_row > 1) & (kcat > 0) & ~iso & ~rev),
        "zero": first(~iso & ~rev & (per_row == 1) & (kcat > 0)),
    }


def _column(model, rxn_id: str) -> dict:
    rxn = model.reactions.get_by_id(rxn_id)
    lower, lower_kind = (0.0, "-inf") if rxn.lower_bound == float("-inf") else (float(rxn.lower_bound), "finite")
    upper, upper_kind = (0.0, "+inf") if rxn.upper_bound == float("inf") else (float(rxn.upper_bound), "finite")
    return {
        "reaction": rxn_id,
        "lower_bound": lower,
        "lower_kind": lower_kind,
        "upper_bound": upper,
        "upper_kind": upper_kind,
        "coefficients": sorted(
            ({"metabolite": m.id, "coefficient": float(c)} for m, c in rxn.metabolites.items()),
            key=lambda e: e["metabolite"],
        ),
    }


def _kcat_edits(inputs: dict) -> dict:
    model = fresh(inputs)
    picks = _pick_edit_rows(model.ec)
    factors = {"isozyme": 2.0, "reversible": 3.0, "complex": 0.5, "zero": 0.0}
    chosen = {name: model.ec.rxns[i] for name, i in picks.items() if i >= 0}
    before = {name: float(model.ec.kcat[i]) for name, i in picks.items() if i >= 0}
    for name, i in picks.items():
        if i >= 0:
            model.ec.kcat[i] = float(model.ec.kcat[i]) * factors[name]
    apply_kcat_constraints(model)
    sol = model.optimize()
    return {
        "picked": {name: str(r) for name, r in chosen.items()},
        "kcat_before": before,
        "columns": [_column(model, str(r)) for r in chosen.values()],
        "growth": num(sol.objective_value),
    }


def run(ctx):
    inputs = ctx["inputs"]
    results = {}
    model = fresh(inputs)
    sol, growth = _growth(model, inputs)
    results["growth"] = run_case(lambda: growth)

    # The usage functions are given the same flux vector on both sides, read
    # from data/growth_fluxes.tsv; the vector the LP returns differs between the two sides.
    table = pd.read_csv(str(inputs["fluxes_file"]), sep="\t")
    shared = pd.Series(table["flux"].values, index=table["reaction"].astype(str).values)
    usage = enzyme_usage(model, shared)
    results["usage"] = run_case(lambda: _usage_rows(usage))
    results["report"] = run_case(lambda: _report(model, usage, inputs))
    solved_usage = enzyme_usage(model, sol.fluxes)
    results["usage_solved"] = run_case(lambda: _usage_summary(solved_usage))
    results["bottlenecks"] = run_case(lambda: _bottlenecks(model, int(inputs["bottlenecks_top"])))
    results["pfba"] = run_case(lambda: _pfba(model))
    results["control_coeffs"] = run_case(lambda: _control_coeffs(model))

    # The tutorial's second condition: the smallest glucose uptake that still grows at 0.25.
    def min_uptake():
        m = fresh(inputs)
        _set_objective(m, str(inputs["c_source"]), 1.0)
        _set_bounds(m, str(inputs["bio_rxn"]), lower=0.25)
        s, g = _growth(m, inputs)
        u = enzyme_usage(m, s.fluxes)
        return {"growth": g, "usage": {"total": _usage_summary(u)["total"]}}

    results["min_uptake"] = run_case(min_uptake)
    results["fseof"] = run_case(lambda: _fseof(fresh(inputs), inputs))
    results["kcat_edits"] = run_case(lambda: _kcat_edits(inputs))
    return results
