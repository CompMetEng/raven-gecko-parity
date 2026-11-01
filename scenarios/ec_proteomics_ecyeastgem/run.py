"""Python side of the genome-scale proteomics scenario.

The stage-4 pipeline of GECKO's full_ecModel tutorial on ecYeastGEM, and loadProtData on the
tutorial proteomics and on an edited file. Every case records the result of the call and a
case that raises is a result. See scenario.yml.
"""

import math
from pathlib import Path

import cobra
import numpy as np

from geckopy import (
    ModelAdapter,
    apply_flux_data_constraints,
    calculate_f_factor,
    constrain_enz_concs,
    fill_enz_concs,
    flexibilize_enz_concs,
    get_conc_control_coeffs,
    load_ec_model,
    load_flux_data,
    load_prot_data,
    set_prot_pool_size,
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
    adapter.params.path = Path(inputs["adapter_python"])
    return adapter


def _solved(model) -> dict:
    """Whether the LP is feasible and its objective; an infeasible LP is a result."""
    sol = model.optimize()
    if sol.status != "optimal":
        return {"feasible": False, "objective": num(0.0)}
    return {"feasible": True, "objective": num(sol.objective_value)}


def fresh(inputs: dict):
    """The tutorial ecModel, growth as the objective."""
    cobra.Configuration().solver = str(inputs["python_solver"])
    adapter = _adapter(inputs)
    model = load_ec_model(str(inputs["model_file"]), adapter)
    model.solver = str(inputs["python_solver"])
    model.tolerance = float(inputs["python_tolerance"])
    model.adapter = adapter
    model.objective = {model.reactions.get_by_id(str(inputs["bio_rxn"])): 1.0}
    model.objective_direction = "max"
    return model


def _prot_rows(prot) -> dict:
    values = np.asarray(prot.abundances, dtype=float)
    if values.ndim == 1:
        values = values.reshape(-1, 1)
    return {
        "ids": [str(u) for u in prot.uniprot_ids],
        "abundances": [[num(v) for v in row] for row in values],
    }


def _exchange(model, inputs) -> list[dict]:
    return [
        {
            "reaction": str(r),
            "lower_bound": num(model.reactions.get_by_id(str(r)).lower_bound),
            "upper_bound": num(model.reactions.get_by_id(str(r)).upper_bound),
        }
        for r in inputs["exchange_rxns"]
    ]


def _usage_bounds(model) -> list[dict]:
    return sorted(
        (
            {
                "reaction": r.id,
                "lower_bound": num(r.lower_bound),
                "upper_bound": num(r.upper_bound),
            }
            for r in model.reactions
            if r.id.startswith("usage_prot_")
        ),
        key=lambda e: e["reaction"],
    )


def _concs(model) -> list[dict]:
    return [{"protein": str(e), "conc": num(c)} for e, c in zip(model.ec.enzymes, model.ec.concs)]


def _loader_cases(inputs: dict) -> dict:
    real = str(inputs["prot_data_file"])
    edge = str(inputs["prot_data_edge_file"])
    empty = str(inputs["prot_data_empty_file"])
    cases = {
        "real_default": lambda: _prot_rows(load_prot_data(real, [3])),
        "real_unfiltered": lambda: _prot_rows(load_prot_data(real, [3], filter_data=False)),
        "edge_default": lambda: _prot_rows(load_prot_data(edge, [3, 2])),
        "edge_unfiltered": lambda: _prot_rows(load_prot_data(edge, [3, 2], filter_data=False)),
        "edge_max_missing_vector": lambda: _prot_rows(
            load_prot_data(edge, [3, 2], max_missing=[0.5, 1.0])
        ),
        "edge_strict": lambda: _prot_rows(
            load_prot_data(edge, [3, 2], min_val=0.5, max_rsd=0.3, cut_lowest=0.0, add_stdevs=0.0)
        ),
        "edge_cut_half": lambda: _prot_rows(load_prot_data(edge, [3, 2], cut_lowest=50.0)),
        "edge_wrong_replicates": lambda: _prot_rows(load_prot_data(edge, [3, 3])),
        "empty_cells_default": lambda: _prot_rows(load_prot_data(empty, [3, 2])),
        "empty_cells_unfiltered": lambda: _prot_rows(load_prot_data(empty, [3, 2], filter_data=False)),
    }
    return {name: run_case(fn) for name, fn in cases.items()}


def _enzyme_concs(inputs: dict) -> dict:
    real = load_prot_data(str(inputs["prot_data_file"]), [3])
    edge_file = str(inputs["prot_data_edge_file"])

    def one(prot, data_col):
        if prot == "edge":
            prot = load_prot_data(edge_file, [3, 2])
        model = fresh(inputs)
        fill_enz_concs(model, prot, data_col=data_col)
        concs = _concs(model)
        constrain_enz_concs(model)
        return {
            "n_measured": int(np.sum(~np.isnan(model.ec.concs))),
            "concs": concs,
            "usage_bounds": _usage_bounds(model),
            "growth": _solved(model),
        }

    return {
        "real": run_case(lambda: one(real, 0)),
        "edge_condition_1": run_case(lambda: one("edge", 0)),
        "edge_condition_2": run_case(lambda: one("edge", 1)),
        "edge_out_of_range": run_case(lambda: one("edge", 2)),
    }


def _f_factor(inputs: dict) -> dict:
    real = load_prot_data(str(inputs["prot_data_file"]), [3])
    model = fresh(inputs)
    subset = [str(e) for e in model.ec.enzymes[:100]]
    return {
        "real": run_case(lambda: num(calculate_f_factor(model, real))),
        "real_unfiltered": run_case(
            lambda: num(calculate_f_factor(model, load_prot_data(str(inputs["prot_data_file"]), [3], filter_data=False)))
        ),
        "edge": run_case(
            lambda: num(calculate_f_factor(model, load_prot_data(str(inputs["prot_data_edge_file"]), [3, 2])))
        ),
        "real_subset": run_case(lambda: num(calculate_f_factor(model, real, enzymes=subset))),
    }


def _flux_data(inputs: dict) -> dict:
    data = load_flux_data(str(inputs["flux_data_file"]))
    parsed = {
        "conds": [str(c) for c in data.conds],
        "p_tot": [num(v) for v in data.p_tot],
        "gr_rate": [num(v) for v in data.gr_rate],
        "exch_fluxes": [[num(v) for v in row] for row in np.atleast_2d(data.exch_fluxes)],
        "exch_mets": [str(m) for m in data.exch_mets],
        "exch_rxn_ids": [str(r) for r in data.exch_rxn_ids],
    }

    def constrained(loose_strict):
        model = fresh(inputs)
        apply_flux_data_constraints(
            model,
            data,
            condition=0,
            max_min_growth="max",
            loose_strict_flux=loose_strict,
            bio_rxn=str(inputs["bio_rxn"]),
            c_source=str(inputs["c_source"]),
        )
        return {"exchange": _exchange(model, inputs), "growth": _solved(model)}

    return {
        "parsed": parsed,
        "loose": run_case(lambda: constrained("loose")),
        "pct_10": run_case(lambda: constrained(10.0)),
    }


def _pipeline(inputs: dict) -> dict:
    prot = load_prot_data(str(inputs["prot_data_file"]), [3])
    data = load_flux_data(str(inputs["flux_data_file"]))
    model = fresh(inputs)
    fill_enz_concs(model, prot)
    constrain_enz_concs(model)
    f = float(calculate_f_factor(model, prot))
    set_prot_pool_size(model, p_tot=float(data.p_tot[0]), f=f)
    apply_flux_data_constraints(
        model,
        data,
        condition=0,
        max_min_growth="max",
        loose_strict_flux="loose",
        bio_rxn=str(inputs["bio_rxn"]),
        c_source=str(inputs["c_source"]),
    )
    out = {
        "f_factor": num(f),
        "pool_ub": num(model.reactions.get_by_id("prot_pool_exchange").upper_bound),
        "growth_before": _solved(model),
    }

    def ccc():
        # The flag that says which enzymes were analysed follows the vertex of the LP; the
        # coefficient is 0 for an enzyme that was not analysed.
        _, coeffs = get_conc_control_coeffs(model, fold_change=2.0, limit=float(inputs["ccc_limit"]))
        rows = sorted(
            (
                {"protein": str(p), "coeff": num(c)}
                for p, c in zip(model.ec.enzymes, coeffs)
            ),
            key=lambda r: r["protein"],
        )
        return rows

    out["control_coeffs"] = run_case(ccc)

    def flex():
        result = flexibilize_enz_concs(
            model,
            exp_growth=float(data.gr_rate[0]),
            fold_change=float(inputs["fold_change"]),
            bio_rxn=str(inputs["bio_rxn"]),
            verbose=False,
        )
        return {
            "flexed": [
                {
                    "protein": str(result.uniprot_ids[i]),
                    "old_conc": num(result.old_concs[i]),
                    "flex_conc": num(result.flex_concs[i]),
                    "ratio_incr": num(result.ratio_incr[i]),
                    "frequence": int(result.frequence[i]),
                }
                for i in range(len(result.uniprot_ids))
            ],
            "usage_bounds": _usage_bounds(model),
            "pool_ub": num(model.reactions.get_by_id("prot_pool_exchange").upper_bound),
            "growth_after": _solved(model),
        }

    out["flexibilize"] = run_case(flex)
    return out


def run(ctx):
    inputs = ctx["inputs"]
    return {
        "loader": _loader_cases(inputs),
        "enzyme_concs": _enzyme_concs(inputs),
        "f_factor": _f_factor(inputs),
        "flux_data": _flux_data(inputs),
        "pipeline": run_case(lambda: _pipeline(inputs)),
    }
