"""Python side of the genome-scale kcat tuning scenario.

sensitivity_tuning and fit_sigma on the full_ecModel tutorial's ecYeastGEM before its tuning.
See scenario.yml.
"""

import math
from pathlib import Path

import cobra

from geckopy import (
    ModelAdapter,
    fit_sigma,
    load_ec_model,
    sensitivity_tuning,
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


def fresh(inputs: dict):
    """The model before tuning with glucose uptake open and the pool as the adapter sets it."""
    cobra.Configuration().solver = str(inputs["python_solver"])
    folder = Path(inputs["adapter_python"])
    adapter = ModelAdapter.from_folder(folder)
    adapter.params.path = folder
    model = load_ec_model(str(inputs["model_file"]), adapter)
    model.solver = str(inputs["python_solver"])
    model.tolerance = float(inputs["python_tolerance"])
    model.adapter = adapter
    rxn = model.reactions.get_by_id(str(inputs["c_source"]))
    rxn.lower_bound = -1000.0
    model.objective = {model.reactions.get_by_id(str(inputs["bio_rxn"])): 1.0}
    model.objective_direction = "max"
    set_prot_pool_size(model)
    return model


def _tuned(model, result) -> dict:
    sol = model.optimize()
    return {
        "tuned": [
            {
                "reaction": str(result.rxns[i]),
                "enzymes": str(result.enzymes[i]),
                "old_kcat": num(result.old_kcat[i]),
                "new_kcat": num(result.new_kcat[i]),
                "source": str(result.source[i]),
            }
            for i in range(len(result.rxns))
        ],
        "final_growth": num(sol.objective_value),
        "kcat": [num(k) for k in model.ec.kcat],
        "notes": [str(n) if n else "" for n in model.ec.notes],
        "source": [str(s) if s else "" for s in model.ec.source],
    }


def _tune(inputs, *, growth, fold=None, ignore=None) -> dict:
    model = fresh(inputs)
    kwargs = {"desired_growth_rate": float(growth), "bio_rxn": str(inputs["bio_rxn"]), "verbose": False}
    if fold is not None:
        kwargs["fold_change"] = float(fold)
    if ignore is not None:
        kwargs["prot_to_ignore"] = list(ignore)
    base = float(model.optimize().objective_value)
    result = sensitivity_tuning(model, **kwargs)
    out = _tuned(model, result)
    out["base_growth"] = num(base)
    return out


def run(ctx):
    inputs = ctx["inputs"]
    results = {}

    default = run_case(lambda: _tune(inputs, growth=inputs["desired_growth"]))
    results["default"] = default

    def ignoring():
        first = default["value"]["tuned"][0]["enzymes"].split(";")[0]
        out = _tune(inputs, growth=inputs["desired_growth"], ignore=[first])
        out["ignored"] = first
        return out

    results["ignore_first"] = run_case(ignoring)
    results["low_growth_fold_2"] = run_case(
        lambda: _tune(inputs, growth=inputs["desired_growth_low"], fold=inputs["fold_change_low"])
    )

    def sigma():
        model = fresh(inputs)
        result = fit_sigma(
            model,
            growth_rate=float(inputs["sigma_growth_rate"]),
            p_tot=float(inputs["sigma_p_tot"]),
            f=float(inputs["sigma_f"]),
            method="grid",
        )
        return {
            "sigma": num(result.sigma),
            "model_pool_ub": num(model.reactions.get_by_id("prot_pool_exchange").upper_bound),
        }

    results["sigma_fitter"] = run_case(sigma)
    return results
