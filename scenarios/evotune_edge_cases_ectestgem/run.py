"""Python side of the evotune edge cases.

Screens and the CMA-ES run use a single-route ecTestGEM ecModel whose protein pool limits growth, with flux
data that caps the uptake of e1. A call that raises records only ``raised``.
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
    fill_eccodes_from_gem,
    load_conventional_gem,
    make_ec_model,
    set_prot_pool_size,
)
from geckopy.adapter.params import EvotuneParams, SourceGroupRule
from geckopy.databases.flux_data import load_flux_data
from geckopy.kcat_tuning.evotune import (
    EvotuneData,
    cmaes_kcat_tuning,
    screen_kcat_leverage,
    select_tunable_mask,
)
from geckopy.kcat_tuning.evotune.parsimony import (
    FrontierPoint,
    best_parsimonious,
    n_changed,
    parsimony_frontier,
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


def _num(x):
    """A float, or the text NaN, Inf or -Inf, which JSON cannot hold."""
    x = float(x)
    if math.isnan(x):
        return "NaN"
    if math.isinf(x):
        return "Inf" if x > 0 else "-Inf"
    return x


def _model(adapter, conv, inputs: dict, sources: dict):
    conv = conv.copy()
    for rxn_id in ("R2", "R4"):
        conv.reactions.get_by_id(rxn_id).bounds = (0.0, 0.0)

    model = make_ec_model(conv, adapter, gecko_light=False)
    fill_eccodes_from_gem(model)
    model.adapter = adapter

    index = {str(r): i for i, r in enumerate(model.ec.rxns)}
    for rxn_id, kcat in inputs["kcats"].items():
        model.ec.kcat[index[str(rxn_id)]] = float(kcat)
    for rxn_id, source in sources.items():
        model.ec.source[index[str(rxn_id)]] = str(source)

    model.solver = str(inputs["python_solver"])
    model.objective = str(inputs["bio_rxn"])
    apply_kcat_constraints(model)
    pool = inputs["pool"]
    set_prot_pool_size(model, p_tot=float(pool["p_tot"]), f=float(pool["f"]), sigma=float(pool["sigma"]))
    return model


def _params(inputs: dict, max_growth_weight: float) -> EvotuneParams:
    ev = inputs["evotune"]
    groups = list(ev["source_groups"])
    return EvotuneParams(
        source_groups={g: SourceGroupRule(sources=[g]) for g in groups},
        sigma0_log_source=dict(zip(groups, [float(s) for s in ev["sigma0_log_source"]])),
        sigma0_log_default=float(ev["sigma0_log_default"]),
        max_growth_weight=float(max_growth_weight),
        max_generations=int(ev["max_generations"]),
        rmse_threshold=float(ev["rmse_threshold"]),
        prior_penalty_weight=0.0,
        tie_isozymes=False,
    )


def _data(inputs: dict, flux_file: str, zero_flux: list) -> EvotuneData:
    return EvotuneData(
        flux_data=load_flux_data(inputs[flux_file]),
        max_grate=load_flux_data(inputs["max_growth_file"]),
        zero_flux=[str(r) for r in zero_flux],
    )


def _rows(screen) -> dict:
    rows = [
        {
            "rxn": str(r.rxn_id),
            "n_isozymes": int(r.n_isozymes),
            "source_group": str(r.source_group),
            "kcat0": _num(r.kcat0),
            "leverage": _num(r.leverage),
            "sigma0_log": _num(r.sigma0_log),
            "rank": _num(r.rank),
            "cum_leverage_share": _num(r.cum_leverage_share),
        }
        for r in screen.itertuples()
    ]
    return {
        "raised": False,
        "order": [r["rxn"] for r in rows],
        "rows": sorted(rows, key=lambda r: r["rxn"]),
    }


def _screen(adapter, conv, inputs, sources, params, data):
    model = _model(adapter, conv, inputs, sources)
    try:
        screen = screen_kcat_leverage(
            model, params=params, evotune_data=data, bio_rxn=str(inputs["bio_rxn"]), n_proc=1,
        )
    except Exception:
        return {"raised": True}, model, None
    return _rows(screen), model, screen


def _mask_ids(model, mask) -> list:
    return sorted(str(r) for r, keep in zip(model.ec.rxns, mask) if keep)


def _helpers(inputs: dict) -> dict:
    kcat = np.array(inputs["helper_kcat"], dtype=float)
    kcat0 = np.array(inputs["helper_kcat0"], dtype=float)
    sigma = np.array(inputs["helper_sigma0_log"], dtype=float)

    def score(vec):
        return float(np.mean(np.abs(np.log(vec / kcat))))

    points = parsimony_frontier(score, kcat, kcat0, sigma)
    best = best_parsimonious(points)

    nan_points = []
    for i, (rmse, n) in enumerate(zip(inputs["helper_points_rmse"], inputs["helper_points_n_changed"])):
        if i + 1 == int(inputs["helper_nan_index"]):
            rmse = float("nan")
        nan_points.append(FrontierPoint(
            threshold=float(i), kcat=kcat0.copy(), n_changed=int(n), mean_dev=0.0, rmse=float(rmse),
        ))
    try:
        nan_best = {"raised": False, "threshold": float(best_parsimonious(nan_points).threshold)}
    except Exception:
        nan_best = {"raised": True}

    return {
        "n_changed": int(n_changed(kcat, kcat0)),
        "frontier": [
            {
                "threshold": float(p.threshold),
                "n_changed": int(p.n_changed),
                "mean_dev": _num(p.mean_dev),
                "rmse": _num(p.rmse),
                "kcat": [_num(k) for k in p.kcat],
            }
            for p in points
        ],
        "best_threshold": float(best.threshold),
        "best_with_nan": nan_best,
    }


def run(ctx):
    inputs = ctx["inputs"]
    cobra.Configuration().solver = str(inputs["python_solver"])

    adapter = _adapter(inputs)
    conv = load_conventional_gem(adapter)
    ev = inputs["evotune"]

    params = _params(inputs, ev["max_growth_weight"])
    heavy = _params(inputs, ev["heavy_max_growth_weight"])
    data = _data(inputs, "flux_two_columns_file", [])

    results = {}
    control, model, screen = _screen(adapter, conv, inputs, inputs["control_sources"], params, data)
    results["screen_control"] = control
    # The model has no metabolite formulas; compute_excarbon warns and gives every exchange a carbon weight of 1.
    results["screen_no_formulas"] = control
    results["screen_source_suffix"] = _screen(adapter, conv, inputs, inputs["suffix_sources"], params, data)[0]
    results["screen_max_growth_weight"] = _screen(adapter, conv, inputs, inputs["control_sources"], heavy, data)[0]

    one = "flux_one_column_file"
    for name, zero in (("known", inputs["zero_flux_known"]), ("unknown", inputs["zero_flux_unknown"])):
        results[f"screen_zero_flux_{name}"] = _screen(
            adapter, conv, inputs, inputs["control_sources"], params, _data(inputs, one, zero),
        )[0]

    for share in inputs["target_impact_shares"]:
        key = "select_mask_" + str(float(share)).replace(".", "p")
        results[key] = _mask_ids(model, select_tunable_mask(model, screen, target_impact_share=float(share)))

    flat = _model(adapter, conv, inputs, inputs["control_sources"])
    mask = np.array([str(r) in set(inputs["flat_mask"]) for r in flat.ec.rxns], dtype=bool)
    try:
        tuned = cmaes_kcat_tuning(
            flat, params=params, evotune_data=data, bio_rxn=str(inputs["bio_rxn"]),
            tunable_mask=mask, seed=int(inputs["seed"]), n_proc=1, verbose=False,
        )
        results["cmaes_flat"] = {
            "raised": False,
            "rxns": sorted(str(r) for r in tuned.rxns),
            "n_generations": int(tuned.n_generations),
            "converged": bool(tuned.converged),
            "kcat_unchanged": bool(np.allclose(tuned.new_kcat, tuned.old_kcat)),
        }
    except Exception:
        results["cmaes_flat"] = {"raised": True}

    results["helpers"] = _helpers(inputs)
    return results
