"""Python side of the light Human-GEM ecModel analysis scenario.

Solving, mapping and subsetting on the tutorial's gecko-light ecModel.
See scenario.yml.
"""

import math
import re
import tempfile
from pathlib import Path

import cobra
import numpy as np
from raven_toolbox.io import read_yaml_model

from geckopy import (
    ModelAdapter,
    get_subset_ec_model,
    load_conventional_gem,
    load_ec_model,
    map_rxns_to_conv,
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
    folder = Path(inputs["adapter_python"])
    adapter = ModelAdapter.from_folder(folder)
    adapter.params.path = folder
    adapter.params.conv_gem = folder / "models" / "human-GEM.yml"
    return adapter


def _repaired_model_file(inputs: dict, workdir: str) -> str:
    """A copy of the tutorial ecModel without its empty `enzymes: !!omap` entries.

    The file has 52 ec reactions with no enzyme written that way, which geckopy's reader
    rejects (see ec_model_load_light_humangem). The copy is read by both sides.
    """
    text = (Path(inputs["adapter_python"]) / "models" / str(inputs["model_file"])).read_text(encoding="utf-8")
    text = re.sub(r"^ *- enzymes: !!omap\n(?= {0,4}\S)", "", text, flags=re.MULTILINE)
    path = Path(workdir) / "repaired.yml"
    path.write_text(text, encoding="utf-8", newline="")
    return str(path)


def _fresh(inputs, adapter):
    model = load_ec_model(_REPAIRED, adapter)
    model.solver = str(inputs["python_solver"])
    model.tolerance = float(inputs["python_tolerance"])
    model.adapter = adapter
    model.objective = {model.reactions.get_by_id(adapter.params.bio_rxn): 1.0}
    model.objective_direction = "max"
    return model


def _lp_growth(model) -> float:
    return float(model.optimize().objective_value)


def _subset(inputs, adapter, ec_model) -> dict:
    small = read_yaml_model(str(inputs["subset_file"]))
    small.solver = str(inputs["python_solver"])
    small.tolerance = float(inputs["python_tolerance"])
    small_growth = _lp_growth(small)
    sub = get_subset_ec_model(ec_model, small)
    sub.solver = str(inputs["python_solver"])
    sub.tolerance = float(inputs["python_tolerance"])
    bio = sub.reactions.get_by_id(adapter.params.bio_rxn)
    sub.objective = {bio: 1.0}
    sub.objective_direction = "max"
    return {
        "small_growth": num(small_growth),
        "ec_growth": num(_lp_growth(ec_model)),
        "subset_growth": num(_lp_growth(sub)),
        "counts": {
            "reactions": len(sub.reactions),
            "metabolites": len(sub.metabolites),
            "genes": len(sub.genes),
            "ec_rxns": len(sub.ec.rxns),
            "enzymes": len(sub.ec.enzymes),
        },
        "reactions": sorted(r.id for r in sub.reactions),
        "ec_rxns": [str(r) for r in sub.ec.rxns],
        "enzymes": [str(e) for e in sub.ec.enzymes],
    }


_REPAIRED = ""


def run(ctx):
    with tempfile.TemporaryDirectory(prefix="parity_light_analysis_") as workdir:
        global _REPAIRED
        _REPAIRED = _repaired_model_file(ctx["inputs"], workdir)
        return _run(ctx)


def _run(ctx):
    inputs = ctx["inputs"]
    cobra.Configuration().solver = str(inputs["python_solver"])
    adapter = _adapter(inputs)
    conv = load_conventional_gem(adapter)
    bio_id = adapter.params.bio_rxn
    model = _fresh(inputs, adapter)
    results = {}

    sol = model.optimize()
    results["growth"] = run_case(
        lambda: {
            "objective": num(sol.objective_value),
            "pool": num(sol.fluxes["prot_pool_exchange"]),
        }
    )

    def pool_min():
        m = _fresh(inputs, adapter)
        optimum = float(m.optimize().objective_value)
        bio = m.reactions.get_by_id(bio_id)
        bio.lower_bound = 0.99 * optimum
        m.objective = {m.reactions.get_by_id("prot_pool_exchange"): -1.0}
        m.objective_direction = "max"
        s = m.optimize()
        return {"pool_usage": num(-s.objective_value), "growth_lower_bound": num(bio.lower_bound)}

    results["min_pool"] = run_case(pool_min)

    def mapping():
        # A fixed flux vector, (i * 7919 mod 1013) / 1013 - 0.5 for the i-th reaction (from 1),
        # given to both sides: the LP solution differs between the two (see scenario.yml).
        n = len(model.reactions)
        vector = np.array([((i * 7919) % 1013) / 1013.0 - 0.5 for i in range(1, n + 1)])
        mapped = map_rxns_to_conv(model, conv, vector)
        return {
            "mapped": [{"reaction": r.id, "flux": num(v)} for r, v in zip(conv.reactions, np.ravel(mapped.mapped_flux))],
            "usage": sorted(
                ({"label": str(u), "flux": num(v)} for u, v in zip(mapped.usage_enz, np.ravel(mapped.enz_usage_flux))),
                key=lambda r: r["label"],
            ),
        }

    results["mapping"] = run_case(mapping)
    results["subset"] = run_case(lambda: _subset(inputs, adapter, _fresh(inputs, adapter)))

    def pool_sizes():
        m = _fresh(inputs, adapter)
        default = float(set_prot_pool_size(m))
        default_ub = float(m.reactions.get_by_id("prot_pool_exchange").upper_bound)
        halved = float(set_prot_pool_size(m, p_tot=float(adapter.params.p_tot) * float(inputs["pool_scale"])))
        return {
            "default_ub": num(default_ub),
            "halved_ub": num(m.reactions.get_by_id("prot_pool_exchange").upper_bound),
            "halved_growth": num(_lp_growth(m)),
        }

    results["pool_sizes"] = run_case(pool_sizes)

    return results
