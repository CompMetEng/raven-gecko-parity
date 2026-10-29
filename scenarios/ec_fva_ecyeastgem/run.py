"""Python side of the genome-scale ecFVA scenario.

ec_fva on ecYeastGEM constrained by the tutorial's flux data, once with growth free below
the measured rate and once with growth fixed at it. See scenario.yml.
"""

import math
from pathlib import Path

import cobra

from geckopy import (
    ModelAdapter,
    apply_flux_data_constraints,
    ec_fva,
    load_conventional_gem,
    load_ec_model,
    load_flux_data,
)


def num(x) -> dict:
    """A number as a record, so that NaN is a value both sides can state."""
    x = float(x)
    return {"is_nan": bool(math.isnan(x)), "value": 0.0 if math.isnan(x) else x}


def _rows(df) -> list[dict]:
    return [
        {"reaction": str(rxn_id), "min_flux": num(row.min_flux), "max_flux": num(row.max_flux)}
        for rxn_id, row in df.iterrows()
    ]


def run(ctx):
    inputs = ctx["inputs"]
    cobra.Configuration().solver = str(inputs["python_solver"])
    folder = Path(inputs["adapter_python"])
    adapter = ModelAdapter.from_folder(folder)
    adapter.params.path = folder
    adapter.params.conv_gem = folder / "models" / "yeast-GEM.yml"
    conv = load_conventional_gem(adapter)

    flux_data = load_flux_data(str(inputs["flux_data_file"]))
    flux_data.gr_rate[0] = float(inputs["growth_rate"])

    def build():
        model = load_ec_model(str(inputs["model_file"]), adapter)
        model.solver = str(inputs["python_solver"])
        model.adapter = adapter
        apply_flux_data_constraints(
            model,
            flux_data,
            condition=0,
            max_min_growth="max",
            loose_strict_flux="loose",
            bio_rxn=str(inputs["bio_rxn"]),
            c_source=str(inputs["c_source"]),
        )
        return model

    n_proc = int(inputs["python_n_proc"])
    model = build()
    # RAVEN's solveLP gives Gurobi FeasibilityTol and OptimalityTol of 1e-9; cobrapy's default
    # is 1e-6. The first two cases use RAVEN's value, the third the cobrapy default.
    default_tolerance = float(model.tolerance)
    model.tolerance = float(inputs["python_tolerance"])
    bio = model.reactions.get_by_id(str(inputs["bio_rxn"]))
    free = {
        "growth_bounds": [num(bio.lower_bound), num(bio.upper_bound)],
        "fva": _rows(ec_fva(model, conv, progress=False, n_proc=n_proc)),
    }

    bio.lower_bound = bio.upper_bound
    fixed = {
        "growth_bounds": [num(bio.lower_bound), num(bio.upper_bound)],
        "fva": _rows(ec_fva(model, conv, progress=False, n_proc=n_proc)),
    }

    model.tolerance = default_tolerance
    fixed_default = {
        "growth_bounds": [num(bio.lower_bound), num(bio.upper_bound)],
        "fva": _rows(ec_fva(model, conv, progress=False, n_proc=n_proc)),
    }
    return {"free": free, "fixed": fixed, "fixed_cobra_tolerance": fixed_default}
