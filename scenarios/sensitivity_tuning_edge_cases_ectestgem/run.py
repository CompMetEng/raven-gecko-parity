"""Python side of the sensitivity-tuning edge cases.

One checkpoint per case. A tuning run records the reactions it tuned (in the order the function returns
them and sorted by id), the final growth and the ec table of every reaction (kcat, source, notes). A call
that raises records only ``raised``.
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
    sensitivity_tuning,
    set_prot_pool_size,
    truncate_values,
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


def _model(adapter, conv, inputs: dict, kcats: dict, *, bounds=None, light=False, pool=None):
    conv = conv.copy()
    for rxn_id in ("R2", "R4"):
        conv.reactions.get_by_id(rxn_id).bounds = (0.0, 0.0)
    for rxn_id, (lower, upper) in (bounds or {}).items():
        conv.reactions.get_by_id(rxn_id).bounds = (lower, upper)

    model = make_ec_model(conv, adapter, gecko_light=light)
    fill_eccodes_from_gem(model)
    model.adapter = adapter

    for i, rxn_id in enumerate(model.ec.rxns):
        base = str(rxn_id).split("_", 1)[1] if light else str(rxn_id)
        if base in kcats:
            model.ec.kcat[i] = float(kcats[base])
            model.ec.source[i] = "manual"

    model.solver = str(inputs["python_solver"])
    model.objective = str(inputs["bio_rxn"])
    apply_kcat_constraints(model)
    if pool is None:
        set_prot_pool_size(model)
    else:
        set_prot_pool_size(model, p_tot=float(pool["p_tot"]), f=float(pool["f"]), sigma=float(pool["sigma"]))
    return model


def _ec_table(model) -> list:
    order = sorted(range(len(model.ec.rxns)), key=lambda i: str(model.ec.rxns[i]))
    return [
        {
            "reaction": str(model.ec.rxns[i]),
            "kcat": _num(model.ec.kcat[i]),
            "source": str(model.ec.source[i]),
            "notes": str(model.ec.notes[i]),
        }
        for i in order
    ]


def _tuned(result) -> list:
    return [
        {
            "reaction": str(result.rxns[i]),
            "enzymes": str(result.enzymes[i]),
            "old_kcat": _num(result.old_kcat[i]),
            "new_kcat": _num(result.new_kcat[i]),
            "source": str(result.source[i]),
        }
        for i in range(len(result.rxns))
    ]


def _tune(model, inputs: dict, desired: float, prot_to_ignore=None) -> dict:
    try:
        result = sensitivity_tuning(
            model,
            desired_growth_rate=float(desired),
            bio_rxn=str(inputs["bio_rxn"]),
            prot_to_ignore=prot_to_ignore,
            verbose=False,
        )
    except Exception:
        return {"raised": True}
    tuned = _tuned(result)
    return {
        "raised": False,
        "tuned": tuned,
        "tuned_sorted": sorted(tuned, key=lambda t: t["reaction"]),
        "final_growth": _num(model.optimize().objective_value),
        "ec": _ec_table(model),
    }


def _truncate(inputs: dict) -> list:
    out = []
    for text in inputs["truncate_inputs"]:
        try:
            value = float(truncate_values(np.array([float(text)]))[0])
            out.append({"input": str(text), "raised": False, "output": _num(value)})
        except Exception:
            out.append({"input": str(text), "raised": True})
    return out


def run(ctx):
    inputs = ctx["inputs"]
    cobra.Configuration().solver = str(inputs["python_solver"])

    adapter = _adapter(inputs)
    conv = load_conventional_gem(adapter)

    results = {}

    model = _model(adapter, conv, inputs, inputs["order_kcats"])
    results["order"] = _tune(model, inputs, inputs["order_desired_growth"])
    results["rerun"] = _tune(model, inputs, inputs["rerun_desired_growth"])

    for name in ("whole", "long"):
        model = _model(adapter, conv, inputs, inputs[f"notes_{name}_kcats"])
        results[f"notes_{name}"] = _tune(model, inputs, inputs[f"notes_{name}_desired_growth"])

    model = _model(adapter, conv, inputs, inputs["ceiling_kcats"], pool=inputs["ceiling_pool"])
    results["ceiling"] = _tune(model, inputs, inputs["ceiling_desired_growth"])

    for name, ignored in (
        ("unknown_id", inputs["ignore_unknown_id"]),
        ("lower_case_id", inputs["ignore_lower_case_id"]),
    ):
        model = _model(adapter, conv, inputs, inputs["ignore_kcats"])
        results[f"ignore_{name}"] = _tune(model, inputs, inputs["ignore_desired_growth"], [str(ignored)])

    model = _model(adapter, conv, inputs, inputs["order_kcats"], bounds={"R5": (0.0, 0.0)})
    results["zero_growth"] = _tune(model, inputs, inputs["zero_growth_desired_growth"])

    model = _model(adapter, conv, inputs, inputs["light_kcats"], light=True)
    light = _tune(model, inputs, inputs["light_desired_growth"])
    results["light"] = {
        "raised": light["raised"],
        "n_tuned": 0 if light["raised"] else len(light["tuned"]),
    }

    results["truncate_ties"] = _truncate(inputs)
    return results
