"""Python side of the standard-kcat subsystem scenario.

assign_standard_kcat on ecTestGEM plus a reaction without genes, with the subsystem
assignments described in scenario.yml. A raise is a result.
"""

import importlib.util
import math
import sys
from pathlib import Path

import cobra

from geckopy import (
    ModelAdapter,
    assign_standard_kcat,
    fill_eccodes_from_gem,
    load_conventional_gem,
    load_uniprot_tsv,
    make_ec_model,
)


def num(x) -> dict:
    """A number as a record, so that NaN is a value both sides can state."""
    x = float(x)
    return {"is_nan": bool(math.isnan(x)), "value": 0.0 if math.isnan(x) else x}


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


def _add_r6(conv: cobra.Model) -> None:
    """R6, a reaction without genes."""
    r6 = cobra.Reaction("R6", name="R6")
    r6.lower_bound = 0.0
    r6.upper_bound = 1000.0
    r6.add_metabolites({conv.metabolites.get_by_id("m1c"): -1.0, conv.metabolites.get_by_id("m2c"): 1.0})
    conv.add_reactions([r6])


def _case(inputs: dict, subsystems: dict | None) -> dict:
    """The standard kcat step with `subsystems` (reaction id -> name, None for none at all)."""
    adapter = _adapter(inputs)
    conv = load_conventional_gem(adapter)
    _add_r6(conv)
    if subsystems is not None:
        for rxn in conv.reactions:
            rxn.subsystem = subsystems.get(rxn.id, subsystems["*"])
    model = make_ec_model(conv, adapter, gecko_light=False)
    fill_eccodes_from_gem(model)
    for i in range(len(model.ec.rxns)):
        model.ec.kcat[i] = 10.0
        model.ec.source[i] = "manual"
    uniprot = load_uniprot_tsv(inputs["uniprot_path"], id_type="taxonomy_id", auto_download=False)
    before = {str(r) for r in model.ec.rxns}
    assign_standard_kcat(model, uniprot)
    return {
        "new_rxns": sorted({str(r) for r in model.ec.rxns} - before),
        "rows": sorted(
            (
                {
                    "reaction": str(r),
                    "kcat": num(model.ec.kcat[i]),
                    "source": str(model.ec.source[i]),
                }
                for i, r in enumerate(model.ec.rxns)
            ),
            key=lambda e: e["reaction"],
        ),
    }


def _attempt(fn) -> dict:
    try:
        return {"raised": False, "value": fn()}
    except Exception as error:
        print(f"case raised: {type(error).__name__}: {error}")
        return {"raised": True, "value": {}}


def run(ctx):
    inputs = ctx["inputs"]
    return {
        "uniform": _attempt(lambda: _case(inputs, {"*": "Alpha"})),
        "case_variant": _attempt(lambda: _case(inputs, {"*": "Alpha", "R3": "alpha"})),
        "own_name": _attempt(lambda: _case(inputs, {"*": "Alpha", "R6": "Beta"})),
        "empty_subsystem": _attempt(lambda: _case(inputs, {"*": "Alpha", "R6": ""})),
        "no_subsystem": _attempt(lambda: _case(inputs, None)),
    }
