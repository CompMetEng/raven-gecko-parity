"""Python side of the ec.source format scenario.

parse_kcat_source and format_kcat_source on the strings and calls of ``inputs``, and the source
that set_kcat_for_reactions writes on ecTestGEM. Every result holds only text and numbers so that
the MATLAB side can state the same document: a number is a ``{is_nan, value}`` record, a detail is
a list of ``{key, kind, number, text}`` sorted by key.
"""

import importlib.util
import math
import sys
from pathlib import Path

from geckopy import (
    ModelAdapter,
    fill_eccodes_from_gem,
    format_kcat_source,
    load_conventional_gem,
    make_ec_model,
    parse_kcat_source,
    set_kcat_for_reactions,
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
    x = float("nan") if x is None else float(x)
    return {"is_nan": bool(math.isnan(x)), "value": 0.0 if math.isnan(x) else x}


def run_case(fn) -> dict:
    """Result of one case; a raise is a result on both sides, not a failure."""
    try:
        value = fn()
    except Exception:
        return {"raised": True, "value": {}}
    return {"raised": False, "value": value}


def parse_record(text: str) -> dict:
    info = parse_kcat_source(text)
    detail = []
    for key in sorted(info.detail):
        value = info.detail[key]
        if isinstance(value, str):
            detail.append({"key": key, "kind": "text", "number": num(None), "text": value})
        else:
            detail.append({"key": key, "kind": "number", "number": num(value), "text": ""})
    return {
        "origin": info.origin,
        "detail": detail,
        "tuned": list(info.tuned),
        "prior": num(info.prior),
    }


def format_record(spec: dict) -> dict:
    return {
        "text": format_kcat_source(
            spec.get("origin", ""),
            spec.get("detail") or None,
            spec.get("tuned") or (),
            spec.get("prior"),
        )
    }


def set_kcat_record(inputs: dict) -> dict:
    spec = inputs["set_kcat"]
    adapter = _adapter(inputs)
    conv = load_conventional_gem(adapter)
    model = make_ec_model(conv, adapter, gecko_light=False)
    model.solver = str(inputs["python_solver"])
    fill_eccodes_from_gem(model)
    model.adapter = adapter
    index = list(model.ec.rxns).index(str(spec["rxn"]))
    model.ec.source[index] = str(spec["source_before"])
    others = [s for i, s in enumerate(model.ec.source) if i != index]
    set_kcat_for_reactions(model, [str(spec["rxn"])], float(spec["kcat"]))
    return {
        "source": str(model.ec.source[index]),
        "kcat": num(model.ec.kcat[index]),
        "others_unchanged": bool(
            others == [s for i, s in enumerate(model.ec.source) if i != index]
        ),
    }


def run(ctx):
    inputs = ctx["inputs"]
    cases = {}
    for i, text in enumerate(inputs["parse"]):
        cases[f"parse_{i:02d}"] = run_case(lambda t=text: parse_record(t))
    for i, spec in enumerate(inputs["format"]):
        cases[f"format_{i:02d}"] = run_case(lambda s=spec: format_record(s))
    cases["set_kcat_manual"] = run_case(lambda: set_kcat_record(inputs))
    return cases
