"""Python side of the measured-zero and isozyme-coverage scenario.

get_isozyme_coverage and zero_unmeasured_isozymes on ecTestGEM with ec.concs written per case
(NaN not measured, 0 measured absent), constrain_enz_concs and flexibilize_enz_concs with a
measured zero, and level_complex_concs with a zero subunit. A case that raises is a result
(`raised`), so the two sides can disagree about raising.
"""

import importlib.util
import logging
import math
import re
import sys
from pathlib import Path

import cobra
import numpy as np

from geckopy import (
    ModelAdapter,
    apply_kcat_constraints,
    constrain_enz_concs,
    fill_eccodes_from_gem,
    fill_enz_concs,
    flexibilize_enz_concs,
    get_isozyme_coverage,
    level_complex_concs,
    load_conventional_gem,
    load_prot_data,
    make_ec_model,
    set_prot_pool_size,
    zero_unmeasured_isozymes,
)

PROTEINS = ("P1", "P2", "P3", "P4", "P5")


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


def run_case(fn) -> dict:
    """Result of one case; a raise is a result on both sides, not a failure."""
    try:
        value = fn()
    except Exception:
        return {"raised": True, "value": {}}
    return {"raised": False, "value": value}


def build(inputs: dict, concs: str, *, constrain: bool = True, layout: str | None = None):
    """The R2/R4-blocked single-route ecTestGEM with the concentrations of the case."""
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
    for rxn_id, kcat in inputs["kcats"].items():
        i = index[str(rxn_id)]
        model.ec.kcat[i] = float(kcat)
        model.ec.source[i] = "manual"

    prot_data = load_prot_data(Path(adapter.params.path) / "data" / "proteomics.tsv", [1])
    fill_enz_concs(model, prot_data)
    if layout is not None:
        set_layout(model, inputs["layouts"][layout])
    enzymes = list(model.ec.enzymes)
    for protein, value in inputs["concs"][concs].items():
        model.ec.concs[enzymes.index(protein)] = np.nan if value is None else float(value)
    if constrain:
        constrain_enz_concs(model)
    apply_kcat_constraints(model)
    set_prot_pool_size(model)
    return model


def set_layout(model, layout: dict) -> None:
    """Replace the rows of ec.rxn_enz_mat by the copies of the layout."""
    rxns = list(model.ec.rxns)
    enzymes = list(model.ec.enzymes)
    matrix = np.zeros(model.ec.rxn_enz_mat.shape)
    for rxn, subunits in layout.items():
        for enzyme, copies in (subunits or {}).items():
            matrix[rxns.index(str(rxn)), enzymes.index(str(enzyme))] = float(copies)
    model.ec.rxn_enz_mat = type(model.ec.rxn_enz_mat)(matrix)


def concs_list(model) -> list:
    enzymes = list(model.ec.enzymes)
    return [{"protein": p, "conc": num(model.ec.concs[enzymes.index(p)])} for p in PROTEINS]


def usage_ub_list(model) -> list:
    return [
        {"protein": p, "ub": num(model.reactions.get_by_id("usage_prot_" + p).upper_bound)}
        for p in PROTEINS
    ]


def growth(model) -> dict:
    solution = model.optimize()
    return num(solution.objective_value if solution.status == "optimal" else 0.0)


def coverage_record(model) -> dict:
    coverage, zeroable = get_isozyme_coverage(model)
    rows = [
        {
            "rxn": str(row["rxn"]),
            "n_entries": int(row["n_entries"]),
            "n_measured": int(row["n_measured"]),
            "n_unmeasured": int(row["n_unmeasured"]),
            "measured": sorted(str(e) for e in row["measured_enzymes"]),
            "unmeasured": sorted(str(e) for e in row["unmeasured_enzymes"]),
            "shared": sorted(str(e) for e in row["shared_enzymes"]),
        }
        for _, row in coverage.iterrows()
    ]
    rows.sort(key=lambda r: r["rxn"])
    return {"rows": rows, "zeroable": sorted(str(e) for e in zeroable)}


def zero_record(model, **kwargs) -> dict:
    changed = zero_unmeasured_isozymes(model, **kwargs)
    records = [
        {"uniprot": str(u), "new_ub": num(v)}
        for u, v in zip(changed["uniprot"], changed["new_upper_bound"])
    ]
    records.sort(key=lambda r: r["uniprot"])
    return {"changed": records, "usage_ub": usage_ub_list(model), "concs": concs_list(model)}


def level_zero(inputs: dict) -> dict:
    model = build(inputs, "zero_subunit", constrain=False, layout="pair")
    result = level_complex_concs(model)
    return {
        "concs": concs_list(model),
        "adjusted": sorted(str(e) for e in result.adjusted),
        "added": sorted(str(e) for e in result.added),
    }


def flex_zero(inputs: dict) -> dict:
    model = build(inputs, "zero_p5")
    records = []

    class _Capture(logging.Handler):
        def emit(self, record):
            records.append(record.getMessage())

    handler = _Capture(level=logging.WARNING)
    logger = logging.getLogger("geckopy")
    logger.addHandler(handler)
    try:
        result = flexibilize_enz_concs(
            model,
            exp_growth=float(inputs["exp_growth"]),
            fold_change=float(inputs["fold_change"]),
            iter_per_enzyme=int(inputs["iter_per_enzyme"]),
            bio_rxn="R5",
            verbose=False,
        )
    finally:
        logger.removeHandler(handler)
    listed = []
    for message in records:
        found = re.search(r"enzymes (.*?) have a concentration of 0", message, re.IGNORECASE)
        if found:
            listed.extend(found.group(1).replace(" ", "").split(","))
    flexed = sorted(str(e) for e in result.uniprot_ids)
    return {
        "warned_enzymes": sorted(listed),
        "flexed": flexed,
        "p5_ub": num(model.reactions.get_by_id("usage_prot_P5").upper_bound),
        "final_growth": growth(model),
    }


def run(ctx):
    inputs = ctx["inputs"]
    limit = float(inputs["detection_limit"])
    cases = {}

    for name, concs, layout in (
        ("coverage_partial", "partial", None),
        ("coverage_measured_zero_entry", "measured_zero_entry", None),
        ("coverage_zero_subunit", "zero_subunit", None),
        ("coverage_all_measured", "all_measured", None),
        ("coverage_none_measured", "none_measured", None),
        ("coverage_shared", "partial", "shared"),
    ):
        cases[name] = run_case(
            lambda c=concs, l=layout: coverage_record(build(inputs, c, layout=l)))

    cases["zero_method_zero"] = run_case(
        lambda: zero_record(build(inputs, "partial"), method="zero"))
    cases["zero_detection_limit_default"] = run_case(
        lambda: zero_record(build(inputs, "partial"), method="detection_limit"))
    cases["zero_detection_limit_explicit"] = run_case(
        lambda: zero_record(build(inputs, "partial"), method="detection_limit", limit=limit))

    def not_raised():
        model = build(inputs, "partial")
        model.reactions.get_by_id("usage_prot_P3").upper_bound = 0.05
        return zero_record(model, method="detection_limit", limit=limit)
    cases["zero_bound_not_raised"] = run_case(not_raised)

    cases["zero_shared_enzyme_unchanged"] = run_case(
        lambda: zero_record(build(inputs, "partial", layout="shared"), method="zero"))
    cases["zero_no_partial_group"] = run_case(
        lambda: zero_record(build(inputs, "all_measured"), method="zero"))
    cases["zero_unknown_method"] = run_case(
        lambda: zero_record(build(inputs, "partial"), method="bogus"))
    cases["zero_negative_limit"] = run_case(
        lambda: zero_record(build(inputs, "partial"), method="detection_limit", limit=-1.0))
    cases["zero_no_concentration_no_limit"] = run_case(
        lambda: zero_record(build(inputs, "none_measured"), method="detection_limit"))

    def constrain_zero():
        model = build(inputs, "zero_p5")
        return {"usage_ub": usage_ub_list(model), "growth": growth(model)}
    cases["constrain_measured_zero"] = run_case(constrain_zero)

    cases["flexibilize_measured_zero"] = run_case(lambda: flex_zero(inputs))
    cases["level_zero_subunit"] = run_case(lambda: level_zero(inputs))
    return cases
