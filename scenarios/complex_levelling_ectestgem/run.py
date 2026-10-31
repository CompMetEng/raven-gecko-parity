"""Python side of the complex-levelling scenario.

complex_families, level_complex_concs and flexibilize_enz_concs(per_complex=...) on ecTestGEM with
complexes written into ec.rxn_enz_mat. A case that raises is a result (`raised`), so the two sides
can disagree about raising.
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
    complex_families,
    constrain_enz_concs,
    fill_enz_concs,
    fill_eccodes_from_gem,
    flexibilize_enz_concs,
    level_complex_concs,
    load_conventional_gem,
    load_prot_data,
    make_ec_model,
    set_prot_pool_size,
)

PROTEINS = ["P1", "P2", "P3", "P4", "P5"]


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


def build(inputs: dict):
    """The R2/R4-blocked single-route ecTestGEM of flexibilize_edge_cases_ectestgem, concentrations
    filled and not yet constrained."""
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
    return model


def finish(model) -> None:
    constrain_enz_concs(model)
    apply_kcat_constraints(model)
    set_prot_pool_size(model)


def set_layout(model, layout: dict) -> None:
    """Replace the rows of ec.rxn_enz_mat by the copies of the layout."""
    rxns = list(model.ec.rxns)
    enzymes = list(model.ec.enzymes)
    matrix = np.zeros(model.ec.rxn_enz_mat.shape)
    for rxn, subunits in layout.items():
        for enzyme, copies in (subunits or {}).items():
            matrix[rxns.index(str(rxn)), enzymes.index(str(enzyme))] = float(copies)
    model.ec.rxn_enz_mat = type(model.ec.rxn_enz_mat)(matrix)


def set_concs(model, **concs) -> None:
    enzymes = list(model.ec.enzymes)
    for protein, value in concs.items():
        model.ec.concs[enzymes.index(protein)] = float(value)


def concs_list(model) -> list:
    enzymes = list(model.ec.enzymes)
    return [{"protein": p, "conc": num(model.ec.concs[enzymes.index(p)])} for p in PROTEINS]


def usage_ub_list(model) -> list:
    return [
        {"protein": p, "ub": num(model.reactions.get_by_id("usage_prot_" + p).upper_bound)}
        for p in PROTEINS
    ]


def families_record(families, model) -> dict:
    enzymes = list(model.ec.enzymes)
    complexes = []
    for index, cplx in enumerate(families.complexes):
        owned = set(families.owned(index))
        complexes.append({
            "subunits": [
                {"enzyme": e, "copies": num(cplx[e]), "owned": bool(e in owned)}
                for e in sorted(cplx)
            ],
        })
    complexes.sort(key=lambda c: "\x01".join(s["enzyme"] for s in c["subunits"]))
    members = {}
    for index, cplx in enumerate(families.complexes):
        for e in families.owned(index):
            members[e] = "\x01".join(sorted(families.owned(index)))
    return {
        "complexes": complexes,
        "owner_members": [{"enzyme": e, "members": members.get(e, "")} for e in sorted(enzymes)],
        "shared": sorted(families.shared),
    }


def families_case(inputs: dict, layout: str) -> dict:
    model = build(inputs)
    set_layout(model, inputs["layouts"][layout])
    return families_record(complex_families(model), model)


def level_case(inputs: dict, layout: str, *, nan=(), **kwargs) -> dict:
    model = build(inputs)
    set_layout(model, inputs["layouts"][layout])
    for protein in nan:
        set_concs(model, **{protein: float("nan")})
    before = concs_list(model)
    result = level_complex_concs(model, **kwargs)
    after = concs_list(model)
    finish(model)
    return {
        "concs_before": before,
        "concs": after,
        "usage_ub": usage_ub_list(model),
        "adjusted": sorted(str(e) for e in result.adjusted),
        "added": sorted(str(e) for e in result.added),
        "shared": sorted(str(e) for e in result.shared),
        "skipped": sorted(sorted(str(e) for e in c) for c in result.skipped),
    }


def flex_case(inputs: dict, *, statistic, exp_growth: float, per_complex: bool) -> dict:
    """statistic None leaves the measured concentrations of the complex as they are."""
    model = build(inputs)
    set_layout(model, inputs["layouts"]["binding"])
    if statistic is not None:
        level_complex_concs(model, statistic=statistic)
    finish(model)
    start = float(model.optimize().objective_value)
    result = flexibilize_enz_concs(
        model,
        exp_growth=exp_growth,
        fold_change=float(inputs["fold_change"]),
        iter_per_enzyme=int(inputs["iter_per_enzyme"]),
        bio_rxn="R5",
        per_complex=per_complex,
        verbose=False,
    )
    ids = [str(p) for p in result.uniprot_ids]
    flexed = sorted(
        (
            {
                "protein": ids[i],
                "old_conc": num(result.old_concs[i]),
                "flex_conc": num(result.flex_concs[i]),
                "ratio_incr": num(result.ratio_incr[i]),
                "frequence": num(result.frequence[i]),
            }
            for i in range(len(ids))
        ),
        key=lambda r: r["protein"],
    )
    final = model.optimize()
    optimal = final.status == "optimal"
    return {
        "concs": concs_list(model),
        "start_growth": num(start),
        "flexed": flexed,
        "flex_set": sorted(ids),
        "flex_order": ids,
        "final_is_optimal": bool(optimal),
        "final_growth": num(final.objective_value if optimal else 0.0),
        "usage_ub": usage_ub_list(model),
    }


def run(ctx):
    inputs = ctx["inputs"]
    cobra.Configuration().solver = str(inputs["python_solver"])
    cases = {}

    for layout in ("pair", "trio", "shared", "isozyme", "binding"):
        cases["families_" + layout] = run_case(lambda layout=layout: families_case(inputs, layout))

    for layout in ("pair", "trio"):
        for statistic in ("median", "mean"):
            cases[f"level_{layout}_{statistic}"] = run_case(
                lambda layout=layout, statistic=statistic: level_case(inputs, layout, statistic=statistic))

    cases["level_shared_subunit"] = run_case(lambda: level_case(inputs, "shared"))
    cases["level_isozyme_variant"] = run_case(lambda: level_case(inputs, "isozyme"))
    cases["level_unmeasured_fill"] = run_case(
        lambda: level_case(inputs, "trio", nan=("P3",), fill_missing=True))
    cases["level_unmeasured_nofill"] = run_case(
        lambda: level_case(inputs, "trio", nan=("P3",), fill_missing=False))
    cases["level_complex_unmeasured"] = run_case(
        lambda: level_case(inputs, "trio", nan=("P1", "P2", "P3")))
    cases["level_unknown_statistic"] = run_case(
        lambda: level_case(inputs, "trio", statistic="mode"))

    targets = {
        "measured": (None, float(inputs["exp_growth"])),
        "median": ("median", float(inputs["exp_growth"])),
        "mean": ("mean", float(inputs["exp_growth_mean"])),
    }
    for label, (statistic, target) in targets.items():
        for per_complex in (False, True):
            name = f"flex_{label}_{'per_complex' if per_complex else 'default'}"
            cases[name] = run_case(
                lambda s=statistic, t=target, p=per_complex: flex_case(
                    inputs, statistic=s, exp_growth=t, per_complex=p))

    return cases
