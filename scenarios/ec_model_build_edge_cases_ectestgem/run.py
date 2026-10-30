"""Python side of the makeEcModel edge-case scenario.

One checkpoint per case. Each case edits a fresh copy of the conventional ecTestGEM, builds a
full and a light ecModel and reports a summary of each; a model the builder refuses is
reported as `raised: true` with an empty summary, since which exception it was is not the
finding. The summary sorts every list, because the order of ec.rxns after an unusual edit is
not something either side promises.
"""

from pathlib import Path

import cobra

from geckopy import ModelAdapter, load_conventional_gem, make_ec_model


def _bound(value: float) -> tuple[float, str]:
    value = float(value)
    if value == float("inf"):
        return 0.0, "+inf"
    if value == float("-inf"):
        return 0.0, "-inf"
    return value, "finite"


def _adapter(inputs: dict) -> ModelAdapter:
    adapter = ModelAdapter.from_folder(inputs["adapter_python"])
    fixture = Path(inputs["fixture_dir"])
    adapter.params.path = fixture
    adapter.params.conv_gem = fixture / "models" / "testModel.xml"
    return adapter


def _edit(case: str, model: "cobra.Model", adapter: ModelAdapter) -> None:
    """Apply a case's edit in place; the same edit is made in the .m file."""
    if case == "partial_complex":
        model.reactions.get_by_id("R3").gene_reaction_rule = "G4 and GX"
    elif case == "nested_gpr":
        model.reactions.get_by_id("R2").gene_reaction_rule = "G1 and (G2 or G3)"
    elif case == "neg_bounds":
        model.reactions.get_by_id("R3").bounds = (-10.0, -1.0)
    elif case == "multi_sink":
        sink = cobra.Reaction("SINK2", name="SINK2", lower_bound=-1000.0, upper_bound=1000.0)
        model.add_reactions([sink])
        sink.add_metabolites({model.metabolites.get_by_id("m1c"): -1.0,
                              model.metabolites.get_by_id("m2c"): -1.0})
    elif case == "rxn_id_rev_suffix":
        model.reactions.get_by_id("R1").id = "R1_REV"
    elif case == "rxn_id_exp_suffix":
        model.reactions.get_by_id("R1").id = "R1_EXP_1"
    elif case == "rxn_id_usage_prefix":
        model.reactions.get_by_id("R5").id = "usage_prot_P5"
    elif case == "met_id_prot_prefix":
        model.metabolites.get_by_id("m2c").id = "prot_m2"
    elif case == "enzyme_comp_is_id":
        model.compartments = {"c": "cytosol", "e": "extra"}
    elif case == "enzyme_comp_unknown":
        adapter.params.enzyme_comp = "zz"
    elif case == "no_genes":
        for rxn in model.reactions:
            rxn.gene_reaction_rule = ""
        cobra.manipulation.delete.remove_genes(model, list(model.genes), remove_reactions=False)
    elif case != "baseline":
        raise ValueError(f"unknown case {case}")


def _summary(model) -> dict:
    reactions = []
    for rxn in sorted(model.reactions, key=lambda r: r.id):
        lower, lower_kind = _bound(rxn.lower_bound)
        upper, upper_kind = _bound(rxn.upper_bound)
        reactions.append({
            "id": rxn.id, "lower_bound": lower, "lower_kind": lower_kind,
            "upper_bound": upper, "upper_kind": upper_kind,
            "objective_coefficient": float(rxn.objective_coefficient),
        })
    stoichiometry = sorted(
        ({"reaction": rxn.id, "metabolite": met.id, "coefficient": float(coef)}
         for rxn in model.reactions for met, coef in rxn.metabolites.items()),
        key=lambda e: (e["reaction"], e["metabolite"]),
    )
    ec = model.ec
    matrix = ec.rxn_enz_mat
    if hasattr(matrix, "toarray"):
        matrix = matrix.toarray()
    coupling = sorted(
        ({"reaction": ec.rxns[i], "enzyme": ec.enzymes[j], "coefficient": float(matrix[i][j])}
         for i in range(len(ec.rxns)) for j in range(len(ec.enzymes)) if matrix[i][j] != 0),
        key=lambda e: (e["reaction"], e["enzyme"]),
    )
    enzymes = sorted(
        ({"enzyme": str(e), "gene": str(g), "mw": float(m)}
         for e, g, m in zip(ec.enzymes, ec.genes, ec.mw)),
        key=lambda e: e["enzyme"],
    )
    return {
        "gecko_light": bool(ec.gecko_light),
        "n_reactions": len(model.reactions),
        "n_metabolites": len(model.metabolites),
        "reactions": reactions,
        "metabolites": sorted(met.id for met in model.metabolites),
        "stoichiometry": stoichiometry,
        "ec_rxns": sorted(str(r) for r in ec.rxns),
        "ec_enzymes": enzymes,
        "ec_coupling": coupling,
    }


def _build(case: str, inputs: dict, light: bool) -> dict:
    adapter = _adapter(inputs)
    model = load_conventional_gem(adapter).copy()
    _edit(case, model, adapter)  # a failing edit is a harness error, not a result
    try:
        ec_model = make_ec_model(model, adapter, gecko_light=light)
        return {"raised": False, "model": _summary(ec_model)}
    except Exception:  # noqa: BLE001 --- whether it raised is the finding, not the exception type
        return {"raised": True, "model": {}}


def run(ctx):
    inputs = ctx["inputs"]
    return {
        case: {"full": _build(case, inputs, False), "light": _build(case, inputs, True)}
        for case in inputs["cases"]
    }
