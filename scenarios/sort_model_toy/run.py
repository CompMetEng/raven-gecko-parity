"""Python side of the sortModel scenario on a hand-written model.

Each case runs sort_model with the given flags on a fresh model and records the stoichiometry of the
reactions whose result does not depend on the random search, the order of the metabolites, the order of
the reactions when no search ran, the position of the reaction without a subsystem and the order score of
every subsystem before and after.
"""

import cobra

from raven_toolbox.utils import sort_model, subsystem_order_scores


def build_model(inputs):
    model = cobra.Model("sort_toy")
    metabolites = {
        str(spec["id"]): cobra.Metabolite(str(spec["id"]), name=str(spec["name"]), compartment=str(spec["compartment"]))
        for spec in inputs["mets"]
    }
    model.add_metabolites(list(metabolites.values()))
    for spec in inputs["reactions"]:
        reaction = cobra.Reaction(str(spec["id"]))
        model.add_reactions([reaction])
        reaction.add_metabolites({metabolites[str(e["met"])]: float(e["coefficient"]) for e in spec["stoichiometry"]})
        reaction.bounds = (float(spec["lb"]), float(spec["ub"]))
        reaction.subsystem = str(spec["subsystem"])
    return model


def run_case(inputs, case):
    model = build_model(inputs)
    before = subsystem_order_scores(model)
    kwargs = {
        "sort_reversible": bool(case["sort_reversible"]),
        "sort_met_name": bool(case["sort_met_name"]),
        "sort_reaction_order": bool(case["sort_reaction_order"]),
    }
    if case["seed"] is not None:
        kwargs["seed"] = int(case["seed"])
    sort_model(model, **kwargs)
    searched = kwargs["sort_reaction_order"]
    coefficients = {
        r.id: sorted(
            ({"met": m.id, "coefficient": float(c)} for m, c in r.metabolites.items() if c != 0),
            key=lambda e: e["met"],
        )
        for r in model.reactions
        if not (searched and r.reversibility)
    }
    out = {
        "scores_before": before,
        "scores_after": subsystem_order_scores(model),
        "metabolite_order": [m.id for m in model.metabolites],
        "reactions_sorted": sorted(r.id for r in model.reactions),
        "coefficients": coefficients,
        "no_subsystem_positions": {r.id: i + 1 for i, r in enumerate(model.reactions) if not r.subsystem},
    }
    if not searched:
        out["reaction_order"] = [r.id for r in model.reactions]
    return out


def run(ctx):
    inputs = ctx["inputs"]
    return {case["name"]: run_case(inputs, case) for case in inputs["cases"]}
