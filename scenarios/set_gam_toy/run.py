"""Python side of the setGAM scenario on a hand-written model.

Each case sets the cofactors of a biomass reaction to one or more values on a fresh copy of the model and
records every coefficient of the biomass reaction, whether its magnitude equals the value exactly, and the
bounds of the maintenance reaction.
"""

import cobra

from raven_toolbox.biomass import set_gam


def build_model(mets, reactions):
    model = cobra.Model("gam_toy")
    metabolites = {
        str(spec["id"]): cobra.Metabolite(str(spec["id"]), name=str(spec["name"]), compartment=str(spec["compartment"]))
        for spec in mets
    }
    model.add_metabolites(list(metabolites.values()))
    for spec in reactions:
        reaction = cobra.Reaction(str(spec["id"]))
        model.add_reactions([reaction])
        reaction.add_metabolites({metabolites[str(e["met"])]: float(e["coefficient"]) for e in spec["stoichiometry"]})
        reaction.bounds = (float(spec["lb"]), float(spec["ub"]))
    return model


def run_one(inputs, case, value):
    model = build_model(inputs["mets"], case["reactions"])
    kwargs = {"biomass_rxn": str(case["biomass_rxn"]), "cofactor_met_names": [str(n) for n in case["cofactors"]]}
    if case["ngam_rxn"] is not None:
        kwargs["ngam_rxn"] = str(case["ngam_rxn"])
    if case["ngam_value"] is not None:
        kwargs["ngam_value"] = float(case["ngam_value"])
    try:
        set_gam(model, float(value), **kwargs)
    except (KeyError, ValueError):
        return {"value": float(value), "error": True}
    biomass = model.reactions.get_by_id(str(case["biomass_rxn"]))
    coefficients = sorted(
        (
            {
                "met": met.id,
                "coefficient": float(coef),
                "abs_equals_value": abs(coef) == abs(float(value)),
            }
            for met, coef in biomass.metabolites.items()
            if coef != 0
        ),
        key=lambda e: e["met"],
    )
    ngam = model.reactions.get_by_id("ATPM")
    return {
        "value": float(value),
        "coefficients": coefficients,
        "ngam_lb": float(ngam.lower_bound),
        "ngam_ub": float(ngam.upper_bound),
    }


def run(ctx):
    inputs = ctx["inputs"]
    return {
        case["name"]: {"runs": [run_one(inputs, case, value) for value in case["values"]]}
        for case in inputs["cases"]
    }
