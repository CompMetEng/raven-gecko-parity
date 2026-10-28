"""Python side of the rescaleModelForINIT scenario."""

import cobra

from raven_toolbox.init.prep import rescale_for_init


def _build(reactions):
    model = cobra.Model("toy")
    mets = {}
    for spec in reactions:
        rxn = cobra.Reaction(str(spec["id"]))
        model.add_reactions([rxn])
        coeffs = {}
        for name, coef in spec["coefficients"].items():
            if name not in mets:
                mets[name] = cobra.Metabolite(name, compartment="c")
            coeffs[mets[name]] = float(coef)
        rxn.add_metabolites(coeffs)
        rxn.bounds = (float(spec["lb"]), float(spec["ub"]))
    return model


def run(ctx):
    inputs = ctx["inputs"]
    results = {}
    for case in inputs["cases"]:
        model = _build(inputs["reactions"])
        rescale_for_init(model, float(case["max_stoich_diff"]))
        results[case["name"]] = [
            {
                "id": rxn.id,
                "coefficients": [
                    {"met": met.id, "coefficient": float(coef)}
                    for met, coef in sorted(rxn.metabolites.items(), key=lambda kv: kv[0].id)
                ],
                "lb": float(rxn.lower_bound),
                "ub": float(rxn.upper_bound),
            }
            for rxn in model.reactions
        ]
    return results
