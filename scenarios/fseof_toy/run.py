"""Python side of the FSEOF scenario on a hand-written model."""

import cobra

from raven_toolbox.analysis import fseof


def build_model(inputs):
    model = cobra.Model("fseof_toy")
    metabolites = {
        str(m["id"]): cobra.Metabolite(str(m["id"]), name=str(m["name"]), compartment=str(m["compartment"]))
        for m in inputs["mets"]
    }
    model.add_metabolites(list(metabolites.values()))
    for spec in inputs["reactions"]:
        reaction = cobra.Reaction(str(spec["id"]))
        model.add_reactions([reaction])
        reaction.add_metabolites({metabolites[str(e["met"])]: float(e["coefficient"]) for e in spec["stoichiometry"]})
        reaction.bounds = (float(spec["lb"]), float(spec["ub"]))
        reaction.objective_coefficient = float(spec["c"])
    return model


def run(ctx):
    inputs = ctx["inputs"]
    results = {}
    for case in inputs["cases"]:
        model = build_model(inputs)
        result = fseof(
            model,
            str(inputs["target"]),
            biomass_rxn=str(inputs["biomass"]),
            n_steps=int(case["iterations"]),
            max_fraction=float(case["coefficient"]),
            correlation_threshold=float(case["threshold"]),
        )
        table = result.targets
        out = {"amplify": [], "knockdown": [], "knockout": [], "slopes": {}}
        for row in table.itertuples():
            out[row.target_type].append(row.reaction)
            out["slopes"][row.reaction] = round(float(row.slope), 4)
        for key in ("amplify", "knockdown", "knockout"):
            out[key] = sorted(out[key])
        results[case["name"]] = out
    return results
