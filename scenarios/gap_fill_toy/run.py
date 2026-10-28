"""Python side of the connectivity gap-filling scenario on hand-written models.

Each case is a draft, a template and optional template reaction scores; the result is
what ``connect_blocked_reactions`` adds, connects and cannot connect. Defaults are used
for everything the case does not name, so a difference in a default shows up here.
"""

import cobra

from raven_toolbox.gapfilling import connect_blocked_reactions

cobra.Configuration().processes = 1


def toy_model(model_id, reactions):
    model = cobra.Model(model_id)
    metabolites = {}
    for spec in reactions:
        reaction = cobra.Reaction(str(spec["id"]))
        model.add_reactions([reaction])
        stoichiometry = {}
        for side, sign in (("substrates", -1), ("products", 1)):
            for name in spec[side] or []:
                if name not in metabolites:
                    metabolites[name] = cobra.Metabolite(str(name), name=str(name), compartment="c")
                stoichiometry[metabolites[name]] = sign
        reaction.add_metabolites(stoichiometry)
        reaction.bounds = (float(spec["lb"]), float(spec["ub"]))
    return model


def run(ctx):
    inputs = ctx["inputs"]
    cobra.Configuration().solver = str(inputs["python_solver"])
    results = {}
    for case in inputs["cases"]:
        draft = toy_model("draft", case["draft"])
        template = toy_model("template", case["template"])
        scores = case.get("scores") or []
        kwargs = {"allow_net_production": bool(case["net_production"])}
        if scores:
            kwargs["scores"] = {rxn["id"]: float(s) for rxn, s in zip(case["template"], scores)}
        try:
            result = connect_blocked_reactions(draft, template, **kwargs)
        except RuntimeError as error:
            results[case["name"]] = {"error": str(error)}
            continue
        results[case["name"]] = {
            "added_reactions": sorted(result.added_reactions),
            "newly_connected": sorted(result.newly_connected),
            "cannot_connect": sorted(result.cannot_connect),
        }
    return results
