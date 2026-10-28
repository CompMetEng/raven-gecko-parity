"""Python side of the growth-floor and LP gap-filling scenario on hand-written models.

Each case is a draft with an objective, a template and a minimum growth; the result is
which template reactions ``fill_gaps_milp`` (or ``fill_gaps_fast_lp``) adds.
"""

import cobra

from raven_toolbox.gapfilling import fill_gaps_fast_lp, fill_gaps_milp

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
        draft.objective = {
            draft.reactions.get_by_id(str(rid)): float(c)
            for rid, c in zip(case["objective"], case["objective_coeffs"])
        }
        try:
            if case["algorithm"] == "fastlp":
                result = fill_gaps_fast_lp(draft, template, verbose=False, n_proc=1)
                results[case["name"]] = {
                    "added_reactions": sorted(result.added_reactions),
                    "cannot_connect": sorted(result.cannot_connect),
                }
            else:
                result = fill_gaps_milp(
                    draft,
                    template,
                    min_growth=float(case["min_growth"]),
                    big_m=float(case["big_m"]),
                    verbose=False,
                )
                results[case["name"]] = {
                    "added_reactions": sorted(result.added_reactions),
                    "solved": result.exit_status in ("optimal", "feasible"),
                }
        except (RuntimeError, ValueError) as error:
            results[case["name"]] = {"error": str(error)}
    return results
