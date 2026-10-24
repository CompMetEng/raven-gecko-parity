"""Python side of the generateNewIds scenario on lists of existing identifiers."""

import cobra

from raven_toolbox.utils import generate_new_ids


def run(ctx):
    results = {}
    for case in ctx["inputs"]["cases"]:
        model = cobra.Model("ids")
        if case["type"] == "rxns":
            model.add_reactions([cobra.Reaction(str(i)) for i in case["ids"]])
        else:
            model.add_metabolites([cobra.Metabolite(str(i)) for i in case["ids"]])
        try:
            ids = generate_new_ids(
                model,
                "reactions" if case["type"] == "rxns" else "metabolites",
                str(case["prefix"]),
                quantity=int(case["quantity"]),
                num_length=int(case["num_length"]),
            )
        except ValueError:
            results[case["name"]] = {"error": True}
            continue
        results[case["name"]] = {"ids": list(ids)}
    return results
