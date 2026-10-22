"""Python side of the reporter-metabolites scenario on a hand-written model.

The result of each case is the list of tests the function returns, each with the metabolites in the
order returned and their scores. A standard deviation that is NaN (one neighbouring gene) is reported
as a flag plus a zero, because MATLAB's jsonencode writes NaN as null.
"""

import math

import cobra

from raven_toolbox.analysis import reporter_metabolites


def toy_model(reactions):
    model = cobra.Model("reporter")
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
        if spec["rule"]:
            reaction.gene_reaction_rule = str(spec["rule"])
    return model


def run(ctx):
    inputs = ctx["inputs"]
    results = {}
    for case in inputs["cases"]:
        model = toy_model(inputs["reactions"])
        pvalues = dict(zip(case["genes"], (math.nan if p is None else float(p) for p in case["pvalues"]), strict=True))
        fold_changes = None
        if case["fold_changes"]:
            fold_changes = dict(zip(case["genes"], (float(f) for f in case["fold_changes"]), strict=True))
        tests = []
        for result in reporter_metabolites(model, pvalues, gene_fold_changes=fold_changes):
            table = result.table
            std = [float(v) for v in table["std_z"]]
            tests.append(
                {
                    "test": result.test,
                    "mets": [str(m) for m in table["metabolite"]],
                    "z_scores": [float(v) for v in table["z_score"]],
                    "p_values": [float(v) for v in table["p_value"]],
                    "n_genes": [int(v) for v in table["n_genes"]],
                    "mean_z": [float(v) for v in table["mean_z"]],
                    "std_z": [0.0 if math.isnan(v) else v for v in std],
                    "std_is_nan": [math.isnan(v) for v in std],
                }
            )
        results[case["name"]] = {"n_results": len(tests), "tests": tests}
    return results
