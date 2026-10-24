"""Python side of the contractModel scenario: remove_duplicate_reactions on hand-written models."""

import warnings

import cobra

from raven_toolbox.manipulation import remove_duplicate_reactions


def clauses(rule):
    """Top-level OR clauses of a gene rule, without redundant outer parentheses, sorted."""
    out, depth, current = [], 0, ""
    tokens = str(rule).replace("(", " ( ").replace(")", " ) ").split()
    for token in tokens:
        if token == "(":
            depth += 1
        elif token == ")":
            depth -= 1
        if token == "or" and depth == 0:
            out.append(current)
            current = ""
        else:
            current += token
    if current:
        out.append(current)
    cleaned = []
    for clause in out:
        while clause.startswith("(") and clause.endswith(")") and clause.count("(") == 1:
            clause = clause[1:-1]
        cleaned.append(clause)
    return sorted(set(cleaned))


def build_model(mets, reactions):
    model = cobra.Model("contract_toy")
    metabolites = {
        str(m["id"]): cobra.Metabolite(str(m["id"]), name=str(m["name"]), compartment=str(m["compartment"]))
        for m in mets
    }
    model.add_metabolites(list(metabolites.values()))
    for spec in reactions:
        reaction = cobra.Reaction(str(spec["id"]))
        model.add_reactions([reaction])
        reaction.add_metabolites({metabolites[str(e["met"])]: float(e["coefficient"]) for e in spec["stoichiometry"]})
        reaction.bounds = (float(spec["lb"]), float(spec["ub"]))
        reaction.objective_coefficient = float(spec["c"])
        reaction.gene_reaction_rule = str(spec["rule"])
    return model


def run(ctx):
    inputs = ctx["inputs"]
    results = {}
    for case in inputs["cases"]:
        model = build_model(inputs["mets"], case["reactions"])
        restrict = None
        if case["mets"]:
            restrict = [model.metabolites.get_by_id(str(m)) for m in case["mets"]]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            removed = remove_duplicate_reactions(model, restrict_to_metabolites=restrict)
        results[case["name"]] = {
            "reactions": [
                {
                    "id": r.id,
                    "lb": float(r.lower_bound),
                    "ub": float(r.upper_bound),
                    "c": float(r.objective_coefficient),
                    "rule": clauses(r.gene_reaction_rule),
                }
                for r in sorted(model.reactions, key=lambda r: r.id)
            ],
            "removed": sorted(removed),
        }
    return results
