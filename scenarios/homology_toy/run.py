"""Python side of the getModelFromHomology scenario on one-reaction templates.

Templates convert a_c to b_c; only the gene rules differ. The orthologs are fixed, so the
model-building logic is the only variable. The result per case is the reactions of the
draft with their gene rules as sorted DNF clauses.
"""

import cobra
import pandas as pd

from raven_toolbox.manipulation import gpr_to_dnf
from raven_toolbox.reconstruction.homology import get_model_from_homology, make_ortholog_hits


def template_model(spec):
    model = cobra.Model(str(spec["id"]))
    a = cobra.Metabolite("a_c", name="a", compartment="c")
    b = cobra.Metabolite("b_c", name="b", compartment="c")
    for reaction_spec in spec["reactions"]:
        reaction = cobra.Reaction(str(reaction_spec["id"]))
        reaction.add_metabolites({a: -1, b: 1})
        model.add_reactions([reaction])
        reaction.gene_reaction_rule = str(reaction_spec["rule"])
    return model


def run(ctx):
    results = {}
    for case in ctx["inputs"]["cases"]:
        templates = [template_model(spec) for spec in case["templates"]]
        frames = []
        for template in templates:
            pairs = [
                (str(p["source"]), str(p["target"]))
                for p in case["pairs"]
                if str(p["template"]) == template.id
            ]
            frames.append(make_ortholog_hits(pairs, template.id, str(case["target"])))
        hits = pd.concat(frames, ignore_index=True)
        order = [str(t) for t in case["preferred_order"]] or None
        try:
            draft = get_model_from_homology(
                templates,
                hits,
                str(case["target"]),
                complex_policy=str(case["complex_policy"]),
                preferred_order=order,
            ).model
        except (RuntimeError, ValueError) as error:
            results[case["name"]] = {"error": str(error)}
            continue
        reactions = []
        for reaction in sorted(draft.reactions, key=lambda r: r.id):
            clauses = sorted(sorted(clause) for clause in gpr_to_dnf(reaction.gpr))
            reactions.append({"reaction": reaction.id, "clauses": clauses})
        results[case["name"]] = {"n_reactions": len(reactions), "reactions": reactions}
    return results
