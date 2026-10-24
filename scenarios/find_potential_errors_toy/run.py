"""Python side of the findPotentialErrors scenario."""

import cobra

from raven_toolbox.utils import find_non_dnf_grrules


def run(ctx):
    rules = [str(rule) for rule in ctx["inputs"]["rules"]]
    model = cobra.Model("toy")
    for i, rule in enumerate(rules, start=1):
        rxn = cobra.Reaction(f"R{i}")
        model.add_reactions([rxn])
        rxn.gene_reaction_rule = rule
    return {"flagged": sorted(issue.reaction_id for issue in find_non_dnf_grrules(model))}
