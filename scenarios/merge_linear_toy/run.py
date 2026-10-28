"""Python side of the mergeLinear scenario on hand-written models.

The result is the partition of the original reactions into merged groups, the reactions of
the group id returned for every original reaction, the reactions of the reduced model and the
reactions reversed by the merge, all sorted.
"""

import cobra

from raven_toolbox.init import merge_linear


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
    results = {}
    for case in ctx["inputs"]["cases"]:
        reduced, orig_ids, group_ids, reversed_flags = merge_linear(toy_model("toy", case["model"]))
        members = {}
        for rxn_id, group in zip(orig_ids, group_ids, strict=True):
            if group:
                members.setdefault(group, []).append(rxn_id)
        partition = sorted((sorted(ids) for ids in members.values()), key=lambda ids: ids[0])
        results[case["name"]] = {
            "groups": [{"members": ids} for ids in partition],
            "group_ids": [
                {"reaction": rxn_id, "group_id": int(group)}
                for rxn_id, group in sorted(zip(orig_ids, group_ids, strict=True))
            ],
            "reactions": sorted(rxn.id for rxn in reduced.reactions),
            "reversed_reactions": sorted(
                rxn_id for rxn_id, flag in zip(orig_ids, reversed_flags, strict=True) if flag
            ),
        }
    return results
