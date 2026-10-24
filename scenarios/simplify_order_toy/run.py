"""Python side of the simplification-order scenario on hand-written models.

Each case is a small model, a list of simplifyModel modes and a list of reserved
reactions. The result is the model after the call, with reactions, metabolites and
genes in model order, because the surviving reaction of a merge and the position of
the others are part of what is compared.
"""

import cobra

from raven_toolbox.manipulation import simplify_model

cobra.Configuration().processes = 1

_MODES = {
    "groupLinear": "group_linear",
    "constrainReversible": "constrain_reversible",
    "deleteZeroInterval": "delete_zero_interval",
    "deleteDuplicates": "delete_duplicates",
    "deleteInaccessible": "delete_dead_end",
}


def as_list(value):
    return [] if value is None else [str(v) for v in value]


def toy_model(model_id, reactions, extra_metabolites):
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
        if spec["rule"]:
            reaction.gene_reaction_rule = str(spec["rule"])
    model.add_metabolites([cobra.Metabolite(str(m), name=str(m), compartment="c") for m in extra_metabolites or []])
    return model


def run(ctx):
    inputs = ctx["inputs"]
    cobra.Configuration().solver = str(inputs["python_solver"])
    results = {}
    for case in inputs["cases"]:
        model = toy_model("toy", case["reactions"], case["extra_metabolites"])
        flags = {_MODES[mode]: True for mode in as_list(case["modes"])}
        simplify_model(model, reserved=as_list(case["reserved"]), **flags)
        results[case["name"]] = {
            "reactions": [r.id for r in model.reactions],
            "metabolites": [m.id for m in model.metabolites],
            "genes": [g.id for g in model.genes],
            "bounds": [
                {"reaction": r.id, "lower_bound": float(r.lower_bound), "upper_bound": float(r.upper_bound)}
                for r in model.reactions
            ],
            "stoichiometry": [
                {"reaction": r.id, "metabolite": m.id, "coefficient": float(c)}
                for r in model.reactions
                for m, c in sorted(r.metabolites.items(), key=lambda item: item[0].id)
            ],
        }
    return results
