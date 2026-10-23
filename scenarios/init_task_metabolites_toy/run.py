"""Python side of the task-metabolite scenario for prepINITModel.

A four-metabolite model has the dead-end reaction R3 (B to T) as the only reaction of T, and one task names T
as its output. The result records the metabolites and reactions of the prepared reference model, the tasks the
prepared data keeps and the essential reactions.
"""

import cobra

from raven_toolbox.init import prep_init_model
from raven_toolbox.tasks import parse_task_list

cobra.Configuration().processes = 1


def build_model():
    model = cobra.Model("toy")
    model.compartments = {"c": "c", "e": "e"}
    mets = {}
    for name, comp in (("A", "e"), ("B", "c"), ("C", "c"), ("T", "c")):
        mets[name] = cobra.Metabolite(f"{name}_{comp}", name=name, compartment=comp)
    spec = [
        ("EXA", {"A": 1}),
        ("R1", {"A": -1, "B": 1}),
        ("R2", {"B": -1, "C": 1}),
        ("R3", {"B": -1, "T": 1}),
        ("EXC", {"C": -1}),
    ]
    for rid, stoich in spec:
        rxn = cobra.Reaction(rid)
        rxn.bounds = (-1000, 1000) if rid == "R2" else (0, 1000)
        model.add_reactions([rxn])
        rxn.add_metabolites({mets[n]: c for n, c in stoich.items()})
    return model


def run(ctx):
    inputs = ctx["inputs"]
    model = build_model()
    prep = prep_init_model(model, parse_task_list(inputs["tasks"]), ext_comp=str(inputs["ext_comp"]))
    ref = prep.ref_model
    return {
        "metabolites": sorted(f"{m.name}[{m.compartment}]" for m in ref.metabolites),
        "reactions": sorted(r.id for r in ref.reactions),
        "kept_tasks": [str(t.id) for t in prep.tasks],
        "essential_reactions": sorted(prep.essential_rxns),
    }
