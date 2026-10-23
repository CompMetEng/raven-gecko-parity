"""Python side of the task-metabolite scenario for ftINIT.

The metabolite T is used only by the reactions R3 and R5, whose genes have low expression, so the
extraction removes both. The metabolite U has only the dead-end reaction R4, which the first simplification
removes. Tasks name T and U as outputs. The result records the metabolites and
reactions of the extracted model and the tasks the prepared data keeps.
"""

import cobra
import numpy as np

from raven_toolbox.init import TranscrData, ftinit, prep_init_model
from raven_toolbox.tasks import parse_task_list

cobra.Configuration().processes = 1


def build_model():
    model = cobra.Model("toy")
    model.compartments = {"c": "c", "e": "e"}
    mets = {}
    for name, comp in (("A", "e"), ("B", "c"), ("C", "c"), ("T", "c"), ("U", "c")):
        mets[name] = cobra.Metabolite(f"{name}_{comp}", name=name, compartment=comp)
    spec = [
        ("EXA", {"A": 1}, ""),
        ("R1", {"A": -1, "B": 1}, ""),
        ("R2", {"B": -1, "C": 1}, "G2"),
        ("R3", {"B": -1, "T": 1}, "G3"),
        ("R5", {"T": -1, "C": 1}, "G5"),
        ("EXC", {"C": -1}, ""),
        ("R4", {"B": -1, "U": 1}, ""),
    ]
    for rid, stoich, rule in spec:
        rxn = cobra.Reaction(rid)
        rxn.bounds = (-1000, 1000) if rid == "R2" else (0, 1000)
        model.add_reactions([rxn])
        rxn.add_metabolites({mets[n]: c for n, c in stoich.items()})
        if rule:
            rxn.gene_reaction_rule = rule
    return model


def run(ctx):
    inputs = ctx["inputs"]
    cobra.Configuration().solver = str(inputs["python_solver"])
    model = build_model()
    prep = prep_init_model(model, parse_task_list(inputs["tasks"]), ext_comp=str(inputs["ext_comp"]))
    genes = sorted(g.id for g in model.genes)
    levels = np.array([[float(inputs["levels"][g])] for g in genes])
    transcr = TranscrData(genes, [str(inputs["tissue"])], [""], levels, float(inputs["threshold"]))
    out = ftinit(prep, str(inputs["tissue"]), None, None, transcr=transcr)
    return {
        "metabolites": sorted(f"{m.name}[{m.compartment}]" for m in out.metabolites),
        "reactions": sorted(r.id for r in out.reactions),
        "kept_tasks": [str(t.id) for t in prep.tasks],
        "reference_metabolites": sorted(f"{m.name}[{m.compartment}]" for m in prep.ref_model.metabolites),
    }
