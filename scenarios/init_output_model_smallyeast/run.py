"""Python side of the prepINITModel and ftINIT scenario on smallYeast with a task list.

The template is read, its boundary reactions are opened and a task list is
prepared. ``ftinit`` then extracts a model from gene levels that alternate between a score of +1 and -1
in sorted gene order. The result records the tasks the prepared data keeps, the merged essential
reactions, the number of reference-model reactions whose bounds or stoichiometry differ from the template,
and the bounds and equation of three reactions in the extracted model, and the number of genes of the
reference model and the extracted model. A task list in which two tasks share an id keeps the one that can be performed.
"""

import cobra
import numpy as np

from raven_toolbox.init import TranscrData, ftinit, prep_init_model
from raven_toolbox.io import read_yaml_model
from raven_toolbox.tasks import parse_task_list

cobra.Configuration().processes = 1


def equation(rxn):
    return sorted(
        ({"met": m.id, "coefficient": float(c)} for m, c in rxn.metabolites.items()), key=lambda e: e["met"]
    )


def run(ctx):
    inputs = ctx["inputs"]
    cobra.Configuration().solver = str(inputs["python_solver"])
    model = read_yaml_model(inputs["model"])
    model.solver = str(inputs["python_solver"])
    for rid in inputs["open_bounds"]:
        model.reactions.get_by_id(str(rid)).upper_bound = float(inputs["open_ub"])
    template = {r.id: (float(r.lower_bound), float(r.upper_bound), equation(r)) for r in model.reactions}

    # A task list in which two tasks share an id: the infeasible one is removed by its position.
    repeated = prep_init_model(model, parse_task_list(inputs["tasks_duplicate"]), ext_comp=str(inputs["ext_comp"]))

    tasks = parse_task_list(inputs["tasks"])
    prep = prep_init_model(model, tasks, ext_comp=str(inputs["ext_comp"]))

    genes = sorted(g.id for g in model.genes)
    levels = np.array(
        [[float(inputs["gene_level_up"] if i % 2 == 0 else inputs["gene_level_down"])] for i in range(len(genes))]
    )
    transcr = TranscrData(genes, [str(inputs["tissue"])], [""], levels, float(inputs["threshold"]))
    out = ftinit(prep, str(inputs["tissue"]), None, None, transcr=transcr)

    n_diff = 0
    for rxn in prep.ref_model.reactions:
        lb, ub, eq = template[rxn.id]
        mine = equation(rxn)
        same_eq = [e["met"] for e in mine] == [e["met"] for e in eq] and all(
            abs(a["coefficient"] - b["coefficient"]) <= 1e-12 for a, b in zip(mine, eq)
        )
        if not (same_eq and float(rxn.lower_bound) == lb and float(rxn.upper_bound) == ub):
            n_diff += 1

    reported = []
    for rid in inputs["reported_reactions"]:
        rid = str(rid)
        if rid in out.reactions:
            r = out.reactions.get_by_id(rid)
            reported.append(
                {
                    "reaction": rid,
                    "in_model": True,
                    "lower_bound": float(r.lower_bound),
                    "upper_bound": float(r.upper_bound),
                    "equation": equation(r),
                }
            )
        else:
            reported.append({"reaction": rid, "in_model": False, "lower_bound": 0.0, "upper_bound": 0.0, "equation": []})

    return {
        "n_tasks_parsed": len(tasks),
        "kept_tasks": [{"id": str(t.id), "description": str(t.description)} for t in prep.tasks],
        "essential_reactions": sorted(prep.essential_rxns),
        "repeated_id_kept_tasks": [{"id": str(t.id), "description": str(t.description)} for t in repeated.tasks],
        "n_ref_reactions": len(prep.ref_model.reactions),
        "n_ref_genes": len(prep.ref_model.genes),
        "n_output_genes": len(out.genes),
        "n_ref_reactions_differing_from_template": n_diff,
        "n_output_reactions": len(out.reactions),
        "output_reactions": sorted(r.id for r in out.reactions),
        "reported": reported,
    }
