"""Python side of the prepINITModel and ftINIT scenario on RAVEN's ten-reaction test model.

Each case prepares the model (``prep_init_model``) and extracts a model from the expression
levels of one tissue (``ftinit`` with ``transcr``, the input RAVEN's ``ftINIT`` takes). The
result records the prepared data and the extracted model; reaction and gene lists are sorted so
that only the content is compared.
"""

import math

import cobra
import numpy as np

from raven_toolbox.init import TranscrData, ftinit, prep_init_model
from raven_toolbox.manipulation import add_reactions_from_equations, gpr_to_dnf
from raven_toolbox.tasks import Task

cobra.Configuration().processes = 1


def run(ctx):
    inputs = ctx["inputs"]
    cobra.Configuration().solver = str(inputs["python_solver"])
    results = {}
    for case in inputs["cases"]:
        results[str(case["name"])] = run_case(inputs, case)
    return results


def build_model(inputs, overrides):
    model = cobra.Model("testModel")
    model.add_metabolites(
        [
            cobra.Metabolite(str(m["id"]), name=str(m["name"]), compartment=str(m["compartment"]))
            for m in inputs["metabolites"]
        ]
    )
    rules = {str(r["id"]): str(r["rule"]) for r in inputs["reactions"]}
    if overrides:
        rules.update({str(k): str(v) for k, v in inputs["rule_overrides"].items()})
    add_reactions_from_equations(
        model,
        [
            {
                "id": str(r["id"]),
                "equation": str(r["equation"]),
                "bounds": (float(r["lb"]), float(r["ub"])),
                "gene_reaction_rule": rules[str(r["id"])],
            }
            for r in inputs["reactions"]
        ],
        mets_by="name",
        allow_new_mets=False,
    )
    model.objective = str(inputs["objective"])
    return model


def make_task(spec):
    return Task(
        id=str(spec["id"]),
        description=str(spec["description"]),
        inputs=[(str(i["metabolite"]), float(i["lb"]), math.inf) for i in as_list(spec["inputs"])],
        outputs=[(str(o["metabolite"]), float(o["lb"]), float(o["ub"])) for o in as_list(spec["outputs"])],
    )


def run_case(inputs, case):
    model = build_model(inputs, bool(case["overrides"]))
    model.solver = str(inputs["python_solver"])
    tasks = [make_task(inputs["task"])] if case["with_task"] else None
    prep = prep_init_model(
        model, tasks, ext_comp=str(inputs["ext_comp"]), spontaneous=[str(r) for r in as_list(case["spontaneous"])]
    )

    genes = sorted(g.id for g in model.genes)
    levels = np.array([[float(inputs["levels"][g])] for g in genes])
    transcr = TranscrData(genes, [str(inputs["tissue"])], [""], levels, float(inputs["threshold"]))
    out = ftinit(
        prep,
        str(inputs["tissue"]),
        None,
        None,
        transcr=transcr,
        series=str(case["series"]),
        remove_genes=bool(case["remove_genes"]),
    )

    masks = prep.masks
    objective = str(inputs["objective"])
    flux = 0.0
    if objective in out.reactions:
        out.objective = objective
        flux = float(out.slim_optimize())
    return {
        "prep": {
            "ref_reactions": sorted(r.id for r in prep.ref_model.reactions),
            "n_min_reactions": len(prep.min_model.reactions),
            "essential_reactions": sorted(prep.essential_rxns),
            "exchange": sorted(masks.exchange),
            "import": sorted(masks.import_rxns),
            "simple_transport": sorted(masks.simple_transport),
            "advanced_transport": sorted(masks.advanced_transport),
            "spontaneous": sorted(masks.spontaneous),
            "no_gpr": sorted(masks.no_gpr),
        },
        "model": {
            "reactions": sorted(r.id for r in out.reactions),
            "genes": sorted(g.id for g in out.genes),
            "rules": [
                {"reaction": r.id, "clauses": sorted(sorted(c) for c in gpr_to_dnf(r.gpr))}
                for r in sorted(out.reactions, key=lambda r: r.id)
            ],
            "has_objective": objective in out.reactions,
            "objective_flux": flux,
        },
    }


def as_list(value):
    return list(value) if value else []
