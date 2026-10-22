"""Python side of the edge-case model scenario.

The model is built from scenario.yml, written to YAML and read back, split into irreversible
reactions, and one gene is removed with the default arguments. A second model, with two
metabolites of one name in one compartment, only records whether it can be written.
"""

import tempfile
from pathlib import Path

import cobra

from raven_toolbox.io import read_yaml_model, write_yaml_model
from raven_toolbox.manipulation import convert_to_irreversible, gpr_to_dnf, remove_genes


def build(inputs, duplicate=False):
    model = cobra.Model("edge", name="edge")
    model.compartments = {"c": "cytosol"}
    metabolites = {}

    def metabolite(name):
        if name not in metabolites:
            display = inputs["duplicate_name"] if duplicate and name in ("X1", "X2") else name
            metabolites[name] = cobra.Metabolite(
                f"{name}_c", name=str(display), compartment="c", formula="C"
            )
        return metabolites[name]

    for spec in inputs["reactions"]:
        reaction = cobra.Reaction(str(spec["id"]))
        model.add_reactions([reaction])
        stoichiometry = {}
        for side, sign in (("substrates", -1), ("products", 1)):
            for name in spec[side] or []:
                stoichiometry[metabolite(str(name))] = sign
        reaction.add_metabolites(stoichiometry)
        reaction.bounds = (float(spec["lb"]), float(spec["ub"]))
        reaction.name = str(spec["name"])
        reaction.objective_coefficient = float(spec["objective"])
        if spec["rule"]:
            reaction.gene_reaction_rule = str(spec["rule"])
    model.add_metabolites([metabolite(str(inputs["orphan_metabolite"]))])
    model.genes.append(cobra.Gene(str(inputs["orphan_gene"])))
    return model


def summary(model):
    reactions = []
    for r in sorted(model.reactions, key=lambda r: r.id):
        reactions.append(
            {
                "reaction": r.id,
                "lower_bound": float(r.lower_bound),
                "upper_bound": float(r.upper_bound),
                "clauses": sorted(sorted(c) for c in gpr_to_dnf(r.gpr)),
                "objective": float(r.objective_coefficient),
            }
        )
    return {
        "reactions": reactions,
        "n_metabolites": len(model.metabolites),
        "genes": sorted(g.id for g in model.genes),
    }


def run(ctx):
    inputs = ctx["inputs"]
    model = build(inputs)
    with tempfile.TemporaryDirectory(prefix="parity_edge_") as workdir:
        path = Path(workdir) / "edge.yml"
        try:
            write_yaml_model(model, path)
            lines = path.read_text(encoding="utf-8").splitlines()
            reread = summary(read_yaml_model(path))
        except (RuntimeError, ValueError) as error:
            lines, reread = [], {"error": str(error)}

    with tempfile.TemporaryDirectory(prefix="parity_edge_") as workdir:
        try:
            write_yaml_model(build(inputs, duplicate=True), Path(workdir) / "dup.yml")
            duplicate_names_written = True
        except (RuntimeError, ValueError):
            duplicate_names_written = False

    irreversible = build(inputs)
    convert_to_irreversible(irreversible)
    removed = build(inputs)
    remove_genes(removed, [str(inputs["removed_gene"])])

    return {
        "duplicate_names_written": duplicate_names_written,
        "yaml_lines": [line for line in lines if not line.lstrip().startswith(("date:", "- date:"))],
        "reread": reread,
        "irreversible": summary(irreversible),
        "removed_gene": summary(removed),
    }
