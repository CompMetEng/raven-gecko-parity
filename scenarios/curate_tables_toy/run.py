"""Python side of the curateModelFromTables scenario.

The model has two empty reactions; the tables give them metabolites, rules and annotations.
Content is recorded by metabolite name and compartment, identifiers separately.
"""

import cobra

from raven_toolbox.curation import batch_curate_from_tsv
from raven_toolbox.manipulation import gpr_to_dnf


def run(ctx):
    inputs = ctx["inputs"]
    return {case["name"]: curate(inputs, case["reactions"]) for case in inputs["cases"]}


def curate(inputs, empty_reactions):
    model = cobra.Model("toy")
    model.compartments = {"c": "c"}
    for rxn_id in empty_reactions:
        model.add_reactions([cobra.Reaction(str(rxn_id))])

    try:
        batch_curate_from_tsv(
            model,
            mets_tsv=inputs["mets_table"],
            genes_tsv=inputs["genes_table"],
            rxns_tsv=inputs["rxns_table"],
            rxns_coeffs_tsv=inputs["coeffs_table"],
            met_id_prefix=str(inputs["met_prefix"]),
            rxn_id_prefix=str(inputs["rxn_prefix"]),
        )
    except (RuntimeError, ValueError) as error:
        return {"error": str(error)}

    def label(metabolite):
        return f"{metabolite.name}[{metabolite.compartment}]"

    metabolites = sorted(
        (
            {
                "metabolite": label(m),
                "formula": str(m.formula or ""),
                "charge": float(m.charge) if m.charge is not None else "NaN",
            }
            for m in model.metabolites
        ),
        key=lambda r: r["metabolite"],
    )
    reactions = []
    for r in model.reactions:
        stoichiometry = sorted(
            ({"metabolite": label(m), "coefficient": float(c)} for m, c in r.metabolites.items()),
            key=lambda e: e["metabolite"],
        )
        reactions.append(
            {
                "name": str(r.name or ""),
                "lower_bound": float(r.lower_bound),
                "upper_bound": float(r.upper_bound),
                "clauses": sorted(sorted(clause) for clause in gpr_to_dnf(r.gpr)),
                "stoichiometry": stoichiometry,
            }
        )
    reactions.sort(key=lambda r: (r["name"], str(r["stoichiometry"]), str(r["clauses"])))
    genes = sorted(
        ({"gene": g.id, "short_name": str(g.name or "")} for g in model.genes), key=lambda g: g["gene"]
    )
    return {
        "metabolites": metabolites,
        "reactions": reactions,
        "genes": genes,
        "n_reactions": len(model.reactions),
        "n_metabolites": len(model.metabolites),
        "reaction_ids": sorted(r.id for r in model.reactions),
        "metabolite_ids": sorted(m.id for m in model.metabolites),
    }
