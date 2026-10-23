"""Python side of the curateModelFromTables edge-case scenario.

The model has one reaction and two metabolites with a formula and a charge. Each case curates it from a
set of tables and records the metabolites (formula and charge), the reactions (bounds and stoichiometry)
and the genes by name, so that identifiers, which each side numbers its own way, are not compared.
"""

from pathlib import Path

import cobra

from raven_toolbox.curation import batch_curate_from_tsv
from raven_toolbox.manipulation import gpr_to_dnf


def build_model():
    model = cobra.Model("toy")
    model.compartments = {"c": "c"}
    x = cobra.Metabolite("x_c", name="x", compartment="c", formula="CH4", charge=1)
    y = cobra.Metabolite("y_c", name="y", compartment="c", formula="C2H6", charge=2)
    rxn = cobra.Reaction("R0")
    rxn.name = "existing"
    model.add_reactions([rxn])
    rxn.add_metabolites({x: -1, y: 1})
    rxn.bounds = (0, 1000)
    return model


def run(ctx):
    inputs = ctx["inputs"]
    return {str(case["name"]): curate(inputs, case) for case in inputs["cases"]}


def curate(inputs, case):
    folder = Path(inputs["tables_dir"])

    def table(key):
        name = str(case[key])
        return None if name == "none" else folder / name

    model = build_model()
    try:
        batch_curate_from_tsv(
            model,
            mets_tsv=table("mets"),
            genes_tsv=table("genes"),
            rxns_tsv=table("rxns"),
            rxns_coeffs_tsv=table("coeffs"),
            met_id_prefix=str(inputs["met_prefix"]),
            rxn_id_prefix=str(inputs["rxn_prefix"]),
        )
    except Exception:
        return {"error": True}

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
    reactions = sorted(
        (
            {
                "name": str(r.name or ""),
                "lower_bound": float(r.lower_bound),
                "upper_bound": float(r.upper_bound),
                "clauses": sorted(sorted(clause) for clause in gpr_to_dnf(r.gpr)),
                "stoichiometry": sorted(
                    ({"metabolite": label(m), "coefficient": float(c)} for m, c in r.metabolites.items()),
                    key=lambda e: e["metabolite"],
                ),
            }
            for r in model.reactions
        ),
        key=lambda r: r["name"],
    )
    genes = sorted(
        ({"gene": g.id, "short_name": str(g.name or "")} for g in model.genes), key=lambda g: g["gene"]
    )
    return {
        "error": False,
        "metabolites": metabolites,
        "reactions": reactions,
        "genes": genes,
        "n_reactions": len(model.reactions),
        "n_metabolites": len(model.metabolites),
    }
