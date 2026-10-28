"""Python side of the curateModelFromTables annotation scenario.

The model has two metabolites with a KEGG compound annotation, two genes with a short name and an NCBI gene
annotation, and two reactions with a KEGG reaction annotation and a confidence score. Each case curates it from
a set of tables and records annotations, charges and scores by name, so that identifiers, which each side
numbers its own way, are not compared.
"""

from pathlib import Path

import cobra
from cobra.core.gene import Gene

from raven_toolbox.curation import batch_curate_from_tsv


def build_model():
    model = cobra.Model("toy")
    model.compartments = {"c": "c"}
    x = cobra.Metabolite("x_c", name="x", compartment="c", formula="CH4", charge=1)
    y = cobra.Metabolite("y_c", name="y", compartment="c", formula="C2H6", charge=2)
    x.annotation["kegg.compound"] = "C1"
    y.annotation["kegg.compound"] = "C2"
    rxn = cobra.Reaction("R0")
    rxn.name = "existing"
    model.add_reactions([rxn])
    rxn.add_metabolites({x: -1, y: 1})
    rxn.bounds = (0, 1000)
    rxn.annotation["kegg.reaction"] = "R1"
    rxn.notes["confidence_score"] = 2
    second = cobra.Reaction("R0b")
    second.name = "existing2"
    model.add_reactions([second])
    second.add_metabolites({x: 1, y: -1})
    second.bounds = (0, 1000)
    second.annotation["kegg.reaction"] = "R2"
    second.notes["confidence_score"] = 4
    for gid, name, ncbi in (("G1", "one", "11"), ("G2", "two", "22")):
        gene = Gene(gid, name=name)
        gene.annotation["ncbigene"] = ncbi
        model.genes.append(gene)
    return model


def error_class(message):
    """Kind of table error, from the words both sides use in their messages."""
    for words, kind in (("finite whole number", "charge"), ("0, 1, true or false", "rev"),
                        ("at least 1", "rxn_index"), ("coefficient", "coefficient"),
                        ("is not a number", "not_a_number")):
        if words in message:
            return kind
    return "other"


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
            overwrite_with_empty=bool(case.get("overwrite_with_empty", True)),
        )
    except Exception as exc:
        return {"error": True, "error_class": error_class(str(exc))}

    def bound(value):
        return "Inf" if value == float("inf") else "-Inf" if value == float("-inf") else float(value)

    def annotation(entity, namespace):
        return str(entity.annotation.get(namespace, ""))

    metabolites = sorted(
        (
            {
                "metabolite": f"{m.name}[{m.compartment}]",
                "formula": str(m.formula or ""),
                "charge": float(m.charge) if m.charge is not None else "NaN",
                "kegg_compound": annotation(m, "kegg.compound"),
            }
            for m in model.metabolites
        ),
        key=lambda r: r["metabolite"],
    )
    reactions = sorted(
        (
            {
                "name": str(r.name or ""),
                "kegg_reaction": annotation(r, "kegg.reaction"),
                "metacyc": annotation(r, "metacyc"),
                "confidence_score": float(r.notes["confidence_score"]) if "confidence_score" in r.notes else "NaN",
                "lb": bound(r.lower_bound),
                "ub": bound(r.upper_bound),
            }
            for r in model.reactions
        ),
        key=lambda r: r["name"],
    )
    genes = sorted(
        (
            {"gene": g.id, "short_name": str(g.name or ""), "ncbigene": annotation(g, "ncbigene")}
            for g in model.genes
        ),
        key=lambda g: g["gene"],
    )
    return {"error": False, "metabolites": metabolites, "reactions": reactions, "genes": genes}
