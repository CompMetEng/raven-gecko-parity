"""Python side of the curateModelFromTables table layout scenario.

The model has two metabolites, x and y in compartment c, and the reactions R0 (x to y) and R0b (y to x). Each case
curates it from a set of table files and records the metabolites (name, compartment, formula, charge) and the
reactions (name, bounds and stoichiometry by metabolite name) so that identifiers, which each side numbers its own
way, are not compared. The cases cover the layout of the files (a trailing tab, an extra cell, spaces before a cell,
quotes, blank lines), rows without content, repeated coefficient rows, reaction indices written in different
notations, reactions that only one of the two reaction tables has, and cells that are refused.
"""

import tempfile
from pathlib import Path

import cobra

from raven_toolbox.curation import batch_curate_from_tsv


def build_model(gene_ids=("G1",)):
    model = cobra.Model("toy")
    model.compartments = {"c": "c"}
    x = cobra.Metabolite("x_c", name="x", compartment="c", formula="CH4", charge=1)
    y = cobra.Metabolite("y_c", name="y", compartment="c", formula="C2H6", charge=2)
    rxn = cobra.Reaction("R0")
    rxn.name = "existing"
    model.add_reactions([rxn])
    rxn.add_metabolites({x: -1, y: 1})
    rxn.bounds = (0, 1000)
    second = cobra.Reaction("R0b")
    second.name = "existing2"
    model.add_reactions([second])
    second.add_metabolites({x: 1, y: -1})
    second.bounds = (0, 1000)
    for gene_id in gene_ids:
        model.genes.append(cobra.Gene(str(gene_id)))
    return model


def error_class(message):
    """Kind of table error, from the words both sides use in their messages."""
    for words, kind in (
        (("needs its own index",), "repeated_rxn_index"),
        (("finite whole number",), "charge"),
        (("0, 1, true or false",), "rev"),
        (("rxnIdx is empty",), "rxn_index_missing"),
        (("at least 1",), "rxn_index"),
        (("no matching row", "not matched"), "reaction_not_matched"),
        (("has no compartment", "must match model.comps"), "compartment"),
        (("cells but the header has",), "extra_cell"),
        (("not a valid gene rule",), "gene_rule"),
        (("more than once",), "repeated_gene"),
        (("Not all genes", "not in the model or in genes_df"), "unknown_gene"),
        (("coefficient",), "coefficient"),
        (("is not a number",), "not_a_number"),
    ):
        if any(word in message for word in words):
            return kind
    return "other"


def run(ctx):
    inputs = ctx["inputs"]
    results = {str(case["name"]): curate(inputs, case) for case in inputs["cases"]}
    # One new reaction per gene rule, written to a table of its own; the model holds every gene the rules name.
    with tempfile.TemporaryDirectory() as scratch:
        for number, rule in enumerate(inputs["gene_rules"], start=1):
            rxns = Path(scratch) / f"rxns_{number}.tsv"
            rxns.write_text(f"rxnIdx\trxnNames\tgrRules\n1\tfresh\t{rule}\n", encoding="utf-8")
            case = {"name": f"gene_rule_{number:02d}", "mets": "none", "rxns": str(rxns), "coeffs": "coeffs_fresh.tsv"}
            results[case["name"]] = curate(inputs, case)
    return results


def curate(inputs, case):
    folder = Path(inputs["tables_dir"])

    def table(key):
        name = str(case.get(key, "none"))
        if name == "none":
            return None
        return Path(name) if Path(name).is_absolute() else folder / name

    model = build_model(inputs["model_genes"])
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
    except Exception as exc:
        return {"error": True, "error_class": "" if case.get("ignore_error_class") else error_class(str(exc))}

    def bound(value):
        return "Inf" if value == float("inf") else "-Inf" if value == float("-inf") else float(value)

    metabolites = sorted(
        (
            {
                "metabolite": f"{m.name}[{m.compartment}]",
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
                "lb": bound(r.lower_bound),
                "ub": bound(r.upper_bound),
                "genes": sorted(g.id for g in r.genes),
                "stoichiometry": sorted(
                    ({"metabolite": f"{m.name}[{m.compartment}]", "coefficient": float(c)}
                     for m, c in r.metabolites.items()),
                    key=lambda e: e["metabolite"],
                ),
            }
            for r in model.reactions
        ),
        key=lambda r: r["name"],
    )
    return {
        "error": False,
        "metabolites": metabolites,
        "reactions": reactions,
        "model_genes": sorted(g.id for g in model.genes),
    }
