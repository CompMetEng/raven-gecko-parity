"""Python side of the compartment-merging scenario.

``merge_compartments`` is handed a name-based grouping key rather than left on
its default id-suffix one, because RAVEN groups by ``metNames``. Merged
metabolites are compared by id, name and ChEBI annotation. smallYeast is merged
first, then each toy case, a model built from the scenario inputs. See
scenario.yml.

Reactions removed by the merge are reported as the difference between the model
before and after, not as the lists each function returns: RAVEN reports
``deletedRxns`` for the deleteRxnsWithOneMet mode only, and removes the
reactions that cancelled to nothing separately without reporting them, so the
returned lists are not the same quantity on the two sides. The difference is.
"""

import cobra

from raven_toolbox.io import read_yaml_model
from raven_toolbox.manipulation import merge_compartments


def run(ctx):
    inputs = ctx["inputs"]
    flag = bool(inputs["delete_single_metabolite_reactions"])
    results = _merge(inputs, read_yaml_model(inputs["model"]), flag)
    results["toy_cases"] = [
        {
            "name": str(case["name"]),
            **_merge(inputs, _toy_model(case), bool(case.get("delete_single_metabolite_reactions", flag))),
        }
        for case in inputs["toy_cases"]
    ]
    return results


def _toy_model(case):
    model = cobra.Model(str(case["name"]))
    mets = {}
    for entry in case["metabolites"]:
        met = cobra.Metabolite(str(entry["id"]), name=str(entry["name"]), compartment=str(entry["comp"]))
        if entry["chebi"]:
            met.annotation = {"chebi": str(entry["chebi"])}
        mets[met.id] = met
    model.add_metabolites(list(mets.values()))
    reactions = []
    for entry in case["reactions"]:
        rxn = cobra.Reaction(str(entry["id"]), lower_bound=float(entry["lb"]), upper_bound=float(entry["ub"]))
        pairs = zip(entry["mets"], entry["coefficients"], strict=True)
        rxn.add_metabolites({mets[str(m)]: float(c) for m, c in pairs})
        reactions.append(rxn)
    model.add_reactions(reactions)
    return model


def _merge(inputs, model, delete_single_metabolite_reactions):
    before_rxns = {rxn.id for rxn in model.reactions}
    before_species = {met.name for met in model.metabolites}

    merged, _deleted, deduplicated = merge_compartments(
        model,
        merged_id=str(inputs["merged_compartment"]),
        base_metabolite=lambda met: met.name,
        drop_single_metabolite_reactions=delete_single_metabolite_reactions,
        deduplicate_reactions=True,
    )

    after_rxns = {rxn.id for rxn in merged.reactions}

    return {
        "n_reactions_before": len(before_rxns),
        "n_reactions_after": len(after_rxns),
        "n_species_before": len(before_species),
        "n_metabolites_before": len(model.metabolites),
        "n_metabolites_after": len(merged.metabolites),
        "compartments": sorted(merged.compartments),
        "reactions": sorted(after_rxns),
        "removed_reactions": sorted(before_rxns - after_rxns),
        "deduplicated_reactions": sorted(deduplicated),
        "metabolites": sorted(met.name for met in merged.metabolites),
        "metabolite_ids": [
            {"id": met.id, "name": str(met.name), "chebi": _chebi(met)}
            for met in sorted(merged.metabolites, key=lambda m: m.id)
        ],
        "bounds": [
            {
                "reaction": rxn.id,
                "lower_bound": float(rxn.lower_bound),
                "upper_bound": float(rxn.upper_bound),
            }
            for rxn in sorted(merged.reactions, key=lambda r: r.id)
        ],
        # This is where the coefficient summing shows up, the part of the
        # merge with any arithmetic in it.
        "stoichiometry": _stoichiometry(merged),
    }


def _chebi(met):
    value = met.annotation.get("chebi", "")
    return ";".join(str(v) for v in value) if isinstance(value, list) else str(value)


def _stoichiometry(model):
    entries = [
        {"reaction": rxn.id, "metabolite": met.id, "coefficient": float(coeff)}
        for rxn in model.reactions
        for met, coeff in rxn.metabolites.items()
    ]
    entries.sort(key=lambda e: (e["reaction"], e["metabolite"]))
    return entries
