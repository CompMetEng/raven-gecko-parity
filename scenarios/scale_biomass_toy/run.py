"""Python side of the biomass rescaling and fraction scenario on hand-written models.

Each case applies one operation to a fresh model: rescale a pseudoreaction by a factor, scale a component to a
mass fraction (optionally balancing a second component), or compute the fractions. A case records the
coefficients of every reaction, and the fractions where the operation produces them; a refused operation is
recorded as an error.
"""

import cobra

from raven_toolbox.biomass import BiomassComponent, BiomassConfig, rescale_pseudoreaction, scale_biomass, sum_biomass


def build_model(inputs, case):
    model = cobra.Model("biomass_toy")
    formula_overrides = {str(o["met"]): str(o["formula"]) for o in case["formula_overrides"]}
    charge_overrides = {str(o["met"]): o["charge"] for o in case["charge_overrides"]}
    metabolites = {}
    for spec in inputs["models"][case["model"]]:
        met_id = str(spec["id"])
        charge = None if case["omit_charges"] else spec["charge"]
        if met_id in charge_overrides:
            charge = charge_overrides[met_id]
        metabolites[met_id] = cobra.Metabolite(
            met_id,
            name=str(spec["name"]),
            formula=formula_overrides.get(met_id, str(spec["formula"])),
            charge=None if charge is None else int(charge),
            compartment=str(spec["compartment"]),
        )
    model.add_metabolites(list(metabolites.values()))
    for spec in inputs["reactions"]:
        reaction = cobra.Reaction(str(spec["id"]), name=str(spec["name"]))
        model.add_reactions([reaction])
        stoichiometry = {metabolites[str(e["met"])]: float(e["coefficient"]) for e in spec["stoichiometry"]}
        if case["proton_in_protein"] and spec["id"] == "P1":
            stoichiometry[metabolites["h_c"]] = -0.1
        reaction.add_metabolites(stoichiometry)
        reaction.bounds = (float(spec["lb"]), float(spec["ub"]))
    return model


def build_config(inputs, case):
    components = [
        BiomassComponent(str(c["name"]), str(c["pseudoreaction_name"]), str(c["mass_strategy"]))
        for c in inputs["configs"][case["config"]]
    ]
    return BiomassConfig("P1", str(case["proton_met"]), components)


def coefficients(model):
    return sorted(
        (
            {"reaction": reaction.id, "met": met.id, "coefficient": float(coef)}
            for reaction in model.reactions
            for met, coef in reaction.metabolites.items()
            if coef != 0
        ),
        key=lambda e: (e["reaction"], e["met"]),
    )


def run_case(inputs, case):
    model = build_model(inputs, case)
    config = build_config(inputs, case)
    try:
        if case["op"] == "rescale":
            rescale_pseudoreaction(model, config, str(case["component"]), float(case["factor"]))
            return {"coefficients": coefficients(model)}
        if case["op"] == "scale":
            kwargs = {}
            if case["balance_out"] is not None:
                kwargs["balance_out"] = str(case["balance_out"])
            if case["target_total"] is not None:
                kwargs["target_total"] = float(case["target_total"])
            scale_biomass(model, config, str(case["component"]), float(case["new_value"]), **kwargs)
            return {"coefficients": coefficients(model), "fractions": sum_biomass(model, config)}
        return {"fractions": sum_biomass(model, config)}
    except (KeyError, ValueError):
        return {"error": True}


def run(ctx):
    inputs = ctx["inputs"]
    return {case["name"]: run_case(inputs, case) for case in inputs["cases"]}
