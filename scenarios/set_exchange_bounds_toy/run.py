"""Python side of the setExchangeBounds scenario on a hand-written model.

Each case applies set_exchange_bounds to a fresh copy of the model and records the bounds of every reaction and the
metabolites that were not used.
"""

import warnings

import cobra

from raven_toolbox.manipulation import set_exchange_bounds


def build_model(inputs, reactions):
    model = cobra.Model("exchange_toy")
    metabolites = {
        str(spec["id"]): cobra.Metabolite(str(spec["id"]), name=str(spec["name"]), compartment=str(spec["compartment"]))
        for spec in inputs["mets"]
    }
    model.add_metabolites(list(metabolites.values()))
    model.compartments = {str(k): str(v) for k, v in inputs["compartments"].items()}
    for spec in reactions:
        reaction = cobra.Reaction(str(spec["id"]))
        model.add_reactions([reaction])
        reaction.add_metabolites({metabolites[str(e["met"])]: float(e["coefficient"]) for e in spec["stoichiometry"]})
        reaction.bounds = (float(spec["lb"]), float(spec["ub"]))
    return model


def run(ctx):
    inputs = ctx["inputs"]
    results = {}
    for case in inputs["cases"]:
        model = build_model(inputs, case["reactions"])
        kwargs = {"close_others": bool(case["close_others"]), "media_only": bool(case["media_only"])}
        mets = None if case["mets"] is None else [str(m) for m in case["mets"]]
        lb, ub = case["lb"], case["ub"]
        lb = [float(x) for x in lb] if isinstance(lb, list) else (None if lb is None else float(lb))
        ub = [float(x) for x in ub] if isinstance(ub, list) else (None if ub is None else float(ub))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                unused = set_exchange_bounds(model, mets, lb, ub, **kwargs)
            except (ValueError, KeyError):
                results[case["name"]] = {"error": True}
                continue
        results[case["name"]] = {
            "bounds": [
                {"reaction": r.id, "lb": float(r.lower_bound), "ub": float(r.upper_bound)}
                for r in sorted(model.reactions, key=lambda r: r.id)
            ],
            "unused": sorted(str(m) if isinstance(m, str) else m.id for m in unused),
        }
    return results
