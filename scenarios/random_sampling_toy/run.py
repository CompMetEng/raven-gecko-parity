"""Python side of the randomSampling scenario: properties of the sample matrix that do not depend on the generator."""

import warnings

import cobra
import numpy as np

from raven_toolbox.analysis import random_sampling


def build_model(inputs):
    model = cobra.Model("sampling_toy")
    metabolites = {
        str(m["id"]): cobra.Metabolite(str(m["id"]), name=str(m["name"]), compartment=str(m["compartment"]))
        for m in inputs["mets"]
    }
    model.add_metabolites(list(metabolites.values()))
    for spec in inputs["reactions"]:
        reaction = cobra.Reaction(str(spec["id"]))
        model.add_reactions([reaction])
        reaction.add_metabolites({metabolites[str(e["met"])]: float(e["coefficient"]) for e in spec["stoichiometry"]})
        reaction.bounds = (float(spec["lb"]), float(spec["ub"]))
    return model


def run(ctx):
    inputs = ctx["inputs"]
    n = int(inputs["n_samples"])
    results = {}
    for method in inputs["methods"]:
        model = build_model(inputs)
        kwargs = {"method": str(method), "seed": 1}
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                samples = random_sampling(model, n, **kwargs).samples
        except Exception:
            results[method] = {"error": True}
            continue
        ids = [r.id for r in model.reactions]
        values = samples[ids].to_numpy(dtype=float).T
        lb = np.array([r.lower_bound for r in model.reactions])[:, None]
        ub = np.array([r.upper_bound for r in model.reactions])[:, None]
        s_matrix = cobra.util.create_stoichiometric_matrix(model)
        mean = values.mean(axis=1)
        i1, i2 = ids.index("R1"), ids.index("R2")
        results[method] = {
            "n_reactions": int(values.shape[0]),
            "n_samples": int(values.shape[1]),
            "steady_state": bool(np.abs(s_matrix @ values).max() < 1e-6),
            "max_residual": f"{np.abs(s_matrix @ values).max():.0e}",
            "steady_state_1e4": bool(np.abs(s_matrix @ values).max() < 1e-4),
            "within_bounds": bool(((values >= lb - 1e-6) & (values <= ub + 1e-6)).all()),
            "mean_near_centroid": bool(abs(mean[i1] - float(inputs["centroid"])) < float(inputs["band"])),
            "parallel_means_close": bool(abs(mean[i1] - mean[i2]) < float(inputs["band"])),
        }
    return results
