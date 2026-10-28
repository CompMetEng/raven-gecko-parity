"""Python side of the sampleMaxVolEllipse scenario."""

import numpy as np

from raven_toolbox.analysis.flux_sampling import max_volume_ellipsoid


def run(ctx):
    results = {}
    for case in ctx["inputs"]["cases"]:
        a = np.array(case["A"], dtype=float)
        b = np.array(case["b"], dtype=float)
        x0 = None if case["x0"] is None else np.array(case["x0"], dtype=float)
        x, e, converged = max_volume_ellipsoid(a, b, x0)
        results[case["name"]] = {
            "centre": [float(v) for v in x],
            "E": [[float(v) for v in row] for row in e],
            "E_Et": [[float(v) for v in row] for row in e @ e.T],
            "converged": bool(converged),
        }
    return results
