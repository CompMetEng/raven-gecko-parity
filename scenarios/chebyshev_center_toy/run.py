"""Python side of the sampleChebyshevCenter scenario."""

import numpy as np

from raven_toolbox.analysis.flux_sampling import _chebyshev_center


def run(ctx):
    results = {}
    for case in ctx["inputs"]["cases"]:
        A = np.array(case["A"], dtype=float)
        b = np.array(case["b"], dtype=float)
        x0 = _chebyshev_center(A, b)
        radius = float(np.min((b - A @ x0) / np.linalg.norm(A, axis=1)))
        record = {"radius": radius}
        if case["unique_centre"]:
            record["centre"] = [float(v) for v in x0]
        results[case["name"]] = record
    return results
