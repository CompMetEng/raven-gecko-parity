"""Python side of the UniProt localisation scenario on a local export.

``load_uniprot`` is the parser behind ``fetch_uniprot_localization``; RAVEN's counterpart is
``parseScores`` with the ``uniprot`` predictor, which ``getUniProtScores`` uses after the download.
"""

import warnings

from raven_toolbox.localization import DEFAULT_COMPARTMENT_MAP, load_uniprot


def run(ctx):
    inputs = ctx["inputs"]
    results = {}
    for case in inputs["cases"]:
        compartment_map = DEFAULT_COMPARTMENT_MAP if case["map"] == "default" else inputs[case["map"]]
        id_column = None if case["id_column"] is None else str(inputs[case["id_column"]])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            scores = load_uniprot(
                inputs[case["file"]], compartment_map=dict(compartment_map), id_column=id_column
            ).df
        genes = [str(g) for g in scores.index]
        compartments = [str(c) for c in scores.columns]
        matrix = [[float(v) for v in row] for row in scores.to_numpy()]
        results[case["name"]] = {
            "genes": genes,
            "compartments": compartments,
            "scores": matrix,
            "by_gene": {
                gene: {comp: matrix[i][j] for j, comp in enumerate(compartments) if matrix[i][j] != 0}
                for i, gene in enumerate(genes)
            },
        }
    return results
