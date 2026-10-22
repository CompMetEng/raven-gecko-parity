"""Python side of the HPA branch of the scoreModel scenario."""

import pandas as pd

from raven_toolbox.init.hpa import HPAData, hpa_gene_scores
from raven_toolbox.io import read_yaml_model


def run(ctx):
    inputs = ctx["inputs"]
    model = read_yaml_model(inputs["model"])
    genes = sorted(gene.id for gene in model.genes)[:6]
    df = pd.DataFrame(
        [
            {
                "gene_id": genes[int(row["gene"])],
                "gene_name": genes[int(row["gene"])],
                "tissue": str(row["tissue"]),
                "celltype": str(row["celltype"]),
                "level": str(row["level"]),
                "reliability": "Supportive",
            }
            for row in inputs["rows"]
        ]
    )
    hpa = HPAData(df=df)
    results = {}
    for case in inputs["cases"]:
        celltype = str(case["celltype"]) or None
        scores = hpa_gene_scores(
            hpa, str(case["tissue"]), celltype, multiple_celltype=str(case["multiple"])
        )
        results[case["name"]] = [
            {
                "gene": gene,
                "has_score": gene in scores,
                "score": float(scores.get(gene, 0.0)),
            }
            for gene in genes
        ]
    return results
