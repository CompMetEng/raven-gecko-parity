"""Python side of the HPA branch of the scoreModel scenario."""

import pandas as pd

from raven_toolbox.init.hpa import HPAData, hpa_gene_scores
from raven_toolbox.io import read_yaml_model


def run(ctx):
    inputs = ctx["inputs"]
    model = read_yaml_model(inputs["model"])
    genes = sorted(gene.id for gene in model.genes)[:6]

    def build(rows):
        return HPAData(df=pd.DataFrame(
            [
                {
                    "gene_id": genes[int(row["gene"])],
                    "gene_name": genes[int(row["gene"])],
                    "tissue": str(row["tissue"]),
                    "celltype": str(row["celltype"]),
                    "level": str(row["level"]),
                    "reliability": "Supportive",
                }
                for row in rows
            ]
        ))

    results = {}
    for case in inputs["cases"]:
        celltype = str(case["celltype"]) or None
        rows = list(inputs["rows"]) + list(case.get("extra_rows", []))
        table = None
        if case.get("table"):
            table = {}
            for entry in case["table"]:
                table[str(entry["name"])] = float(entry["score"])
        try:
            scores = hpa_gene_scores(
                build(rows), str(case["tissue"]), celltype,
                level_scores=table, multiple_celltype=str(case["multiple"]),
            )
        except ValueError:
            results[case["name"]] = {"error": True}
            continue
        results[case["name"]] = [
            {
                "gene": gene,
                "has_score": gene in scores,
                "score": float(scores.get(gene, 0.0)),
            }
            for gene in genes
        ]
    return results
