"""Python side of the scoreModel scenario with zero and missing expression levels."""

import math

from raven_toolbox.init import gene_scores_from_expression, score_reactions_from_genes
from raven_toolbox.io import read_yaml_model


def run(ctx):
    inputs = ctx["inputs"]
    model = read_yaml_model(inputs["model"])
    genes = sorted(gene.id for gene in model.genes)
    results = {}
    for case in inputs["cases"]:
        expression = {
            gene: (math.nan if level is None else float(level))
            for gene, level in zip(genes[:6], case["levels"], strict=True)
        }
        gene_scores = gene_scores_from_expression(expression, float(case["threshold"]))
        reaction_scores = score_reactions_from_genes(
            model,
            gene_scores,
            isozyme_scoring=str(inputs["isozyme_scoring"]),
            complex_scoring=str(inputs["complex_scoring"]),
            no_gene_score=float(inputs["no_gene_score"]),
        )
        results[case["name"]] = {
            "gene_scores": [
                {
                    "gene": gene,
                    "has_score": gene in gene_scores and not math.isnan(gene_scores[gene]),
                    "score": 0.0
                    if gene not in gene_scores or math.isnan(gene_scores[gene])
                    else float(gene_scores[gene]),
                }
                for gene in genes[:6]
            ],
            "reaction_scores": [
                {"reaction": rxn, "score": float(reaction_scores[rxn])}
                for rxn in sorted(reaction_scores)
            ],
        }
    return results
