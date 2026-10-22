"""Python side of the hitsFromSuppOrth scenario."""

import warnings

from raven_toolbox.reconstruction.homology import hits_from_supporth


def run(ctx):
    inputs = ctx["inputs"]
    results = {}
    for case in inputs["cases"]:
        kwargs = {
            "min_support": int(case["min_support"]),
            "require_tools": [str(t) for t in case["require_tools"]],
            "min_identity": case["min_identity"],
            "max_evalue": case["max_evalue"],
        }
        if case["use_gene_map"]:
            kwargs["protein_to_gene"] = inputs["gene_map"]
        mapping = case["species_to_model"]
        if mapping is not None:
            mapping = {str(k): (None if v is None else str(v)) for k, v in mapping.items()}
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                hits = hits_from_supporth(inputs["table"], str(case["organism"]), mapping, **kwargs)
            except ValueError:
                results[case["name"]] = {"error": True}
                continue
        pairs = sorted(
            f"{r.from_id}|{r.to_id}|{r.from_gene}|{r.to_gene}|{float(r.bitscore):.10g}|{float(r.evalue):g}|{float(r.identity):g}"
            for r in hits.itertuples()
        )
        forward = hits.iloc[: len(hits) // 2]
        info = sorted(
            f"{r.to_id}|{r.from_gene}|{r.to_gene}|{r.support}|{r.supported_by}" for r in forward.itertuples()
        )
        results[case["name"]] = {"pairs": pairs, "support": info}
    return results
