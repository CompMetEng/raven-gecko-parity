"""Python side of the parseHPA and parseHPArna scenario: records, independent of the table layout."""

import math

from raven_toolbox.init import parse_hpa, parse_hpa_rna


def run(ctx):
    inputs = ctx["inputs"]
    results = {}
    for case in inputs["protein"]:
        try:
            data = parse_hpa(f"{inputs['data_dir']}/{case['file']}")
        except (ValueError, KeyError):
            results[case["name"]] = {"error": True}
            continue
        df = data.df
        results[case["name"]] = {
            "genes": sorted(set(df["gene_id"])),
            "records": sorted(
                f"{r.gene_id}|{r.tissue}|{r.celltype}|{r.level}|{r.reliability}" for r in df.itertuples()
            ),
        }
    for case in inputs["rna"]:
        try:
            data = parse_hpa_rna(f"{inputs['data_dir']}/{case['file']}")
        except (ValueError, KeyError):
            results[case["name"]] = {"error": True}
            continue
        df = data.df
        results[case["name"]] = {
            "genes": sorted(set(df["gene_id"])),
            "tissues": sorted(set(df["tissue"])),
            "records": sorted(
                f"{r.gene_id}|{r.gene_name}|{r.tissue}|{float(r.tpm):.10g}"
                for r in df.itertuples()
                if not math.isnan(float(r.tpm))
            ),
        }
    return results
