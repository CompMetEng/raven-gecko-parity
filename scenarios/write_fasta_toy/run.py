"""Python side of the writeFasta scenario."""

import tempfile
from pathlib import Path

from raven_toolbox.localization import write_fasta


def run(ctx):
    inputs = ctx["inputs"]
    records = {str(r["header"]): str(r["sequence"]) for r in inputs["records"]}
    results = {}
    with tempfile.TemporaryDirectory() as tmp:
        for case in inputs["cases"]:
            path = Path(tmp) / f"{case['name']}.fa"
            write_fasta(records, path, wrap=int(case["width"]))
            results[case["name"]] = {"text": path.read_text(encoding="utf-8")}
    return results
