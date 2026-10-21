"""Python side of the YAML float-digits scenario.

Sets bounds on three reactions, writes the model and returns the lines of the
file that hold a lower_bound or upper_bound. The file goes to a temporary
directory: what the writer prints is the point, not an artefact to keep.
"""

import tempfile
from pathlib import Path

from raven_toolbox.io import read_yaml_model, write_yaml_model


def run(ctx):
    inputs = ctx["inputs"]
    model = read_yaml_model(inputs["model"])

    for rxn_id, bound in dict(inputs["bounds"]).items():
        rxn = model.reactions.get_by_id(str(rxn_id))
        # Widen first so neither assignment can pass through lb > ub.
        rxn.bounds = (float(bound["lb"]), float(bound["ub"]))

    with tempfile.TemporaryDirectory(prefix="parity_yaml_float_") as workdir:
        target = Path(workdir) / "model.yml"
        write_yaml_model(model, target)
        text = target.read_text(encoding="utf-8").replace("\r\n", "\n")

    lines = [
        line
        for line in text.split("\n")
        if "lower_bound" in line or "upper_bound" in line
    ]
    return {"n_bound_lines": len(lines), "bound_lines": lines}
