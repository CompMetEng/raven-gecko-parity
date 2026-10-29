"""Python side of the light Human-GEM flux variability scenario.

ec_fva on the tutorial's gecko-light ecModel with growth fixed at its optimum. See
scenario.yml.
"""

import math
import re
import tempfile
from pathlib import Path

import cobra

from geckopy import ModelAdapter, ec_fva, load_conventional_gem, load_ec_model


def num(x) -> dict:
    """A number as a record, so that NaN is a value both sides can state."""
    x = float(x)
    return {"is_nan": bool(math.isnan(x)), "value": 0.0 if math.isnan(x) else x}


def run(ctx):
    inputs = ctx["inputs"]
    cobra.Configuration().solver = str(inputs["python_solver"])
    folder = Path(inputs["adapter_python"])
    adapter = ModelAdapter.from_folder(folder)
    adapter.params.path = folder
    adapter.params.conv_gem = folder / "models" / "human-GEM.yml"
    conv = load_conventional_gem(adapter)

    # The tutorial ecModel has 52 ec reactions written as an `enzymes: !!omap` entry with no
    # content, which geckopy's reader rejects (see ec_model_load_light_humangem). Both sides
    # read a copy without them.
    text = (folder / "models" / str(inputs["model_file"])).read_text(encoding="utf-8")
    text = re.sub(r"^ *- enzymes: !!omap\n(?= {0,4}\S)", "", text, flags=re.MULTILINE)
    with tempfile.TemporaryDirectory(prefix="parity_light_fva_") as workdir:
        path = Path(workdir) / "repaired.yml"
        path.write_text(text, encoding="utf-8", newline="")
        model = load_ec_model(path, adapter)

    model.solver = str(inputs["python_solver"])
    model.tolerance = float(inputs["python_tolerance"])
    bio_id = adapter.params.bio_rxn
    model.objective = {model.reactions.get_by_id(bio_id): 1.0}
    model.objective_direction = "max"
    optimum = float(model.optimize().objective_value)
    bio = model.reactions.get_by_id(bio_id)
    bio.lower_bound = optimum
    bio.upper_bound = optimum

    df = ec_fva(model, conv, progress=False, n_proc=int(inputs["python_n_proc"]))
    stride = int(inputs["fva_stride"])
    rows = [
        {"reaction": str(r), "min_flux": num(row.min_flux), "max_flux": num(row.max_flux)}
        for i, (r, row) in enumerate(df.iterrows())
        if i % stride == 0
    ]
    return {"optimum": num(optimum), "n_reactions": len(df), "sample": rows}
