"""Rewrite data/ecModel.python.yml and data/ecModel.matlab.yml from the results of
ec_model_io_ectestgem.

Usage: python regenerate_fixtures.py <python result json> <matlab result json>

Each file is the ``written.lines`` checkpoint of one side of ec_model_io_ectestgem: what
save_ec_model (Python) and saveEcModel (MATLAB) write for the freshly built ecTestGEM
ecModel. Run it after a change to either writer; the scenario refuses to run against a
fixture its own writer no longer reproduces.
"""

import json
import sys
from pathlib import Path

data = Path(__file__).parent / "data"
for name, source in (("python", sys.argv[1]), ("matlab", sys.argv[2])):
    lines = json.loads(Path(source).read_text(encoding="utf-8"))["results"]["written"]["lines"]
    with open(data / f"ecModel.{name}.yml", "w", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(lines) + "\n")
