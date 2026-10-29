"""Write data/growth_fluxes.tsv: the growth-maximising flux vector of the tutorial ecYeastGEM.

The vector is the solution of the LP that run.py calls `growth` (glucose uptake open, growth as
the objective), solved by geckopy with Gurobi, written as `reaction<TAB>flux` sorted by reaction
id. Both sides of the scenario read this file, so that getEnzymeUsage, reportEnzymeUsage and
mapRxnsToConv are given the same fluxes: the LP has many optimal vertices and the two
implementations land on different ones (see scenario.yml).

    python regenerate_fluxes.py <GECKO>/tutorials/full_ecModel data/growth_fluxes.tsv
"""

import sys
from pathlib import Path

import cobra

from geckopy import ModelAdapter, load_ec_model

folder = Path(sys.argv[1])
cobra.Configuration().solver = "gurobi"
adapter = ModelAdapter.from_folder(folder)
adapter.params.path = folder
model = load_ec_model("ecYeastGEM.yml", adapter)
model.reactions.get_by_id("r_1714").lower_bound = -1000.0
model.objective = {model.reactions.get_by_id("r_4041"): 1.0}
model.objective_direction = "max"
fluxes = model.optimize().fluxes
with open(sys.argv[2], "w", encoding="utf-8", newline="") as handle:
    handle.write("reaction\tflux\n")
    for reaction in sorted(fluxes.index):
        handle.write(f"{reaction}\t{float(fluxes[reaction])!r}\n")
