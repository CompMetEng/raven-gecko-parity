"""Section 5, Python half: smallYeast with an empty metabolite, reaction and gene name, written by raven-toolbox."""
import sys, warnings
from pathlib import Path
from raven_toolbox.io import read_yaml_model, write_yaml_model

out = Path(sys.argv[1]); out.mkdir(parents=True, exist_ok=True)
src = Path(sys.argv[2])
log = {}
for case in ("met", "rxn", "gene", "space"):
    m = read_yaml_model(src)
    if case == "met":
        m.metabolites[0].name = ""
    elif case == "rxn":
        m.reactions[0].name = ""
    elif case == "gene":
        m.genes[0].name = ""
    else:
        m.metabolites[0].name = " "
    f = out / f"empty_{case}.yml"
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        write_yaml_model(m, f)
    log[case] = {"file": f.name, "warnings": [str(x.message) for x in w]}
    print(case, log[case])
