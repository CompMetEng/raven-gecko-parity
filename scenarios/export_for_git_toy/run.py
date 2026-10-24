"""Python side of the exportForGit scenario.

A four-reaction model is exported twice: as txt only, to compare the table, and with the
default formats, to compare which files are written. Equations are recorded apart from the
rest of each row, because the two sides spell them differently.
"""

import tempfile
from pathlib import Path

import cobra

from raven_toolbox.io import export_for_git
from raven_toolbox.io.metadata import set_model_metadata


def build(reactions):
    model = cobra.Model("exportToy")
    model.compartments = {"c": "cytosol"}
    metabolites = {}
    for spec in reactions:
        reaction = cobra.Reaction(str(spec["id"]))
        model.add_reactions([reaction])
        stoichiometry = {}
        for side, sign in (("substrates", -1), ("products", 1)):
            for name in spec[side] or []:
                if name not in metabolites:
                    metabolites[name] = cobra.Metabolite(
                        f"{name}_c", name=str(name), compartment="c", formula="C"
                    )
                stoichiometry[metabolites[name]] = sign
        reaction.add_metabolites(stoichiometry)
        reaction.bounds = (float(spec["lb"]), float(spec["ub"]))
        reaction.objective_coefficient = float(spec["objective"])
        reaction.name = str(spec["id"])
        if spec["rule"]:
            reaction.gene_reaction_rule = str(spec["rule"])
    set_model_metadata(model, version="1.0.0", date="2026-10-22")
    return model


def flag_checks(model, flag_cases):
    """Whether a call with one flag value is refused as a type error; text is refused, a bool is not."""
    checks = {}
    with tempfile.TemporaryDirectory(prefix="parity_flags_") as workdir:
        for case in flag_cases:
            try:
                export_for_git(model, Path(workdir) / str(case["name"]), formats=("txt",),
                               **{str(case["flag"]): case["value"]})
            except TypeError:
                checks[str(case["name"])] = {"rejected": True}
            else:
                checks[str(case["name"])] = {"rejected": False}
    return checks


def run(ctx):
    model = build(ctx["inputs"]["reactions"])
    with tempfile.TemporaryDirectory(prefix="parity_export_") as workdir:
        root = export_for_git(model, Path(workdir) / "txt", formats=("txt",))
        lines = (root / "txt" / "model.txt").read_text(encoding="utf-8").splitlines()
        dependencies = (root / "dependencies.txt").read_text(encoding="utf-8").splitlines()

        default_root = export_for_git(model, Path(workdir) / "all")
        written = sorted(
            f"{p.parent.name}/{p.name}" for p in default_root.rglob("*") if p.is_file()
        )

    header = lines[0].split("\t")
    rows = []
    for line in lines[1:]:
        fields = line.split("\t")
        rows.append(
            {
                "reaction": fields[0],
                "equation": fields[1],
                "rule": fields[2],
                "lower_bound": fields[3],
                "upper_bound": fields[4],
                "objective": fields[5],
            }
        )
    rows.sort(key=lambda r: r["reaction"])
    return {
        "header": header,
        "rows": rows,
        "libsbml_recorded": any(line.split("\t")[0] == "libSBML" for line in dependencies),
        "files_written": written,
        "flag_checks": flag_checks(model, ctx["inputs"]["flag_cases"]),
    }
