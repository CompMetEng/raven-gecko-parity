"""Python side of the ecModel cross-load scenario.

``load_ec_model`` reads the file ``save_ec_model`` wrote and the file ``saveEcModel`` wrote
(both vendored under data/), and each loaded model is summarised. See scenario.yml for why
the files are fixtures and how they are kept current. ``model.ec.concs`` is left out of the
summary for the reason given in ec_model_io_ectestgem.
"""

import shutil
import tempfile
from pathlib import Path

from geckopy import ModelAdapter, load_conventional_gem, load_ec_model, make_ec_model, save_ec_model

_DATE_PREFIX = "  - date:"


def _bound(value: float) -> tuple[float, str]:
    """An infinite bound as a class and a zero --- see ec_model_expansion_ectestgem."""
    value = float(value)
    if value == float("inf"):
        return 0.0, "+inf"
    if value == float("-inf"):
        return 0.0, "-inf"
    return value, "finite"


def _adapter(inputs: dict) -> ModelAdapter:
    adapter = ModelAdapter.from_folder(inputs["adapter_python"])
    fixture = Path(inputs["fixture_dir"])
    adapter.params.path = fixture
    adapter.params.conv_gem = fixture / "models" / "testModel.xml"
    return adapter


def _names(value) -> list[str]:
    """An annotation or subsystem value as a sorted list of strings, whatever its shape."""
    if value is None or value == "":
        return []
    if isinstance(value, str):
        return [value]
    return sorted(str(v) for v in value if str(v) != "")


def _sbo(entity) -> str:
    return ";".join(_names(entity.annotation.get("sbo")))


def _ec_codes(rxn) -> str:
    value = rxn.annotation.get("ec-code")
    if value is None:
        return ""
    return ";".join([value] if isinstance(value, str) else [str(v) for v in value])


def _gene_associations(model) -> list[dict]:
    entries = [{"reaction": rxn.id, "gene": gene.id} for rxn in model.reactions for gene in rxn.genes]
    return sorted(entries, key=lambda e: (e["reaction"], e["gene"]))


def _stoichiometry(model) -> list[dict]:
    entries = [
        {"reaction": rxn.id, "metabolite": met.id, "coefficient": float(coefficient)}
        for rxn in model.reactions
        for met, coefficient in rxn.metabolites.items()
    ]
    return sorted(entries, key=lambda e: (e["reaction"], e["metabolite"]))


def _reactions(model) -> list[dict]:
    records = []
    for rxn in sorted(model.reactions, key=lambda r: r.id):
        lower, lower_kind = _bound(rxn.lower_bound)
        upper, upper_kind = _bound(rxn.upper_bound)
        records.append(
            {
                "id": rxn.id,
                "lower_bound": lower,
                "lower_kind": lower_kind,
                "upper_bound": upper,
                "upper_kind": upper_kind,
                "objective_coefficient": float(rxn.objective_coefficient),
                "subsystem": _names(rxn.subsystem),
                "ec_codes": _ec_codes(rxn),
                "sbo": _sbo(rxn),
            }
        )
    return records


def _metabolites(model) -> list[dict]:
    records = []
    for met in sorted(model.metabolites, key=lambda m: m.id):
        records.append(
            {
                "id": met.id,
                "name": str(met.name or ""),
                "compartment": str(met.compartment or ""),
                "charge_set": met.charge is not None,
                "charge": float(met.charge or 0.0),
                "sbo": _sbo(met),
            }
        )
    return records


def _genes(model) -> list[dict]:
    return [
        {"id": gene.id, "name": str(gene.name or ""), "sbo": _sbo(gene)}
        for gene in sorted(model.genes, key=lambda g: g.id)
    ]


def _ec_data(ec) -> dict:
    matrix = ec.rxn_enz_mat
    if hasattr(matrix, "toarray"):
        matrix = matrix.toarray()
    coupling = [
        {"reaction": ec.rxns[i], "enzyme": ec.enzymes[j], "coefficient": float(matrix[i][j])}
        for i in range(len(ec.rxns))
        for j in range(len(ec.enzymes))
        if matrix[i][j] != 0
    ]
    return {
        # Not sorted: the expansion order is the result, as in ec_model_expansion_ectestgem.
        "rxns": [str(r) for r in ec.rxns],
        "genes": [str(g) for g in ec.genes],
        "enzymes": [str(e) for e in ec.enzymes],
        "mw": [float(m) for m in ec.mw],
        "sequence": [str(s) for s in ec.sequence],
        "eccodes": [str(c) if c else "" for c in ec.eccodes],
        "kcat": [float(k) for k in ec.kcat],
        "coupling": sorted(coupling, key=lambda e: (e["reaction"], e["enzyme"])),
    }


def _summary(model) -> dict:
    meta = model.notes.get("metaData") or {}
    return {
        "loaded": True,
        "model_id": str(model.id or ""),
        "description": str(meta.get("description") or ""),
        "geckopy_version": str(meta.get("geckopy_version") or ""),
        "n_reactions": len(model.reactions),
        "n_metabolites": len(model.metabolites),
        "n_genes": len(model.genes),
        "reactions": _reactions(model),
        "metabolites": _metabolites(model),
        "genes": _genes(model),
        "gene_associations": _gene_associations(model),
        "stoichiometry": _stoichiometry(model),
        "gecko_light": bool(model.ec.gecko_light),
        "ec": _ec_data(model.ec),
    }


def _load(path: Path) -> dict:
    try:
        return _summary(load_ec_model(path))
    except Exception as error:  # a reader that refuses the file is the result
        return {"loaded": False, "error": type(error).__name__}


def _without_date(path: Path) -> list[str]:
    lines = path.read_text(encoding="utf-8").replace("\r\n", "\n").split("\n")
    return [line for line in lines if not line.startswith(_DATE_PREFIX)]


def run(ctx):
    inputs = ctx["inputs"]
    data = Path(inputs["scenario_root"]) / "data"
    adapter = _adapter(inputs)
    ec_model = make_ec_model(load_conventional_gem(adapter), adapter, gecko_light=False)

    with tempfile.TemporaryDirectory(prefix="parity_ecmodel_cross_") as workdir:
        own = Path(workdir) / "ecModel.yml"
        save_ec_model(ec_model, own)
        if _without_date(own) != _without_date(data / "ecModel.python.yml"):
            raise RuntimeError(
                "save_ec_model no longer writes data/ecModel.python.yml; regenerate the fixtures "
                "with regenerate_fixtures.py"
            )
        # A copy in a scratch folder, so that a loader that wants to write next to its file
        # cannot touch the vendored fixtures.
        python_file = Path(workdir) / "from_python.yml"
        matlab_file = Path(workdir) / "from_matlab.yml"
        shutil.copyfile(data / "ecModel.python.yml", python_file)
        shutil.copyfile(data / "ecModel.matlab.yml", matlab_file)
        return {
            "file_written_by_python": _load(python_file),
            "file_written_by_matlab": _load(matlab_file),
        }
