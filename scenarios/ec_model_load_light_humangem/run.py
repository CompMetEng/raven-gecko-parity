"""Python side of the light Human-GEM ecModel loading scenario.

``load_ec_model`` reads the light ecModel of GECKO's light_ecModel tutorial as published and
as repaired, and the repaired model is summarised field by field. See scenario.yml.
"""

import re
import tempfile
import warnings
from pathlib import Path

from geckopy import ModelAdapter, load_ec_model


def _bound(value: float) -> tuple[float, str]:
    """An infinite bound as a class and a zero, since JSON has no infinity."""
    value = float(value)
    if value == float("inf"):
        return 0.0, "+inf"
    if value == float("-inf"):
        return 0.0, "-inf"
    return value, "finite"


def _number(value) -> tuple[float, str]:
    """A number that may be NaN as a class and a zero, since JSON has no NaN."""
    value = float(value)
    if value != value:
        return 0.0, "nan"
    return _bound(value)


def _names(value) -> list[str]:
    if value is None or value == "":
        return []
    if isinstance(value, str):
        return [value]
    return sorted(str(v) for v in value if str(v) != "")


def _ec_codes(rxn) -> str:
    value = rxn.annotation.get("ec-code")
    if value is None:
        return ""
    return ";".join([value] if isinstance(value, str) else [str(v) for v in value])


def _reactions(model) -> list[dict]:
    records = []
    for rxn in model.reactions:
        lower, lower_kind = _bound(rxn.lower_bound)
        upper, upper_kind = _bound(rxn.upper_bound)
        records.append(
            {
                "id": rxn.id,
                "name": str(rxn.name or ""),
                "lower_bound": lower,
                "lower_kind": lower_kind,
                "upper_bound": upper,
                "upper_kind": upper_kind,
                "objective_coefficient": float(rxn.objective_coefficient),
                "subsystem": _names(rxn.subsystem),
                "ec_codes": _ec_codes(rxn),
                "genes": sorted(g.id for g in rxn.genes),
            }
        )
    return records


def _gene_rules(model) -> list[dict]:
    return [{"id": rxn.id, "rule": str(rxn.gene_reaction_rule or "")} for rxn in model.reactions]


def _metabolites(model) -> list[dict]:
    records = []
    for met in model.metabolites:
        records.append(
            {
                "id": met.id,
                "name": str(met.name or ""),
                "compartment": str(met.compartment or ""),
                "formula": str(met.formula or ""),
                "charge_set": met.charge is not None,
                "charge": float(met.charge or 0.0),
            }
        )
    return records


def _stoichiometry(model) -> list[dict]:
    entries = [
        {"reaction": rxn.id, "metabolite": met.id, "coefficient": float(coefficient)}
        for rxn in model.reactions
        for met, coefficient in rxn.metabolites.items()
    ]
    return sorted(entries, key=lambda e: (e["reaction"], e["metabolite"]))


def _ec_data(ec) -> dict:
    matrix = ec.rxn_enz_mat.tocoo()
    coupling = [
        {"reaction": ec.rxns[i], "enzyme": ec.enzymes[j], "coefficient": float(v)}
        for i, j, v in zip(matrix.row, matrix.col, matrix.data)
        if v != 0
    ]
    mw = [_number(m) for m in ec.mw]
    concs = [_number(c) for c in ec.concs]
    return {
        "rxns": [str(r) for r in ec.rxns],
        "genes": [str(g) for g in ec.genes],
        "enzymes": [str(e) for e in ec.enzymes],
        "mw": [m[0] for m in mw],
        "mw_kind": [m[1] for m in mw],
        "sequence": [str(s) for s in ec.sequence],
        "concs": [c[0] for c in concs],
        "concs_kind": [c[1] for c in concs],
        "eccodes": [str(c) if c else "" for c in ec.eccodes],
        "kcat": [float(k) for k in ec.kcat],
        "source": [str(s) if s else "" for s in ec.source],
        "notes": [str(s) if s else "" for s in ec.notes],
        "coupling": sorted(coupling, key=lambda e: (e["reaction"], e["enzyme"])),
    }


def _summary(model) -> dict:
    return {
        "counts": {
            "reactions": len(model.reactions),
            "metabolites": len(model.metabolites),
            "genes": len(model.genes),
            "ec_rxns": len(model.ec.rxns),
            "enzymes": len(model.ec.enzymes),
        },
        "model_id": str(model.id or ""),
        "gecko_light": bool(model.ec.gecko_light),
        "reaction_order": [r.id for r in model.reactions],
        "metabolite_order": [m.id for m in model.metabolites],
        "gene_order": [g.id for g in model.genes],
        "compartments": sorted(str(c) for c in model.compartments),
        "reactions": _reactions(model),
        "gene_rules": _gene_rules(model),
        "metabolites": _metabolites(model),
        "stoichiometry": _stoichiometry(model),
        "ec": _ec_data(model.ec),
    }


def _attempt(fn) -> dict:
    """A raise is a result on both sides."""
    try:
        return {"raised": False, "value": fn()}
    except Exception as error:
        print(f"case raised: {type(error).__name__}: {error}")
        return {"raised": True, "value": {}}


def _load(path, adapter):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return load_ec_model(path, adapter)


def run(ctx):
    inputs = ctx["inputs"]
    adapter = ModelAdapter.from_folder(Path(inputs["adapter_python"]))
    source = Path(inputs["adapter_python"]) / "models" / str(inputs["model_file"])
    text = source.read_text(encoding="utf-8")

    with tempfile.TemporaryDirectory(prefix="parity_load_light_") as workdir:
        published = Path(workdir) / "published.yml"
        published.write_text(text, encoding="utf-8", newline="")
        # The file has 52 ec reactions written as an `enzymes: !!omap` entry with no content
        # (the reactions that carry the standard kcat). as_published records what each reader
        # does with the file as it is; `repaired` is the summary of a copy without those empty
        # entries, which both readers can load.
        repaired = re.sub(r"^ *- enzymes: !!omap\n(?= {0,4}\S)", "", text, flags=re.MULTILINE)
        fixed = Path(workdir) / "repaired.yml"
        fixed.write_text(repaired, encoding="utf-8", newline="")
        return {
            "as_published": _attempt(lambda: _summary(_load(published, adapter))["counts"]),
            "repaired": _summary(_load(fixed, adapter)),
        }
