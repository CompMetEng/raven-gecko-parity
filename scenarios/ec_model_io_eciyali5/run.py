"""Python side of the GECKO 3 ecModel scenario on eciYali5-GEM.

load_ec_model reads a real ecModel written by GECKO 3 (reverse protein reactions), the
upgrade is run on the file read without it, and the model is written with save_ec_model and
read back. Every loaded model is summarised field by field. See scenario.yml.
"""

import gzip
import re
import tempfile
import warnings
from pathlib import Path

from geckopy import ModelAdapter, load_ec_model, save_ec_model
from geckopy.upgrade import upgrade_ec_model
from raven_toolbox.io import read_yaml_model


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


def _protein_reactions(model) -> list[dict]:
    """The usage reactions and the pool exchange: direction, bounds and objective weight."""
    out = []
    for rxn in model.reactions:
        if rxn.id.startswith("usage_prot_") or rxn.id == "prot_pool_exchange":
            lower, lower_kind = _bound(rxn.lower_bound)
            upper, upper_kind = _bound(rxn.upper_bound)
            out.append(
                {
                    "id": rxn.id,
                    "lower_bound": lower,
                    "lower_kind": lower_kind,
                    "upper_bound": upper,
                    "upper_kind": upper_kind,
                    "objective_coefficient": float(rxn.objective_coefficient),
                    "coefficients": sorted(
                        ({"metabolite": m.id, "coefficient": float(c)} for m, c in rxn.metabolites.items()),
                        key=lambda e: e["metabolite"],
                    ),
                }
            )
    return sorted(out, key=lambda e: e["id"])


def _attempt(fn) -> dict:
    """A raise is a result on both sides."""
    try:
        return {"raised": False, "value": fn()}
    except Exception as error:
        print(f"case raised: {type(error).__name__}: {error}")
        return {"raised": True, "value": {}}


def _quiet(fn):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return fn()


def run(ctx):
    inputs = ctx["inputs"]
    adapter = ModelAdapter.from_folder(Path(inputs["adapter_python"]))
    results = {}
    with tempfile.TemporaryDirectory(prefix="parity_eciyali_") as workdir:
        work = Path(workdir)
        published = work / "published.yml"
        with gzip.open(inputs["model_gz"], "rb") as source:
            text = source.read().decode("utf-8")
        published.write_text(text, encoding="utf-8", newline="")

        # The file has two constructs that a strict YAML parser rejects. GECKO 3 wrote three
        # kcat notes as ""text"" and 30 reactions with the standard kcat as an `enzymes: !!omap`
        # entry with no content. The cases as_published and notes_repaired record what each
        # reader does with the file as published and with the notes repaired (both readers load
        # both); the other cases read a copy in which both are repaired (doubled quotes become
        # single ones, the empty `enzymes` entries are left out).
        def counts_of(path):
            return _summary(_quiet(lambda: load_ec_model(path, adapter)))["counts"]

        results["as_published"] = _attempt(lambda: counts_of(published))
        notes_fixed = re.sub(r'^( *- notes: )""([^\n]+)""$', r'\1"\2"', text, flags=re.MULTILINE)
        notes_path = work / "notes_repaired.yml"
        notes_path.write_text(notes_fixed, encoding="utf-8", newline="")
        results["notes_repaired"] = _attempt(lambda: counts_of(notes_path))
        repaired = re.sub(r"^ *- enzymes: !!omap\n(?= {0,4}\S)", "", notes_fixed, flags=re.MULTILINE)
        original = work / "eciYali5-GEM.yml"
        original.write_text(repaired, encoding="utf-8", newline="")

        loaded = _quiet(lambda: load_ec_model(original, adapter))
        results["loaded"] = _summary(loaded)
        results["loaded_protein_reactions"] = _protein_reactions(loaded)

        def raw_then_upgrade():
            raw = _quiet(lambda: read_yaml_model(str(original)))
            changed = upgrade_ec_model(raw)
            return {"changed": len(changed), "protein_reactions": _protein_reactions(raw), "after": _summary(raw)}

        results["upgrade_raw"] = _attempt(raw_then_upgrade)

        def round_trip():
            written = work / "written.yml"
            save_ec_model(loaded, written)
            again = _quiet(lambda: load_ec_model(written, adapter))
            return {"reloaded": _summary(again), "protein_reactions": _protein_reactions(again)}

        results["round_trip"] = _attempt(round_trip)

        def second_round_trip():
            first = work / "first.yml"
            second = work / "second.yml"
            save_ec_model(loaded, first)
            once = _quiet(lambda: load_ec_model(first, adapter))
            save_ec_model(once, second)
            twice = _quiet(lambda: load_ec_model(second, adapter))
            return {"stable": _summary(once) == _summary(twice)}

        results["second_round_trip"] = _attempt(second_round_trip)
    return results
