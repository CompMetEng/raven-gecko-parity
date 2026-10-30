"""Python side of the Yarrowia ecModel build scenario.

The ecModel pipeline on iYali5-GEM with the data of the eciYali5-GEM repository, with a
checkpoint after each stage. See scenario.yml.
"""

import math
import tempfile
from pathlib import Path

import cobra
import numpy as np
import pandas as pd

from geckopy import (
    ModelAdapter,
    apply_custom_kcats,
    apply_kcat_constraints,
    apply_kcat_list,
    assign_standard_kcat,
    fill_eccodes_from_database,
    fill_eccodes_from_gem,
    fuzzy_kcat_matching,
    load_brenda_data,
    load_conventional_gem,
    load_phyl_dist,
    load_uniprot_tsv,
    make_ec_model,
    set_prot_pool_size,
)
from geckopy.databases.kegg_loader import load_kegg_tsv


def num(x) -> dict:
    """A number as a record, so that NaN is a value both sides can state."""
    x = float(x)
    return {"is_nan": bool(math.isnan(x)), "value": 0.0 if math.isnan(x) else x}


def _adapter(inputs: dict) -> ModelAdapter:
    adapter = ModelAdapter.from_folder(Path(inputs["adapter_python"]))
    root = Path(inputs["eciyali_dir"])
    adapter.params.path = root
    adapter.params.conv_gem = root / "model" / "iYali5-GEM" / "iYali5.xml"
    return adapter


def _source_token(source: str) -> str:
    """The base provenance token: 'brenda (wc=0, origin=1)' -> 'brenda'."""
    return source.split(" (", 1)[0].lower()


def _kcat_rows(model) -> list[dict]:
    return sorted(
        (
            {
                "reaction": str(r),
                "kcat": num(model.ec.kcat[i]),
                "source": _source_token(str(model.ec.source[i])) if model.ec.source[i] else "",
            }
            for i, r in enumerate(model.ec.rxns)
        ),
        key=lambda row: row["reaction"],
    )


def _attempt(fn) -> dict:
    """A raise is a result on both sides."""
    try:
        return {"raised": False, "value": fn()}
    except Exception as error:
        print(f"case raised: {type(error).__name__}: {error}")
        return {"raised": True, "value": {}}


def _prot_coefficients(model) -> list[dict]:
    records = []
    for rxn in model.reactions:
        for met, coeff in rxn.metabolites.items():
            if met.id.startswith("prot_") and coeff != 0:
                records.append({"reaction": rxn.id, "metabolite": met.id, "coefficient": float(coeff)})
    return sorted(records, key=lambda e: (e["reaction"], e["metabolite"]))


def _eccodes(model) -> list[str]:
    return [str(c) if c else "" for c in model.ec.eccodes]


def _pipeline(inputs: dict, data_root: Path) -> dict:
    """The pipeline with data/ read from data_root."""
    adapter = _adapter(inputs)
    adapter.params.path = data_root
    conv = load_conventional_gem(adapter)
    conv_counts = {
        "reactions": len(conv.reactions),
        "metabolites": len(conv.metabolites),
        "genes": len(conv.genes),
    }

    # MATLAB's makeEcModel reads data/uniprot.tsv and data/kegg.tsv, and uses KEGG for the genes
    # that UniProt misses; make_ec_model uses KEGG only when it is given the database.
    data = data_root / "data"
    uniprot = load_uniprot_tsv(
        data / "uniprot.tsv",
        id_type=str(adapter.params.uniprot.type),
        gene_id_field=str(adapter.params.uniprot.gene_id_field),
        reviewed=bool(adapter.params.uniprot.reviewed),
        auto_download=False,
    )
    kegg = load_kegg_tsv(data / "kegg.tsv")
    model = make_ec_model(conv, adapter, gecko_light=False, uniprot_db=uniprot, kegg_db=kegg)
    model.solver = str(inputs["python_solver"])
    model.adapter = adapter
    expansion = {
        "conv_counts": conv_counts,
        "ec_counts": {
            "ec_rxns": len(model.ec.rxns),
            "enzymes": len(model.ec.enzymes),
            "genes": len(model.ec.genes),
        },
        "ec_rxns": [str(r) for r in model.ec.rxns],
        "enzymes": [str(e) for e in model.ec.enzymes],
        "mw": [num(m) for m in model.ec.mw],
        "n_model_reactions": len(model.reactions),
        "n_model_metabolites": len(model.metabolites),
    }

    # The model has no EC codes of its own. getECfromGEM raises for that; whether the call
    # raises is a result, and the codes after it are recorded either way.
    try:
        fill_eccodes_from_gem(model)
        from_gem_raised = False
    except Exception:
        from_gem_raised = True
    from_gem = _eccodes(model)

    fill_eccodes_from_database(model, uniprot, kegg_db=kegg)
    from_database = _eccodes(model)

    brenda = load_brenda_data(Path(inputs["databases_dir"]) / "brenda-2026.1-1")
    phyl_dist = load_phyl_dist(Path(inputs["adapter_matlab"]).parent / "PhylDist.mat")
    fuzzy_df = fuzzy_kcat_matching(model, brenda, phyl_dist)
    fuzzy_rows = sorted(
        (
            {
                "reaction": str(r.rxn_id),
                "kcat": num(r.kcat),
                "eccode": str(r.eccode) if r.eccode else "",
                "origin": -1 if pd.isna(r.origin) else int(r.origin),
                "wildcard_level": -1 if pd.isna(r.wildcard_level) else int(r.wildcard_level),
            }
            for r in fuzzy_df.itertuples()
        ),
        key=lambda r: r["reaction"],
    )

    apply_kcat_list(model, fuzzy_df, criteria="max")
    after_assign = _kcat_rows(model)

    apply_custom_kcats(model, data / "customKcats.tsv", apply=False)
    after_custom = _kcat_rows(model)

    # The standard-kcat step runs on a copy and its result is a checkpoint, whether it raises
    # or not; the later stages continue from the model without it, so that they compare the
    # same state on both sides.
    def standard_step():
        probe = model.copy()
        before = set(str(r) for r in probe.ec.rxns)
        before_enzymes = set(str(e) for e in probe.ec.enzymes)
        assign_standard_kcat(probe, uniprot)
        return {
            "rows": _kcat_rows(probe),
            "new_rxns": sorted(set(str(r) for r in probe.ec.rxns) - before),
            "new_enzymes": sorted(set(str(e) for e in probe.ec.enzymes) - before_enzymes),
        }

    standard = _attempt(standard_step)

    apply_kcat_constraints(model)
    set_prot_pool_size(model)
    pool_rxn = model.reactions.get_by_id("prot_pool_exchange")

    return {
        "expansion": expansion,
        "eccodes": {
            "from_gem_raised": from_gem_raised,
            "from_gem": from_gem,
            "from_database": from_database,
            "populated_after_database": sum(1 for c in from_database if c),
        },
        "kcats": {
            "fuzzy": fuzzy_rows,
            "after_assign": after_assign,
            "after_custom": after_custom,
            "standard": standard,
            "prot_coefficients": _prot_coefficients(model),
            "nonzero_count": int(np.sum(model.ec.kcat > 0)),
            "pool_ub": num(pool_rxn.upper_bound),
        },
    }


def run(ctx):
    inputs = ctx["inputs"]
    cobra.Configuration().solver = str(inputs["python_solver"])
    root = Path(inputs["eciyali_dir"])
    results = {"as_is": _pipeline(inputs, root)}

    # The UniProt and KEGG snapshots quote the cells that hold several values, such as
    # "1.1.1.25; 2.5.1.19" in the EC column. `unquoted` runs the pipeline on a copy of the data
    # folder with every double quote removed, which separates what depends on the quotes from
    # what does not.
    with tempfile.TemporaryDirectory(prefix="parity_yarrowia_") as workdir:
        copy = Path(workdir)
        (copy / "data").mkdir()
        for name in ("uniprot.tsv", "kegg.tsv", "customKcats.tsv"):
            text = (root / "data" / name).read_text(encoding="utf-8")
            (copy / "data" / name).write_text(text.replace('"', ""), encoding="utf-8", newline="")
        results["unquoted"] = _pipeline(inputs, copy)
    return results
