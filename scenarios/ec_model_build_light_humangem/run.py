"""Python side of the light Human-GEM ecModel build scenario.

The pipeline of GECKO's light_ecModel tutorial, from make_ec_model to set_prot_pool_size,
with a checkpoint after each stage. See scenario.yml.
"""

import math
from pathlib import Path

import cobra
import numpy as np
import yaml

from geckopy import (
    ModelAdapter,
    apply_complex_data,
    apply_kcat_constraints,
    apply_kcat_list,
    assign_standard_kcat,
    fill_eccodes_from_gem,
    fuzzy_kcat_matching,
    load_brenda_data,
    load_conventional_gem,
    load_phyl_dist,
    load_uniprot_tsv,
    make_ec_model,
    merge_dlkcat_and_fuzzy_kcats,
    read_dlkcat_output,
    set_prot_pool_size,
)

from geckopy.databases.kegg_loader import load_kegg_tsv


def num(x) -> dict:
    """A number as a record, so that NaN is a value both sides can state."""
    x = float(x)
    return {"is_nan": bool(math.isnan(x)), "value": 0.0 if math.isnan(x) else x}


def _adapter(inputs: dict) -> ModelAdapter:
    folder = Path(inputs["fixture_dir"])
    adapter = ModelAdapter.from_folder(Path(inputs["adapter_python"]))
    adapter.params.path = folder
    adapter.params.conv_gem = folder / "models" / "human-GEM.yml"
    return adapter


def _patch_eccodes_from_yaml(conv: cobra.Model, yaml_path: Path) -> None:
    """Recover the `eccodes` reaction field of the YAML file, which cobra's reader drops."""
    with open(yaml_path, encoding="utf-8") as f:
        raw = dict(yaml.safe_load(f))
    eccodes_by_id: dict[str, str] = {}
    for row in raw["reactions"]:
        row = dict(row)
        ec = row.get("eccodes")
        if ec:
            eccodes_by_id[row["id"]] = ";".join(ec) if isinstance(ec, list) else str(ec)
    for rxn in conv.reactions:
        if rxn.id in eccodes_by_id:
            rxn.annotation["ec-code"] = eccodes_by_id[rxn.id]


def _source_token(source: str) -> str:
    """The base provenance token: 'brenda (wc=0, origin=1)' -> 'brenda'."""
    return source.split(" (", 1)[0].lower()


def _coupling(ec) -> list[dict]:
    matrix = ec.rxn_enz_mat.tocoo()
    entries = [
        {"reaction": str(ec.rxns[i]), "enzyme": str(ec.enzymes[j]), "coefficient": float(v)}
        for i, j, v in zip(matrix.row, matrix.col, matrix.data)
        if v != 0
    ]
    return sorted(entries, key=lambda e: (e["reaction"], e["enzyme"]))


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


def run(ctx):
    inputs = ctx["inputs"]
    cobra.Configuration().solver = str(inputs["python_solver"])
    adapter = _adapter(inputs)
    conv = load_conventional_gem(adapter)
    conv_counts = {
        "reactions": len(conv.reactions),
        "metabolites": len(conv.metabolites),
        "genes": len(conv.genes),
    }
    _patch_eccodes_from_yaml(conv, adapter.params.conv_gem)

    # MATLAB's makeEcModel reads data/uniprot.tsv and data/kegg.tsv, and uses KEGG for the genes
    # that UniProt misses. make_ec_model reads the data/kegg.tsv of the adapter folder by itself;
    # this fixture keeps its data in another folder, so the databases are passed.
    data = Path(inputs["fixture_dir"]) / "data"
    uniprot_db = load_uniprot_tsv(
        data / "uniprot.tsv",
        id_type=str(adapter.params.uniprot.type),
        gene_id_field=str(adapter.params.uniprot.gene_id_field),
        reviewed=bool(adapter.params.uniprot.reviewed),
        auto_download=False,
    )
    kegg_db = load_kegg_tsv(data / "kegg.tsv")
    model = make_ec_model(conv, adapter, gecko_light=True, uniprot_db=uniprot_db, kegg_db=kegg_db)
    model.solver = str(inputs["python_solver"])
    model.adapter = adapter
    expansion = {
        "conv_counts": conv_counts,
        "ec_counts": {
            "ec_rxns": len(model.ec.rxns),
            "enzymes": len(model.ec.enzymes),
            "genes": len(model.ec.genes),
        },
        "gecko_light": bool(model.ec.gecko_light),
        "ec_rxns": [str(r) for r in model.ec.rxns],
        "enzymes": sorted(str(e) for e in model.ec.enzymes),
        "n_model_reactions": len(model.reactions),
        "n_model_metabolites": len(model.metabolites),
    }

    apply_complex_data(model, apply=False)
    complexes = _coupling(model.ec)

    fill_eccodes_from_gem(model)
    eccodes = [str(c) if c else "" for c in model.ec.eccodes]

    brenda = load_brenda_data(Path(inputs["databases_dir"]) / "brenda-2026.1-1")
    phyl_dist = load_phyl_dist(Path(inputs["adapter_matlab"]).parent / "PhylDist.mat")
    fuzzy_df = fuzzy_kcat_matching(model, brenda, phyl_dist)
    dlkcat_df = read_dlkcat_output(model, Path(inputs["fixture_dir"]) / "data" / "DLKcat.tsv")
    dlkcat_rows = sorted(
        (
            {
                "reaction": str(r.rxn_id),
                "gene": str(r.genes[0]) if r.genes else "",
                "kcat": num(r.kcat),
            }
            for r in dlkcat_df.itertuples()
        ),
        key=lambda r: (r["reaction"], r["gene"], r["kcat"]["value"]),
    )
    merged = merge_dlkcat_and_fuzzy_kcats(dlkcat_df, fuzzy_df)
    merged_rows = sorted(
        (
            {"reaction": str(r.rxn_id), "kcat": num(r.kcat), "source": _source_token(str(r.source))}
            for r in merged.itertuples()
        ),
        key=lambda r: (r["reaction"], r["source"], r["kcat"]["value"]),
    )

    apply_kcat_list(model, merged, criteria="max")
    after_assign = _kcat_rows(model)

    # The standard-kcat step runs on a copy and its result is a checkpoint, whether it raises
    # or not; the later stages continue from the model without it, so that they compare the
    # same state on both sides.
    def standard_step():
        probe = model.copy()
        before = set(str(r) for r in probe.ec.rxns)
        before_enzymes = set(str(e) for e in probe.ec.enzymes)
        assign_standard_kcat(probe, uniprot_db)
        return {
            "rows": _kcat_rows(probe),
            "new_rxns": sorted(set(str(r) for r in probe.ec.rxns) - before),
            "new_enzymes": sorted(set(str(e) for e in probe.ec.enzymes) - before_enzymes),
        }

    standard = _attempt(standard_step)

    apply_kcat_constraints(model)
    pool = set_prot_pool_size(model)
    pool_rxn = model.reactions.get_by_id("prot_pool_exchange")

    return {
        "expansion": expansion,
        "complexes": complexes,
        "eccodes": {
            "populated": sum(1 for c in eccodes if c),
            "per_ec_rxn": eccodes,
        },
        "dlkcat": dlkcat_rows,
        "kcats": {
            "merged": merged_rows,
            "after_assign": after_assign,
            "standard": standard,
            "prot_coefficients": _prot_coefficients(model),
            "nonzero_count": int(np.sum(model.ec.kcat > 0)),
            "pool_ub": num(pool_rxn.upper_bound),
        },
    }
