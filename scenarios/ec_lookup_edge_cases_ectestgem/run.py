"""Python side of the EC-lookup and subsetting edge-case scenario.

Three checkpoints (see scenario.yml). A call that raises is reported as `raised: true` with an
empty result, because which exception it was is not the finding.
"""

from pathlib import Path

import cobra

from geckopy import (
    ModelAdapter,
    copy_ec_to_gem,
    fill_eccodes_from_database,
    fill_eccodes_from_gem,
    get_subset_ec_model,
    load_conventional_gem,
    load_uniprot_tsv,
    make_ec_model,
)
from geckopy.databases.kegg_loader import load_kegg_tsv


def _bound(value: float) -> tuple[float, str]:
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


def _eccodes(model) -> list[dict]:
    """ec.eccodes as one record per ec reaction, sorted by reaction id."""
    ec = model.ec
    records = [{"reaction": str(r), "eccodes": str(c) if c else ""} for r, c in zip(ec.rxns, ec.eccodes)]
    return sorted(records, key=lambda e: e["reaction"])


def _database_check(inputs: dict, check: dict) -> dict:
    adapter = _adapter(inputs)
    model = make_ec_model(load_conventional_gem(adapter), adapter, gecko_light=False)
    root = Path(inputs["variants_dir"]) / str(check["variant"])
    uniprot = load_uniprot_tsv(root / "data" / "uniprot.tsv")
    kegg = load_kegg_tsv(root / "data" / "kegg.tsv")
    subset = [str(r) for r in check["subset"]] or None
    try:
        fill_eccodes_from_database(model, uniprot, kegg_db=kegg, ec_rxns=subset, action=str(check["action"]))
    except Exception:  # noqa: BLE001 --- whether it raised is the finding, not the exception type
        return {"raised": True, "eccodes": []}
    return {"raised": False, "eccodes": _eccodes(model)}


def _eccode_tokens(model, rxn_id: str) -> list[str]:
    value = model.reactions.get_by_id(rxn_id).annotation.get("ec-code")
    if not value:
        return []
    tokens = value.split(";") if isinstance(value, str) else list(value)
    return sorted(t for t in tokens if t)


def _apply_eccodes(inputs: dict, light: bool) -> dict:
    adapter = _adapter(inputs)
    model = make_ec_model(load_conventional_gem(adapter), adapter, gecko_light=light)
    fill_eccodes_from_gem(model)
    try:
        copy_ec_to_gem(model, overwrite=False)
    except Exception:  # noqa: BLE001
        return {"raised": True, "reactions": []}
    records = [{"reaction": r.id, "eccodes": _eccode_tokens(model, r.id)}
               for r in sorted(model.reactions, key=lambda r: r.id)]
    return {"raised": False, "reactions": records}


def _small_gem(variant: str, conv: "cobra.Model") -> "cobra.Model":
    small = conv.copy()
    if variant == "drop_reaction":
        small.remove_reactions(["R3"], remove_orphans=True)
    elif variant == "drop_gene":
        small.reactions.get_by_id("R2").gene_reaction_rule = "G1 and G2"
        cobra.manipulation.delete.remove_genes(small, ["G3"], remove_reactions=False)
    else:
        raise ValueError(f"unknown subset variant {variant}")
    return small


def _summary(model) -> dict:
    reactions = []
    for rxn in sorted(model.reactions, key=lambda r: r.id):
        lower, lower_kind = _bound(rxn.lower_bound)
        upper, upper_kind = _bound(rxn.upper_bound)
        reactions.append({
            "id": rxn.id, "lower_bound": lower, "lower_kind": lower_kind,
            "upper_bound": upper, "upper_kind": upper_kind,
        })
    associations = sorted(
        ({"reaction": rxn.id, "gene": gene.id} for rxn in model.reactions for gene in rxn.genes),
        key=lambda e: (e["reaction"], e["gene"]),
    )
    ec = model.ec
    matrix = ec.rxn_enz_mat
    if hasattr(matrix, "toarray"):
        matrix = matrix.toarray()
    coupling = sorted(
        ({"reaction": ec.rxns[i], "enzyme": ec.enzymes[j], "coefficient": float(matrix[i][j])}
         for i in range(len(ec.rxns)) for j in range(len(ec.enzymes)) if matrix[i][j] != 0),
        key=lambda e: (e["reaction"], e["enzyme"]),
    )
    return {
        "reactions": reactions,
        "metabolites": sorted(met.id for met in model.metabolites),
        "genes": sorted(gene.id for gene in model.genes),
        "gene_associations": associations,
        "ec_rxns": sorted(str(r) for r in ec.rxns),
        "ec_enzymes": sorted(str(e) for e in ec.enzymes),
        "ec_genes": sorted(str(g) for g in ec.genes),
        "ec_coupling": coupling,
    }


def _subset(inputs: dict, variant: str) -> dict:
    adapter = _adapter(inputs)
    conv = load_conventional_gem(adapter)
    big = make_ec_model(conv, adapter, gecko_light=False)
    small = _small_gem(variant, conv)
    try:
        sub = get_subset_ec_model(big, small)
    except Exception:  # noqa: BLE001
        return {"raised": True, "model": {}}
    return {"raised": False, "model": _summary(sub)}


def run(ctx):
    inputs = ctx["inputs"]
    return {
        "database_checks": {str(c["name"]): _database_check(inputs, c) for c in inputs["database_checks"]},
        "apply_eccodes": {
            "full": _apply_eccodes(inputs, False),
            "light": _apply_eccodes(inputs, True),
        },
        "subset": {str(v): _subset(inputs, str(v)) for v in inputs["subset_variants"]},
    }
