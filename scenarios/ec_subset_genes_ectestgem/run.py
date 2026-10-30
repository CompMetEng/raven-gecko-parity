"""Python side of the ecModel subset scenario.

get_subset_ec_model on ecTestGEM with context models that lack a reaction and a gene of its
rule. A raise is a result. See scenario.yml.
"""

import importlib.util
import sys
from pathlib import Path

import cobra

from geckopy import ModelAdapter, get_subset_ec_model, load_conventional_gem, make_ec_model


def _adapter(inputs: dict) -> ModelAdapter:
    base = ModelAdapter.from_folder(inputs["adapter_python"])
    fixture = Path(inputs["fixture_dir"])
    base.params.path = fixture
    base.params.conv_gem = fixture / "models" / "testModel.xml"
    adapter_module_path = Path(inputs["adapter_python"]) / "adapter.py"
    spec = importlib.util.spec_from_file_location("ectestgem_adapter", adapter_module_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module.TestGEMAdapter(base.params)


def _small_gem(variant: str, conv: cobra.Model) -> cobra.Model:
    small = conv.copy()
    if variant == "drop_reaction_and_gene":
        small.remove_reactions(["R2"], remove_orphans=False)
        small.genes.remove(small.genes.get_by_id("G3"))
    elif variant == "drop_reaction":
        small.remove_reactions(["R2"], remove_orphans=True)
    elif variant != "identity":
        raise ValueError(variant)
    return small


def _summary(model) -> dict:
    ec = model.ec
    return {
        "reactions": sorted(r.id for r in model.reactions),
        "genes": sorted(g.id for g in model.genes),
        "gene_associations": sorted(
            f"{r.id}:{g.id}" for r in model.reactions for g in r.genes
        ),
        "ec_rxns": sorted(str(r) for r in ec.rxns),
        "ec_genes": sorted(str(g) for g in ec.genes),
        "ec_enzymes": sorted(str(e) for e in ec.enzymes),
    }


def _subset(inputs: dict, variant: str, light: bool) -> dict:
    adapter = _adapter(inputs)
    conv = load_conventional_gem(adapter)
    big = make_ec_model(conv, adapter, gecko_light=light)
    small = _small_gem(variant, conv)
    try:
        sub = get_subset_ec_model(big, small)
    except Exception as error:
        print(f"case raised: {type(error).__name__}: {error}")
        return {"raised": True, "model": {}}
    return {"raised": False, "model": _summary(sub)}


def run(ctx):
    inputs = ctx["inputs"]
    return {
        flavour: {v: _subset(inputs, v, light) for v in ("identity", "drop_reaction", "drop_reaction_and_gene")}
        for flavour, light in (("full", False), ("light", True))
    }
