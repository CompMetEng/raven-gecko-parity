"""Python side of the defaultCompartmentMap scenario."""

from raven_toolbox.localization import DEFAULT_COMPARTMENT_MAP


def run(ctx):
    pairs = sorted(f"{label}|{code}" for label, code in DEFAULT_COMPARTMENT_MAP.items())
    return {"count": len(DEFAULT_COMPARTMENT_MAP), "pairs": pairs}
