"""Python side of the complex-data edge-case scenario.

Three groups of checkpoints, matching the three ledger rows this scenario covers:

``complexes``
    apply_complex_data on each file of ``inputs.complex_cases``, for both ecModel flavours:
    the non-zero entries of ec.rxn_enz_mat as sorted (reaction, enzyme, coefficient) records.
``f_factor``
    calculate_f_factor on hand-made proteomics tables. NaN is not emitted (it would read as a
    harness difference): a flag and a zero instead.
``fill_concs``
    fill_enz_concs on hand-made tables, as ec.concs paired with the enzyme, sorted by enzyme.

A case whose file or table is missing fails the run; a case that raised is a flag.
"""

from pathlib import Path

import numpy as np

from geckopy import (
    ModelAdapter,
    apply_complex_data,
    calculate_f_factor,
    fill_enz_concs,
    load_conventional_gem,
    make_ec_model,
)
from geckopy.databases.pax_db_loader import ProtData


def _adapter(inputs: dict) -> ModelAdapter:
    """geckopy's ecTestGEM adapter, repointed at the GECKO copy of the fixture."""
    adapter = ModelAdapter.from_folder(inputs["adapter_python"])
    fixture = Path(inputs["fixture_dir"])
    adapter.params.path = fixture
    adapter.params.conv_gem = fixture / "models" / "testModel.xml"
    return adapter


def _coupling(ec) -> list[dict]:
    """Non-zero entries of ec.rxn_enz_mat as sorted records."""
    matrix = ec.rxn_enz_mat
    if hasattr(matrix, "toarray"):
        matrix = matrix.toarray()
    entries = [
        {"reaction": ec.rxns[i], "enzyme": ec.enzymes[j], "coefficient": float(matrix[i][j])}
        for i in range(len(ec.rxns))
        for j in range(len(ec.enzymes))
        if matrix[i][j] != 0
    ]
    return sorted(entries, key=lambda e: (e["reaction"], e["enzyme"]))


def _concs(model) -> list[dict]:
    """model.ec.concs paired with the enzyme, sorted by enzyme; NaN as a flag and a zero."""
    records = []
    for enzyme, conc in zip(model.ec.enzymes, model.ec.concs):
        value = float(conc)
        has_conc = value == value
        records.append({"enzyme": str(enzyme), "has_conc": has_conc, "conc": value if has_conc else 0.0})
    return sorted(records, key=lambda r: r["enzyme"])


def _prot_data(case: dict) -> ProtData:
    # None (a null in the declaration) is a missing value; np.array(dtype=float) reads it as NaN.
    return ProtData(
        uniprot_ids=[str(i) for i in case["ids"]],
        abundances=np.array(case["abundances"], dtype=float),
    )


def run(ctx):
    inputs = ctx["inputs"]
    adapter = _adapter(inputs)
    model = load_conventional_gem(adapter)

    # --- complexes
    complex_dir = Path(inputs["complex_dir"])
    complexes = {}
    for name in sorted(inputs["complex_cases"]):
        path = Path(inputs["fixture_dir"]) / "data" / "ComplexPortal.json" if name == "fixture" else complex_dir / f"{name}.json"
        # A missing file must fail the run, not match the other side's identical failure.
        if not path.is_file():
            raise FileNotFoundError(path)
        per_flavour = {}
        for flavour, light in (("full", False), ("light", True)):
            ec = make_ec_model(model, adapter, gecko_light=light)
            try:
                apply_complex_data(ec, path=path, apply=False)
                per_flavour[flavour] = {"raised": False, "coupling": _coupling(ec.ec)}
            except Exception:  # noqa: BLE001 --- whether one was raised is the finding
                per_flavour[flavour] = {"raised": True, "coupling": []}
        complexes[name] = per_flavour

    # --- f factor
    ec_for_f = make_ec_model(model, adapter, gecko_light=False)
    f_factor = {}
    for name in sorted(inputs["f_factor_cases"]):
        case = inputs["f_factor_cases"][name]
        kwargs = {"enzymes": []} if case["empty_enzymes"] else {}
        try:
            value = float(calculate_f_factor(ec_for_f, _prot_data(case), **kwargs))
            is_nan = value != value
            f_factor[name] = {"raised": False, "is_nan": is_nan, "value": 0.0 if is_nan else value}
        except Exception:  # noqa: BLE001
            f_factor[name] = {"raised": True, "is_nan": False, "value": 0.0}

    # --- fill concs
    fill_concs = {}
    for name in sorted(inputs["fill_cases"]):
        case = inputs["fill_cases"][name]
        ec = make_ec_model(model, adapter, gecko_light=False)
        try:
            fill_enz_concs(ec, _prot_data(case), data_col=int(case["data_col"]) - 1)
            fill_concs[name] = {"raised": False, "concs": _concs(ec)}
        except Exception:  # noqa: BLE001
            fill_concs[name] = {"raised": True, "concs": []}

    return {"complexes": complexes, "f_factor": f_factor, "fill_concs": fill_concs}
