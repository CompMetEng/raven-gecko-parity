"""Python side of the OpenKineticsPredictor (OKP) scenario.

Three checkpoints on one ecModel built from ecTestGEM, all through kcats_from_okp: find_met_smiles
feeding the two write checkpoints (the job input that kcats_from_okp builds, with the submission
stopped before it reaches the network), then the stored result of a job parsed into a kcat list
(use_stored=True). See scenario.yml for why this stops at the CSV boundary rather than the real
OKP REST API, and for the one real MATLAB bug this scenario's write checkpoint found and
confirmed fixed ([GECKO #437]).

The read_output checkpoint compares the source column verbatim: both sides carry the
'OKP-' prefix.
"""

import copy
import importlib.util
import os
import shutil
import sys
import tempfile
from pathlib import Path
from unittest import mock

from geckopy import ModelAdapter, kcats_from_okp, load_conventional_gem, make_ec_model
from geckopy.databases.pubchem import find_met_smiles
from geckopy.gather_kcats.open_kinetics_predictor import OKPError
from geckopy.gather_kcats.open_kinetics_predictor import run as okp_run

METHOD = "CataPro"
JOB_ID = "scenario-job"


def _adapter(inputs: dict) -> ModelAdapter:
    """geckopy's ecTestGEM adapter --- see kcat_chain_ectestgem/run.py for why the
    actual TestGEMAdapter subclass is needed rather than the generic ModelAdapter.
    """
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


def _smiles_snapshot(model) -> dict:
    return {m.id: str(m.annotation.get("smiles", "")) for m in sorted(model.metabolites, key=lambda m: m.id)}


def _csv_pairs(csv_text: str) -> list[str]:
    lines = [line for line in csv_text.strip("\n").split("\n")[1:] if line]
    return sorted(lines)


class _OfflineClient:
    """Stands in for OKPClient: kcats_from_okp has written the job input when it asks the
    client to submit, and the submission stops there, so no key is used and nothing is sent."""

    def __init__(self, *args, **kwargs):
        pass

    def submit(self, *args, **kwargs):
        raise OKPError("offline")


def _project_adapter(adapter, project: Path):
    """A copy of the adapter whose project folder is a scratch folder."""
    copied = copy.deepcopy(adapter)
    copied.params.path = project
    (project / "data").mkdir(parents=True, exist_ok=True)
    return copied


def _job_input(model, adapter, ec_rxns: list[str]) -> list[str]:
    """The pairs of the job input that kcats_from_okp writes for the selected reactions."""
    with tempfile.TemporaryDirectory() as tmp:
        project = Path(tmp)
        scratch = _project_adapter(adapter, project)
        # The API key is resolved before the input is written; a placeholder is enough because
        # the offline client never sends it.
        with mock.patch.dict(os.environ, {"OKP_API_KEY": "offline"}), \
                mock.patch.object(okp_run, "OKPClient", _OfflineClient):
            try:
                kcats_from_okp(model, method=METHOD, ec_rxns=ec_rxns, save=False, adapter=scratch)
            except OKPError:
                pass
        return _csv_pairs((project / "data" / "OKP.csv").read_text(encoding="utf-8"))


def _checkpoint_smiles(model, inputs: dict) -> dict:
    find_met_smiles(model, cache_path=Path(inputs["scenario_root"]) / "data" / "smilesDB.tsv")
    return _smiles_snapshot(model)


def _checkpoint_write(model, adapter) -> dict:
    # ec.rxns order on this fixture is [R2_EXP_1, R2_EXP_2, R2_REV_EXP_1, R2_REV_EXP_2,
    # R3, R5] (confirmed by direct execution, identical on both sides). Requesting
    # {R2_EXP_1, R2_EXP_2, R3} skips the two REV entries sitting in between --- exactly
    # the subset shape that GECKO #437's bug mishandled, so this checkpoint doubles as
    # its regression check.
    with_smiles = _job_input(model, adapter, ["R2_EXP_1", "R2_EXP_2", "R3"])

    # m2 (R5's substrate) deliberately loses its SMILES here: kcats_from_okp always drops
    # an entry without a SMILES from the job input, so R5's entry does not appear at all.
    model.metabolites.get_by_id("m2c").annotation["smiles"] = ""
    with_filter = _job_input(model, adapter, ["R2_EXP_1", "R2_EXP_2", "R3", "R5"])

    return {
        "with_smiles": with_smiles,
        "after_clearing_m2_with_filter": with_filter,
    }


def _checkpoint_read(model, adapter, inputs: dict) -> dict:
    with tempfile.TemporaryDirectory() as tmp:
        project = Path(tmp)
        scratch = _project_adapter(adapter, project)
        data = project / "data"
        shutil.copyfile(inputs["okp_result_path"], data / "OKP_output_catapro.csv")
        (data / "OKP_job_catapro.txt").write_text(
            f"jobId: {JOB_ID}\nmethod: {METHOD}\nresultOf: {JOB_ID}\n", encoding="utf-8")
        df, done = kcats_from_okp(model, method=METHOD, use_stored=True, save=False, adapter=scratch)
    rows = sorted(
        (
            {
                "reaction": str(r.rxn_id),
                "gene": str(r.genes[0]) if r.genes else "",
                "substrate": str(r.substrates[0]) if r.substrates else "",
                "kcat": float(r.kcat),
                # Recorded verbatim: 'OKP-CataPro', 'OKP-BRENDA' on both sides.
                "source": str(r.source),
            }
            for r in df.itertuples()
        ),
        key=lambda r: (r["reaction"], r["gene"]),
    )
    return {"done": bool(done), "rows": rows}


def run(ctx):
    inputs = ctx["inputs"]
    adapter = _adapter(inputs)
    conv = load_conventional_gem(adapter)
    model = make_ec_model(conv, adapter, gecko_light=False)

    smiles_result = _checkpoint_smiles(model, inputs)
    write_result = _checkpoint_write(model, adapter)
    read_result = _checkpoint_read(model, adapter, inputs)

    return {
        "smiles": smiles_result,
        "write": write_result,
        "read_output": read_result,
    }
