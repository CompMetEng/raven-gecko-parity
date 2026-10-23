"""Python side of the kcat review statistics scenario.

Every case builds its own state, runs the call and records what it returned and whether it
raised. A raise is a value (`raised`). Nothing reaches the network.
"""

import shutil
from pathlib import Path

import numpy as np

from geckopy.databases import load_brenda_data, load_brenda_measurements
from geckopy.gather_kcats.review import EcStats, ec_stats_from_brenda, review_assignment

# One set of assignments for every review case. Columns: reaction, kcat, source, EC code,
# leverage. The comment names the check each row is meant to meet.
ROWS = [
    ("R_A", 10.0, "brenda", "1.1.1.1", 9.0),       # equals the EC maximum, five times its median
    ("R_B", 2.0, "custom", "1.1.1.1", 8.0),        # custom value equal to a median
    ("R_C", 6.0, "Custom", "1.1.1.1", 7.0),        # 'Custom' with a capital, a value equal to a maximum
    ("R_D", 5.0, "brenda", "4.1.1.1", 6.0),        # EC with specific-activity rows only
    ("R_F", 40.0, "custom", "3.5.1.n3", 5.5),       # lowercase letter in the EC, value equal to its maximum
    ("R_G1", 77.0, "brenda", "", 4.0),             # the same value on five reactions
    ("R_G2", 77.0, "brenda", "", 3.9),
    ("R_G3", 77.0, "brenda", "", 3.8),
    ("R_G4", 77.0, "brenda", "", 3.7),
    ("R_G5", 77.0, "brenda", "", 3.6),
    ("R_L", 0.001, "brenda", "", 3.0),             # slow
    ("R_M", 50000.0, "STANDARD", "", 2.0),         # 'STANDARD' with capitals, fast
    ("R_N", 20000.0, "standard", "", 1.5),         # standard fallback, fast
    ("R_P_EXP_1", 0.005, "brenda", "", 1.2),       # two isozyme copies sharing one value, slow
    ("R_P_EXP_2", 0.005, "brenda", "", 1.1),
]


def _num(value) -> dict:
    value = float(value)
    if np.isnan(value):
        return {"nan": True, "value": 0.0}
    return {"nan": False, "value": value}


def _stats_records(stats) -> list:
    rows = []
    for ec in sorted(stats):
        s = stats[ec]
        rows.append({"ec": ec, "n_kcat": int(s.n_kcat), "kcat_max": _num(s.kcat_max),
                     "kcat_median": _num(s.kcat_median), "values": sorted(float(v) for v in s.values),
                     "n_sa": int(s.n_sa)})
    return rows


def _hand_stats() -> dict:
    return {
        "1.1.1.1": EcStats(n_kcat=3, kcat_max=10.0, kcat_median=2.0, values=frozenset({2.0, 6.0, 10.0}), n_sa=0),
        "4.1.1.1": EcStats(n_kcat=0, n_sa=1),
        "3.5.1.n3": EcStats(n_kcat=1, kcat_max=40.0, kcat_median=4.0, values=frozenset({4.0, 40.0}), n_sa=0),
    }


def _review(stats, **kwargs) -> dict:
    rxns = [r[0] for r in ROWS]
    kcats = np.array([r[1] for r in ROWS])
    sources = [r[2] for r in ROWS]
    ecs = [r[3] for r in ROWS]
    lev = np.array([r[4] for r in ROWS])
    try:
        found = review_assignment(rxns, kcats, sources, ecs, ec_stats=stats, leverage=lev, **kwargs)
    except Exception:  # noqa: BLE001 --- whether the call raises is the finding, not the type
        return {"raised": True, "findings": []}
    return {"raised": False,
            "findings": [{"rxn": f.rxn_id, "ec": f.ec_code, "source": f.source, "kcat": float(f.kcat),
                          "leverage": float(f.leverage), "n_isozymes": int(f.n_isozymes),
                          "other_reactions": int(f.other_reactions), "checks": list(f.checks)}
                         for f in found]}


def _loader_stats(inputs) -> dict:
    return ec_stats_from_brenda(load_brenda_data(Path(inputs["brenda_dir"])))


def _case_lowercase_mw(inputs) -> dict:
    """A BRENDA folder whose mw.tsv and sa.tsv hold an EC with a lowercase letter."""
    try:
        stats = ec_stats_from_brenda(load_brenda_data(Path(inputs["brenda_lower_mw_dir"])))
    except Exception:  # noqa: BLE001
        return {"raised": True, "stats": []}
    return {"raised": False, "stats": _stats_records(stats)}


def _case_measurements(inputs) -> dict:
    out = {}
    for name in ("meas_tsv", "meas_xz", "meas_both"):
        try:
            table = load_brenda_measurements(Path(inputs[name]))
        except Exception:  # noqa: BLE001
            out[name] = {"raised": True, "n": 0, "release": "", "rows": []}
            continue
        rows = sorted(({"ec": str(r.ec), "kcat": float(r.kcat), "log10": float(r.log10_kcat)}
                       for r in table.itertuples()), key=lambda r: (r["ec"], r["kcat"]))
        out[name] = {"raised": False, "n": len(table), "release": str(table.attrs.get("release", "")),
                     "rows": rows}
    return out


def run(ctx):
    inputs = ctx["inputs"]
    loader = _loader_stats(inputs)
    hand = _hand_stats()
    optional = {"include": ("standard-fallback", "tied-isozymes")}
    return {
        "ec_stats": _stats_records(loader),
        "review_hand_stats": _review(hand),
        "review_hand_stats_optional": _review(hand, **optional),
        "review_hand_stats_top": _review(hand, top=3),
        "review_loader_stats": _review(loader),
        "lowercase_ec_in_mw": _case_lowercase_mw(inputs),
        "measurements": _case_measurements(inputs),
    }
