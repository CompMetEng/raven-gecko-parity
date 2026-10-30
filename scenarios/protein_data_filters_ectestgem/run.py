"""Python side of the protein-data filter scenario.

One checkpoint per case in ``inputs.cases``: ``load_prot_data`` on a vendored replicate table
with one filter switched on at a time. Each case is emitted as ``{raised, rows}`` where rows
are sorted by UniProt id and hold, per condition, a presence flag and a value (NaN is not
emitted: MATLAB's jsonencode writes it as null, which would read as a difference between the
harnesses and not between the loaders).
"""

from pathlib import Path

from geckopy import load_prot_data


def _rows(prot_data) -> list[dict]:
    rows = []
    for uniprot_id, values in zip(prot_data.uniprot_ids, prot_data.abundances):
        values = [float(v) for v in values]
        present = [v == v for v in values]  # NaN != NaN
        rows.append(
            {
                "id": str(uniprot_id),
                "present": present,
                "value": [v if p else 0.0 for v, p in zip(values, present)],
            }
        )
    return sorted(rows, key=lambda r: r["id"])


def run(ctx):
    results = {}
    data_dir = Path(ctx["inputs"]["data_dir"])
    for name in sorted(ctx["inputs"]["cases"]):
        case = ctx["inputs"]["cases"][name]
        path = data_dir / case["file"]
        # A missing table must fail the run: a case that raised because its file is absent
        # would otherwise match the other side's identical failure.
        if not path.is_file():
            raise FileNotFoundError(path)
        repl = [int(r) for r in case["repl"]]
        max_missing = case["max_missing"]
        if isinstance(max_missing, list):
            max_missing = [float(m) for m in max_missing]
        else:
            max_missing = float(max_missing)
        try:
            prot_data = load_prot_data(
                path,
                repl,
                filter_data=bool(case["filter_data"]),
                min_val=float(case["min_val"]),
                max_rsd=float(case["max_rsd"]),
                max_missing=max_missing,
                cut_lowest=float(case["cut_lowest"]),
                add_stdevs=float(case["add_stdevs"]),
            )
            results[name] = {"raised": False, "rows": _rows(prot_data)}
        except Exception:  # noqa: BLE001 --- whether one was raised is the finding, not its type
            results[name] = {"raised": True, "rows": []}
    return results
