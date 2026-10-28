"""Python side of the getINITSteps scenario."""

from raven_toolbox.init import get_init_steps


def run(ctx):
    results = {}
    for case in ctx["inputs"]["cases"]:
        try:
            steps = get_init_steps(str(case["series"]))
        except ValueError:
            results[case["name"]] = {"error": True}
            continue
        results[case["name"]] = {
            "steps": [
                {
                    "how_to_use_prev": step.how_to_use_prev,
                    "ignore_mask": [int(b) for b in step.ignore_mask],
                    "pos_rev_off": bool(step.pos_rev_off),
                    "allow_met_secr": bool(step.allow_met_secr),
                    "runs": [
                        {"mip_gap": float(r["mip_gap"]), "time_limit": float(r["time_limit"]), "abs_gap": float(r["abs_gap"])}
                        for r in step.milp_runs
                    ],
                }
                for step in steps
            ]
        }
    return results
