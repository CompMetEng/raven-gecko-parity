"""Python side of the molecular-weight scenario.

``calculate_mw`` takes a bare string and returns a float --- no model, no adapter, nothing
to resolve. Every sequence in the declaration is fed through it once and the results
returned in declaration order, not sorted: docs/scenarios.md's default is to sort
everything, but this fixture is a flat list of independent, hand-chosen test cases rather
than model output, and keeping input order is what makes a result that stops matching
scenario.yml's per-sequence comments easy to locate --- the same reasoning
``gpr_dnf_rules`` uses for its own rule list.

``calculate_mw`` returns NaN for a sequence with no recognisable residues, and so does
MATLAB's ``calculateMW``. The NaN is emitted as-is, not folded into a has-value/value pair
the way an absent metabolite charge is elsewhere in this repo: the point of this scenario is
that both sides give the same answer for the same input.
"""

from geckopy import calculate_mw


def run(ctx):
    sequences = list(ctx["inputs"]["sequences"])
    return {
        "sequences": [{"sequence": sequence, "mw": float(calculate_mw(sequence))} for sequence in sequences],
    }
