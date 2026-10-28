"""Python side of the GPR normalisation scenario.

Both sides are handed the rule string. ``grRuleToDNF`` and ``isDnfGrRule`` read
it with ``parseGrRule``; ``gpr_to_dnf`` and ``is_dnf`` read it with
``parse_gr_rule``, the same grammar. A malformed rule raises
``RAVEN:badGrRule`` there and :class:`ValueError` here, and is recorded as
unparsed.
"""

from raven_toolbox.manipulation import gpr_to_dnf
from raven_toolbox.utils import is_dnf


def run(ctx):
    rules = list(ctx["inputs"]["rules"])

    # Input order, not sorted: the rules are a fixture whose order is part of
    # the declaration, and keeping it makes a difference easy to locate in
    # scenario.yml.
    analysed = [_analyse(str(rule)) for rule in rules]

    return {
        "n_rules": len(analysed),
        # A rule that fails to parse is the interesting case, so count it
        # rather than leaving it to be spotted in the per-rule list.
        "n_unparsable": sum(1 for row in analysed if not row["parsed"]),
        "n_dnf": sum(1 for row in analysed if row["is_dnf"]),
        "rules": analysed,
    }


def _analyse(rule):
    try:
        clauses = [[str(gene) for gene in clause] for clause in gpr_to_dnf(rule)]
        dnf = bool(is_dnf(rule))
    except ValueError:
        # Every key always present: an unparsable rule on one side only must
        # read as a value difference, not as a structural one.
        return {"rule": rule, "parsed": False, "is_dnf": False, "n_clauses": 0, "clauses": []}

    return {
        "rule": rule,
        "parsed": True,
        "is_dnf": dnf,
        "n_clauses": len(clauses),
        "clauses": clauses,
    }
