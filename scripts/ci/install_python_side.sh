#!/usr/bin/env bash
# Installs the Python side of a pair from the checkouts under repos/ into the current
# Python, and writes the parity.local.toml that points parity at those checkouts.
#
# Usage: scripts/ci/install_python_side.sh raven|gecko
#
# The gecko pair needs four repos, not two: GECKO depends on RAVEN on the MATLAB path,
# and geckopy depends on raven-toolbox on the Python side. raven-toolbox is therefore
# installed for both pairs.
set -euo pipefail

pair="$1"

cat > parity.local.toml <<'EOF'
[pairs.raven.matlab]
path = "repos/RAVEN"
[pairs.raven.python]
path = "repos/raven-toolbox"
EOF
pip install "./repos/raven-toolbox[excel]"

# gurobipy registers "gurobi" as a solver interface with cobra/optlang; the licence
# file from gurobi_licence.sh has no effect unless the package is installed. A Python
# version without a gurobipy wheel still runs the solver-free scenarios, and the
# report marks the solver-dependent ones as ERROR.
pip install gurobipy || echo "::warning::gurobipy is not installable on $(python -V); solver-dependent scenarios will fail"

if [ "$pair" = "gecko" ]; then
  cat >> parity.local.toml <<'EOF'
[pairs.gecko.matlab]
path = "repos/GECKO"
[pairs.gecko.python]
path = "repos/geckopy"
EOF
  # --no-deps: geckopy's own metadata pins raven-toolbox to a git ref, which would
  # otherwise replace the checkout installed above with whatever that ref currently
  # resolves to. Its other dependencies are read from its own pyproject.toml
  # (skipping the raven-toolbox entry), so this list follows geckopy's changes.
  pip install --no-deps ./repos/geckopy
  python -c "
import tomllib
with open('repos/geckopy/pyproject.toml', 'rb') as f:
    deps = tomllib.load(f)['project']['dependencies']
print('\n'.join(d for d in deps if not d.strip().startswith('raven-toolbox')))
" > geckopy-deps.txt
  pip install -r geckopy-deps.txt
fi
