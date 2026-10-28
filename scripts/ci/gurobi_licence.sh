#!/usr/bin/env bash
# Writes a Gurobi licence file from the environment and exports GRB_LICENSE_FILE to
# later workflow steps. Scenarios that need a MILP solver read it; the others ignore it.
#
# Two shapes are accepted because two exist in the organisation: LICENCE_FILE holds a
# licence file verbatim (the GUROBI_EDUARD organisation secret), and WLSACCESSID,
# WLSSECRET and LICENSEID are the three-field form that raven-toolbox's parity-nightly
# workflow uses.
#
# A missing licence produces a notice, not an error: most scenarios are solver-free
# and still run without one.
set -euo pipefail

umask 077
licence="$HOME/gurobi.lic"
if [ -n "${LICENCE_FILE:-}" ]; then
  printf '%s\n' "$LICENCE_FILE" > "$licence"
elif [ -n "${WLSACCESSID:-}" ] && [ -n "${WLSSECRET:-}" ] && [ -n "${LICENSEID:-}" ]; then
  { echo "WLSACCESSID=$WLSACCESSID"
    echo "WLSSECRET=$WLSSECRET"
    echo "LICENSEID=$LICENSEID"
  } > "$licence"
else
  echo "::notice::no Gurobi licence configured; solver-dependent scenarios will fail"
  exit 0
fi

# Report the shape without printing the value. A named-user academic licence is tied
# to one machine and cannot authenticate on a hosted runner, so finding one here
# explains a failure that otherwise looks like a solver bug.
if grep -q '^[[:space:]]*WLSACCESSID=' "$licence"; then
  echo "Gurobi WLS licence installed"
else
  echo "::warning::the configured licence is not a WLS licence. A named-user" \
       "licence is tied to a specific machine and will not authenticate on a" \
       "hosted runner; solver-dependent scenarios will fail."
fi
echo "GRB_LICENSE_FILE=$licence" >> "$GITHUB_ENV"
