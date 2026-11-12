#!/usr/bin/env bash
#
# Fallback download of the full TROVE benchmark archive from Zenodo.
#
# The H2CO inputs and the reference files compared by scripts/H2CO/compare_results.sh
# are tracked in test/benchmarks/H2CO, so run_regression_tests.sh no longer needs this
# script. Zenodo may block automated clients (e.g. HTTP 403 "Access to this resource has
# been restricted due to unusual traffic from your network" on GitHub-hosted runners), so
# fail with a clear message instead of passing an HTML error page to jq.

set -euo pipefail

ZENODO_DOI="10.5281/zenodo.4546263" # Project DOI (redirects to latest version)

# Extract ID from DOI
ZENODO_ID=$(echo $ZENODO_DOI | cut -d"." -f3)
echo $ZENODO_ID

# Follow redirect to latest version record and fetch file list
if ! RECORD=$(curl -sSL --fail --retry 3 --retry-delay 5 "https://zenodo.org/api/records/$ZENODO_ID"); then
  echo "ERROR: Zenodo API request for record $ZENODO_ID failed (Zenodo may be blocking this network)." >&2
  exit 3
fi
FILE_LIST=$(echo "$RECORD" | jq -r '.files | .[] | .links.self')

for f in $FILE_LIST; do
  curl -sSL --fail --retry 3 --retry-delay 5 -OJ "$f"
done
