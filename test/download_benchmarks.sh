#!/usr/bin/env bash

set -e

ZENODO_DOI="10.5281/zenodo.4546263" # Project DOI (redirects to latest version)

# Extract ID from DOI
ZENODO_ID=$(echo $ZENODO_DOI | cut -d"." -f3)
echo $ZENODO_ID

# Follow redirect to latest version record and fetch file list
FILE_LIST=$(curl -sL https://zenodo.org/api/records/$ZENODO_ID | jq -r '.files | .[] | .links.self')

for f in $FILE_LIST; do
  curl -L -OJ $f
done
