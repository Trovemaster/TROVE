#!/usr/bin/env bash

ZENODO_DOI="10.5281/zenodo.4433456" # Project DOI (redirects to latest version)

# Extract ID from DOI
ZENODO_ID=$(echo $ZENODO_DOI | cut -d"." -f3)

# Get correct URL from redirect
REDIRECT_SUFFIX=$(curl -s https://zenodo.org/api/records/$ZENODO_ID | jq -r '.location')
REDIRECT_URL="https://zenodo.org"$REDIRECT_SUFFIX
echo $REDIRECT_URL

# Fetch file list
FILE_LIST=$(curl -s $REDIRECT_URL | jq -r '.files | .[] | .links.self')

for f in $FILE_LIST; do
  curl -OJ $f
done
