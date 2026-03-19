#!/usr/bin/env bash
# Build the taxonium website and deploy to raven.
# Run from the taxonium_website/ directory.
set -euo pipefail

REMOTE_USER="vlad"
REMOTE_HOST="raven"
REMOTE_STAGING="~/taxonium_website_dist/"
REMOTE_WWW="/var/www/taxonium"

echo "==> Building..."
npm run build

echo "==> Syncing to ${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_STAGING}..."
rsync -avzP dist/ "${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_STAGING}"

echo "==> Installing on server..."
ssh "${REMOTE_USER}@${REMOTE_HOST}" "sudo cp -r ${REMOTE_STAGING}* ${REMOTE_WWW}/"

echo "==> Done. Site updated at ${REMOTE_WWW}."
