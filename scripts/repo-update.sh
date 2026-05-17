#!/bin/sh
# Daily updater: fast-forward main from origin and restart the server
# container if anything moved. Run from the updater sidecar via cron.
set -eu

REPO_DIR="${REPO_DIR:-/repo}"
BRANCH="${BRANCH:-main}"
SERVER_SERVICE="${SERVER_SERVICE:-server}"

cd "$REPO_DIR"

GIT="git -c safe.directory=${REPO_DIR}"

echo "[$(date -u +%FT%TZ)] checking $BRANCH"
$GIT fetch --quiet origin "$BRANCH"

LOCAL=$($GIT rev-parse HEAD)
REMOTE=$($GIT rev-parse "origin/$BRANCH")

if [ "$LOCAL" = "$REMOTE" ]; then
    echo "  up to date at $LOCAL"
    exit 0
fi

echo "  updating $LOCAL -> $REMOTE"
$GIT pull --ff-only --quiet origin "$BRANCH"

echo "  rebuilding + restarting $SERVER_SERVICE"
docker compose -f "$REPO_DIR/docker-compose.yml" up -d --build "$SERVER_SERVICE"
echo "  done"
