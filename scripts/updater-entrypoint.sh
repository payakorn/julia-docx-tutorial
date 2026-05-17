#!/bin/sh
# Entrypoint for the updater sidecar.
# Installs git + docker CLI, registers the cron job, then runs crond in the
# foreground so the container stays alive.
set -eu

apk add --no-cache --quiet git openssh-client docker-cli tzdata >/dev/null

CRON_SCHEDULE="${CRON_SCHEDULE:-0 6 * * *}"
LOG_FILE="${LOG_FILE:-/var/log/updater.log}"

mkdir -p "$(dirname "$LOG_FILE")"
touch "$LOG_FILE"

echo "${CRON_SCHEDULE} /repo/scripts/repo-update.sh >> ${LOG_FILE} 2>&1" \
    > /etc/crontabs/root

echo "[updater] cron schedule: ${CRON_SCHEDULE} (TZ=${TZ:-UTC})"
echo "[updater] logging to ${LOG_FILE}"

# Stream the log so `docker logs updater` shows what happened.
tail -F "$LOG_FILE" &

exec crond -f -l 8
