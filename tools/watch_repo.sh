#!/usr/bin/env bash
# watch_repo.sh — Local macOS notifications for remote repo changes
# Usage: ./tools/watch_repo.sh [interval_seconds]
# Requires: terminal-notifier (brew install terminal-notifier)

INTERVAL="${1:-300}"
REPO_DIR="$(cd "$(dirname "$0")/.." && pwd)"

while true; do
    cd "$REPO_DIR" || exit 1
    git fetch origin --quiet 2>/dev/null
    LOCAL=$(git rev-parse HEAD 2>/dev/null)
    REMOTE=$(git rev-parse @{u} 2>/dev/null)
    if [ "$LOCAL" != "$REMOTE" ]; then
        if command -v terminal-notifier &>/dev/null; then
            terminal-notifier -title "Repo Update" \
                -message "Remote has new commits in $(basename "$REPO_DIR")" \
                -sound default
        else
            echo "[$(date)] Remote has new commits — pull to update."
        fi
    fi
    sleep "$INTERVAL"
done
