#!/usr/bin/env bash
# Build the results explorer from benchmark results and serve it locally.
#
# Usage:
#   scripts/explorer.sh [--engine vina|smina|gnina|rdock] [--results DIR] [--port 8000]
#                       [--skip-build] [--open] [--share-url URL]
#
# Reads <results>/<engine>/<arm>/ (default results: results/server/results if present,
# else results), writes explorer/data.js and explorer/poses/<engine>/, then serves
# explorer/ at http://localhost:<port>/ until you press Ctrl+C. Needs uv (scripts/setup.sh
# installs it); the PoseBusters set is downloaded on first use. Works on Linux, macOS and
# Git Bash on Windows.
set -euo pipefail

REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
ENGINE=vina
RESULTS=
PORT=8000
BUILD=1
OPEN=0
SHARE_URL=
while [ $# -gt 0 ]; do
    case $1 in
        --engine) ENGINE=$2; shift 2 ;;
        --results) RESULTS=$2; shift 2 ;;
        --port) PORT=$2; shift 2 ;;
        --skip-build) BUILD=0; shift ;;
        --open) OPEN=1; shift ;;
        --share-url) SHARE_URL=$2; shift 2 ;;
        -h|--help) sed -n '2,12p' "$0"; exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
done

fail() { echo "ERROR: $*" >&2; exit 1; }
command -v uv >/dev/null || fail "uv not found; install it (https://docs.astral.sh/uv/) or run scripts/setup.sh"
cd "$REPO"

if [ "$BUILD" = 1 ]; then
    echo "==> Building explorer data for $ENGINE"
    args=(--engine "$ENGINE")
    [ -n "$RESULTS" ] && args+=(--results "$RESULTS")
    [ -n "$SHARE_URL" ] && args+=(--share-url "$SHARE_URL")
    uv run --frozen python explorer/build_data.py "${args[@]}"
elif [ ! -f explorer/data.js ]; then
    fail "explorer/data.js does not exist; run without --skip-build first"
fi

URL="http://localhost:$PORT/"
echo "==> Serving the explorer at $URL (Ctrl+C to stop)"
if [ "$OPEN" = 1 ]; then
    # Give the server a moment to start, then open the default browser.
    ( sleep 1
      if command -v xdg-open >/dev/null; then xdg-open "$URL"
      elif command -v open >/dev/null; then open "$URL"
      elif command -v cmd.exe >/dev/null; then cmd.exe /c start "" "$URL"
      fi >/dev/null 2>&1 ) &
fi
exec uv run --frozen python -m http.server "$PORT" --bind 127.0.0.1 --directory explorer
