#!/usr/bin/env bash
# Set up DockPrepBench on a Linux x86_64 machine: Python environment and docking engines.
#
# Usage:
#   scripts/setup.sh [--tools-dir DIR] [--no-gnina] [--skip-tests]
#
# Everything is installed in user space under the tools directory (default: <repo>/tools,
# or $DOCKPREPBENCH_TOOLS), including uv's package cache, so nothing lands in your home
# directory. Safe to rerun: installed components are skipped. Afterwards:
#   source <tools-dir>/env.sh
set -euo pipefail

# Pinned versions -------------------------------------------------------------------------
GNINA_URL=https://github.com/gnina/gnina/releases/download/v1.3.3/gnina.cuda12.8.static
GNINA_BYTES=2056131000  # exact size of the release asset, checked after download
SMINA_URL=https://sourceforge.net/projects/smina/files/smina.static/download
RDOCK_URL=https://github.com/CBDD/rDock/releases/download/v24.03.192/rdock-24.03.192_ubuntu-22.04-clang_clang++_x86_64.tar.gz
MICROMAMBA_URL=https://micro.mamba.pm/api/micromamba/linux-64/latest

# Options ---------------------------------------------------------------------------------
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
TOOLS=${DOCKPREPBENCH_TOOLS:-$REPO/tools}
WITH_GNINA=1
RUN_TESTS=1
while [ $# -gt 0 ]; do
    case $1 in
        --tools-dir) TOOLS=$2; shift 2 ;;
        --no-gnina) WITH_GNINA=0; shift ;;
        --skip-tests) RUN_TESTS=0; shift ;;
        -h|--help) sed -n '2,11p' "$0"; exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
done

step() { printf '\n==> %s\n' "$*"; }
fail() { echo "ERROR: $*" >&2; exit 1; }

[ "$(uname -s)" = Linux ] && [ "$(uname -m)" = x86_64 ] \
    || fail "Linux x86_64 is required (the vina bindings and engine binaries are Linux-only)"
mkdir -p "$TOOLS/bin"
TOOLS=$(cd "$TOOLS" && pwd)
export PATH=$TOOLS/bin:$PATH
export UV_CACHE_DIR=${UV_CACHE_DIR:-$TOOLS/uv-cache}
echo "Repository: $REPO"
echo "Tools:      $TOOLS"

# Python environment ----------------------------------------------------------------------
step "uv"
if ! command -v uv >/dev/null; then
    curl -LsSf https://astral.sh/uv/install.sh \
        | env UV_INSTALL_DIR="$TOOLS/bin" UV_NO_MODIFY_PATH=1 sh >/dev/null
fi
uv --version

step "Python environment (uv sync --locked)"
(cd "$REPO" && nice -n 19 uv sync --locked)

# Engines ---------------------------------------------------------------------------------
step "Smina"
if [ ! -x "$TOOLS/smina/smina" ]; then
    mkdir -p "$TOOLS/smina"
    curl -fsSL -o "$TOOLS/smina/smina.part" "$SMINA_URL"
    file "$TOOLS/smina/smina.part" | grep -q ELF || fail "Smina download is not a Linux binary"
    mv "$TOOLS/smina/smina.part" "$TOOLS/smina/smina" && chmod +x "$TOOLS/smina/smina"
fi
ln -sf "$TOOLS/smina/smina" "$TOOLS/bin/smina"

if [ "$WITH_GNINA" = 1 ]; then
    step "Gnina (2 GB; resumes if interrupted)"
    if [ ! -x "$TOOLS/gnina/gnina" ]; then
        mkdir -p "$TOOLS/gnina"
        part=$TOOLS/gnina/gnina.part
        # GitHub's HTTP/2 transfers of this file can drop mid-way; HTTP/1.1 with resume works.
        until [ "$(stat -c %s "$part" 2>/dev/null)" = "$GNINA_BYTES" ]; do
            nice -n 19 curl --http1.1 -fL -C - -o "$part" "$GNINA_URL" || sleep 15
        done
        mv "$part" "$TOOLS/gnina/gnina" && chmod +x "$TOOLS/gnina/gnina"
    fi
    ln -sf "$TOOLS/gnina/gnina" "$TOOLS/bin/gnina"
fi

step "rDock 24.03 (the 2013.1 conda build segfaults on current systems)"
if [ ! -x "$TOOLS/rdock/bin/rbdock" ]; then
    mkdir -p "$TOOLS/rdock"
    curl -fsSL "$RDOCK_URL" | tar -xz -C "$TOOLS/rdock"
fi
for t in rbdock rbcavity; do ln -sf "$TOOLS/rdock/bin/$t" "$TOOLS/bin/$t"; done

step "Open Babel (receptor conversion for rDock)"
export MAMBA_ROOT_PREFIX=$TOOLS/micromamba
if [ ! -x "$TOOLS/bin/micromamba" ]; then
    curl -fsSL "$MICROMAMBA_URL" | tar -xj -C "$TOOLS" bin/micromamba
fi
if [ ! -x "$MAMBA_ROOT_PREFIX/envs/openbabel/bin/obabel" ]; then
    nice -n 19 micromamba create -y -q -n openbabel --override-channels -c conda-forge openbabel >/dev/null
fi
ln -sf "$MAMBA_ROOT_PREFIX/envs/openbabel/bin/obabel" "$TOOLS/bin/obabel"

# Environment file ------------------------------------------------------------------------
step "Writing $TOOLS/env.sh"
cat > "$TOOLS/env.sh" <<EOF
# Source this to use DockPrepBench's tools: source $TOOLS/env.sh
export PATH="$TOOLS/bin:\$PATH"
export UV_CACHE_DIR="$UV_CACHE_DIR"
export RBT_ROOT="$TOOLS/rdock"
export LD_LIBRARY_PATH="$TOOLS/rdock/lib\${LD_LIBRARY_PATH:+:\$LD_LIBRARY_PATH}"
EOF
# shellcheck source=/dev/null
source "$TOOLS/env.sh"

# Check -----------------------------------------------------------------------------------
step "Versions"
smina --version 2>&1 | head -1
[ "$WITH_GNINA" = 1 ] && { gnina --version 2>&1 | head -1; }
rbdock -h 2>&1 | grep -m1 -i "rbdock" || true
obabel -V
if command -v nvidia-smi >/dev/null; then
    echo "GPUs: $(nvidia-smi --query-gpu=name --format=csv,noheader | sort | uniq -c | tr -s ' ' | paste -sd ';')"
fi

if [ "$RUN_TESTS" = 1 ]; then
    step "Tests (including real docking with every installed engine)"
    (cd "$REPO" && nice -n 19 uv run pytest -q -rs)
fi

step "Done. Use the tools with:  source $TOOLS/env.sh"
