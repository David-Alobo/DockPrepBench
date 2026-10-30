#!/usr/bin/env bash
# Publish the results explorer to GitHub Pages.
#
# Builds the explorer's data from local results, then force-pushes the page and its data as the
# only commit on the gh-pages branch. The data is about 270 MB (one pose file per complex), so
# no history is kept: each publish replaces the last. Pages must serve the gh-pages branch
# (repository Settings > Pages > Deploy from a branch > gh-pages, / (root)).
#
# Usage (repository root):
#   scripts/publish_explorer.sh                                  # every engine with results
#   scripts/publish_explorer.sh --engine vina --engine gnina     # passed on to build_data.py
#   scripts/publish_explorer.sh --skip-build                     # publish the data already built
set -euo pipefail

build=1
args=()
for arg in "$@"; do
    if [[ "$arg" == "--skip-build" ]]; then build=0; else args+=("$arg"); fi
done

root=$(git rev-parse --show-toplevel)
cd "$root"
remote=$(git remote get-url origin)
source_commit=$(git rev-parse --short HEAD)

if (( build )); then
    uv run python explorer/build_data.py "${args[@]}"
    # The independent Astex check, if its results are present (results/server/astex/).
    if [[ -d results/server/astex ]]; then
        uv run python explorer/build_data.py --dataset astex "${args[@]}"
    fi
fi
[[ -f explorer/data.js && -d explorer/poses ]] || { echo "No explorer data; build it first." >&2; exit 1; }

site=$(mktemp -d)
trap 'rm -rf "$site"' EXIT
cp explorer/index.html explorer/app.js explorer/lib.js explorer/data.js "$site/"
cp -r explorer/poses "$site/"
if [[ -f explorer/data_astex.js && -d explorer/poses_astex ]]; then
    cp explorer/data_astex.js "$site/" && cp -r explorer/poses_astex "$site/"
fi
touch "$site/.nojekyll"  # serve files as they are, without a Jekyll build

git -C "$site" init -q -b gh-pages
git -C "$site" add -A
git -C "$site" -c user.name="$(git config user.name)" -c user.email="$(git config user.email)" \
    commit -q -m "Explorer built from $source_commit"
git -C "$site" push -q --force "$remote" gh-pages
echo "Published explorer (from $source_commit) to the gh-pages branch."
