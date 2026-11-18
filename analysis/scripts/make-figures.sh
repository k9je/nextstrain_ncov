#!/usr/bin/env bash
# Render the SARS-CoV-2 phylogeographic figures.
#
#   ncov_deme_global        demes = country in SE Asia, region elsewhere
#   ncov_country_global     demes = country everywhere
#   ncov_deme_sea           the deme map cropped to South & Southeast Asia
#   ncov_singapore_global   only the transitions touching Singapore
#   ncov_singapore_sea      the same, cropped
#
# Every map highlights transitions into the focus deme in red and out of it in
# blue, and writes a 330 dpi PNG plus a vector PDF.
#
# Usage:  bash scripts/make-figures.sh [focus-deme]
set -uo pipefail

PYTHON=${PYTHON:-/hpc/home/ansellim/micromamba/envs/denvfig/bin/python}
FOCUS=${1:-Singapore}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$HERE"
OUT=${OUT:-"$HERE/figures"}
mkdir -p "$OUT"

JSON=auspice/ncov_deme.json
[[ -f $JSON ]] || { echo "ERROR: $JSON not built yet"; exit 1; }

# South & Southeast Asia crop. Radii are in degrees, so a cropped map needs a
# much smaller pie scale than a global one.
SEA_EXTENT=(72 142 -12 38)
SEA_ARGS=(--pie-scale 0.24 --line-scale 1.15 --label-fontsize 8.5 --font-scale 1.5)

# Decluttering for the global maps: SARS-CoV-2 yields far more inferred
# transitions than dengue, so keep the strongest lines and merge pie slices
# below 3% into a single grey wedge.
# SARS-CoV-2 has 1736 transition branches (max 176 for one pair, median 5),
# so the dengue-era 1-2/3-9/10+ bins put 45 of 70 lines in the top bin and
# every one came out the same thick black. Rescaled to the actual spread.
GLOBAL_ARGS=(--max-lines 55 --min-slice-fraction 0.03 --line-bins 9,49)
# SARS-CoV-2 pies are coloured by Nextstrain clade, not dengue's major_lineage.
LINEAGE=(--lineage-key clade_nextstrain)

fail=0
run() { echo "--- $1"; shift; $PYTHON -I scripts/plot-deme-map.py "$@" 2>&1 \
          | grep -vi downloadwarning | sed 's/^/    /' || ((fail++)); }

run "global, SE Asia countries + world regions" \
  --auspice-json "$JSON" --focus-deme "$FOCUS" "${LINEAGE[@]}" "${GLOBAL_ARGS[@]}" \
  --output "$OUT/ncov_deme_global.png" \
  --title "SARS-CoV-2 global phylogeography: demes = country in SE Asia, region elsewhere"

run "global, country everywhere" \
  --auspice-json "$JSON" --focus-deme "$FOCUS" "${LINEAGE[@]}" --deme-key country "${GLOBAL_ARGS[@]}" \
  --output "$OUT/ncov_country_global.png" \
  --title "SARS-CoV-2 global phylogeography: demes = country everywhere"

run "cropped to South & Southeast Asia" \
  --auspice-json "$JSON" --focus-deme "$FOCUS" "${LINEAGE[@]}" \
  --extent "${SEA_EXTENT[@]}" "${SEA_ARGS[@]}" --min-slice-fraction 0.03 --line-bins 9,49 \
  --output "$OUT/ncov_deme_sea.png" \
  --title "SARS-CoV-2 phylogeography, South & Southeast Asia: demes = country in SE Asia, region elsewhere"

run "global, $FOCUS transitions only" \
  --auspice-json "$JSON" --focus-deme "$FOCUS" "${LINEAGE[@]}" --focus-only --min-slice-fraction 0.03 --line-bins 9,49 \
  --output "$OUT/ncov_singapore_global.png" \
  --title "SARS-CoV-2: inferred $FOCUS imports (red) and exports (blue), global demes"

run "cropped, $FOCUS transitions only" \
  --auspice-json "$JSON" --focus-deme "$FOCUS" "${LINEAGE[@]}" --focus-only \
  --extent "${SEA_EXTENT[@]}" "${SEA_ARGS[@]}" --min-slice-fraction 0.03 --line-bins 9,49 \
  --output "$OUT/ncov_singapore_sea.png" \
  --title "SARS-CoV-2: inferred $FOCUS imports (red) and exports (blue), South & Southeast Asia"

echo
echo "=== wrote $(ls "$OUT"/*.png 2>/dev/null | wc -l) PNGs / $(ls "$OUT"/*.pdf 2>/dev/null | wc -l) PDFs into $OUT; $fail failure(s) ==="
