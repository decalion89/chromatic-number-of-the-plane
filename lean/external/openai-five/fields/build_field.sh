#!/bin/sh
# Rebuild OpenAI's proof that the plane is not 5-colourable over a smaller field, and print the axioms.
#
# usage: sh build_field.sh FIELD MATH_LEAN_DIR LEAN_BIN LEAN_PATH WORK_DIR [JOBS]
#
# FIELD          constructible, origami or radicals
# MATH_LEAN_DIR  the lean/ directory of a checkout of https://github.com/openai/math at commit
#                fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb (it holds OAI/...); it is not modified
# LEAN_BIN       the lean binary of toolchain v4.34.1
# LEAN_PATH      the Mathlib build, as for ../build_closure.py
# WORK_DIR       a new directory: the patched sources go to WORK_DIR/src, the .olean files to WORK_DIR/build
# JOBS           parallel compilations (default 2)
#
# The last lines printed are the output of `#print axioms` for the statements of the field (for the constructible
# field also the finite form, `ConstructibleFinite.lean`).
set -eu
here=$(cd "$(dirname "$0")" && pwd)
field=$1; math=$2; lean=$3; lp=$4; work=$5; jobs=${6:-2}
case $field in
  constructible) check=ConstructiblePlane ;;
  origami) check=OrigamiPlane ;;
  radicals) check=RadPlane ;;
  *) echo "unknown field: $field (use constructible, origami or radicals)" >&2; exit 2 ;;
esac
if [ -e "$work" ]; then echo "$work exists; give a new directory" >&2; exit 2; fi
mkdir -p "$work/src/Check"
cp -R "$math/OAI" "$work/src/"
(cd "$work/src" && patch -p1 --no-backup-if-mismatch < "$here/$field.patch")
cp "$here/$check.lean" "$work/src/Check/"
python3 "$here/../build_closure.py" "$work/src" "$lean" "$lp" "$work/build" "$jobs" "Check.$check"
cat "$work/build/Check/$check.out"
if [ "$field" = constructible ]; then
  # the finite form, through Mathlib's compactness theorem; it imports the statement file just built
  cp "$here/ConstructibleFinite.lean" "$work/src/Check/"
  (cd "$work/src" && LEAN_PATH="$work/build:$lp" "$lean" -o "$work/build/Check/ConstructibleFinite.olean" \
    Check/ConstructibleFinite.lean) > "$work/build/Check/ConstructibleFinite.out" 2>&1 || true
  cat "$work/build/Check/ConstructibleFinite.out"
fi
