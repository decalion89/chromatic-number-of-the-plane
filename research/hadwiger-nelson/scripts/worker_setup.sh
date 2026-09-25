#!/bin/sh
# Set up a fresh machine to run the growth searches: python packages, kissat and drat-trim, tabu2.
# Idempotent. Tools go to $TOOLS (default $HOME/hn-tools); prints the paths to export.
set -e
TOOLS=${TOOLS:-$HOME/hn-tools}
mkdir -p "$TOOLS"
python3 -m pip install -q numpy scipy sympy python-sat 2>/dev/null || pip install -q numpy scipy sympy python-sat
# The revisions used for the results of release 1.0.0 (research/hadwiger-nelson/requirements-lock.txt).
KISSAT_TAG=rel-4.0.4
DRAT_TRIM_COMMIT=2e3b2dc0ecf938addbd779d42877b6ed69d9a985
if [ ! -x "$TOOLS/kissat/build/kissat" ]; then
  (cd "$TOOLS" && rm -rf kissat && git clone -q --depth 1 --branch "$KISSAT_TAG" https://github.com/arminbiere/kissat.git && cd kissat && ./configure > /dev/null && make -j4 > /dev/null)
fi
if [ ! -x "$TOOLS/drat-trim/drat-trim" ]; then
  (cd "$TOOLS" && rm -rf drat-trim && git clone -q https://github.com/marijnheule/drat-trim.git && cd drat-trim && git checkout -q "$DRAT_TRIM_COMMIT" && make > /dev/null)
fi
HERE=$(cd "$(dirname "$0")" && pwd)
gcc -O2 -o "$HERE/tabu2" "$HERE/tabu2.c"
echo "KISSAT=$TOOLS/kissat/build/kissat"
echo "DRAT_TRIM=$TOOLS/drat-trim/drat-trim"
