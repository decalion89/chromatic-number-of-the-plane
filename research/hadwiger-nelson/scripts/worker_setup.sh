#!/bin/sh
# Set up a fresh machine to run the growth searches: python packages, kissat and drat-trim, tabu2.
# Idempotent. Tools go to $TOOLS (default $HOME/hn-tools); prints the paths to export.
set -e
TOOLS=${TOOLS:-$HOME/hn-tools}
mkdir -p "$TOOLS"
python3 -m pip install -q numpy scipy sympy python-sat 2>/dev/null || pip install -q numpy scipy sympy python-sat
if [ ! -x "$TOOLS/kissat/build/kissat" ]; then
  (cd "$TOOLS" && rm -rf kissat && git clone -q --depth 1 https://github.com/arminbiere/kissat.git && cd kissat && ./configure > /dev/null && make -j4 > /dev/null)
fi
if [ ! -x "$TOOLS/drat-trim/drat-trim" ]; then
  (cd "$TOOLS" && rm -rf drat-trim && git clone -q --depth 1 https://github.com/marijnheule/drat-trim.git && cd drat-trim && make > /dev/null)
fi
HERE=$(cd "$(dirname "$0")" && pwd)
gcc -O2 -o "$HERE/tabu2" "$HERE/tabu2.c"
echo "KISSAT=$TOOLS/kissat/build/kissat"
echo "DRAT_TRIM=$TOOLS/drat-trim/drat-trim"
