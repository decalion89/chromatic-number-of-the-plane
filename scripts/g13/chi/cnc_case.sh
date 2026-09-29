#!/bin/bash
# Cube and conquer for one case formula, with a checked cover and a drat-trim-checked proof per leaf.
#
# usage: cnc_case.sh CASE.cnf DEPTH VARFILE JOBS [TIME]
#   CASE.cnf  a case formula (plan_C.py, g13cnf.py --rigid34 ..., enum_cert.py)
#   DEPTH     number of branching decisions for cuber2.py; VARFILE lists the branching variables, comma separated:
#             work/vars_c0_lex.txt = x(v,0) along lex_order() (colouring formulas),
#             work/vars_s_lex.txt  = s_v along lex_order() (independent-set formulas E(t))
#   JOBS      parallel certification jobs (1 while the machine is shared, 4 when it is free)
#   TIME      kissat limit per leaf in seconds (default 1200)
#   CUBER_ARGS (environment) extra options for cuber2.py, e.g. "--minpos 11" (stop at 11 points of C0)
# Output: CASE.icnf (leaves; 'c closed' = refuted by unit propagation, certified all the same), CASE.certlog, and
# CASE.todo = leaves not VERIFIED (timeouts): split them again with cuber2.py on CASE.cnf + that cube, or rerun
# certify.py --only with a larger TIME.  Leaf formulas and proofs live only while being checked (disk).
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
cnf=$1; depth=$2; vars=$(cat "$3"); jobs=$4; t=${5:-1200}
stem=${cnf%.cnf}
[ -f "$stem.icnf" ] || nice -n 19 python3 $HERE/cuber2.py "$cnf" "$depth" "$stem.icnf" --vars "$vars" $CUBER_ARGS
nice -n 19 python3 $HERE/certify.py --log "$stem.certlog" --jobs "$jobs" --time "$t" --base "$cnf" --cubes "$stem.icnf"
grep -v "drat-trim VERIFIED" "$stem.certlog" > "$stem.todo" || true
echo "$(grep -c 'drat-trim VERIFIED' "$stem.certlog") leaves verified, $(wc -l < "$stem.todo") not (see $stem.todo)"
