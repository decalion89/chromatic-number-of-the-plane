#!/bin/bash
# Re-split one leaf that timed out: its formula (CASE.cnf plus the leaf's literals as unit clauses) becomes a case of
# its own, cut again by cuber2.py (deeper) and certified leaf by leaf.  The sub-cubes cover the leaf (cover checked
# by cuber2.py), so together with the parent tree they still cover every assignment.
# usage: resplit.sh CASE.cnf CUBES.icnf LEAF_INDEX DEPTH VARFILE JOBS [TIME]
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
cnf=$1; cubes=$2; i=$3; depth=$4; vars=$5; jobs=$6; t=${7:-1200}
sub=${cnf%.cnf}_leaf$i.cnf
nice -n 19 python3 - "$cnf" "$cubes" "$i" "$sub" << 'EOF'
import sys
cnf, cubes, i, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
line = [l for l in open(cubes) if l.startswith("a ") or l.startswith("c closed")][i]
lits = line.split()[1:-1] if line.startswith("a ") else line.split()[2:-1]
head, body = open(cnf).read().split("\n", 1)
nv, nc = map(int, head.split()[2:4])
open(out, "w").write(f"p cnf {nv} {nc + len(lits)}\n" + body + "".join(f"{l} 0\n" for l in lits))
EOF
$HERE/cnc_case.sh "$sub" "$depth" "$vars" "$jobs" "$t"
