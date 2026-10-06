#!/bin/bash
# usage: iter_shrink_z.sh START_NAME ROUNDS, in the folder of START_NAME.cnf and W_START_NAME.json: each round refutes the
# current formula with kissat (seed = round), extracts the clausal core with drat-trim and keeps the vertices that occur
# in it (shrink_core.py, certify_z24.py). KISSAT and DRAT_TRIM name the binaries; PFX (default z) names the rounds.
C=$(dirname "$(readlink -f "$0")")
cur=$1
for r in $(seq 1 $2); do
  "$KISSAT" --seed=$r $cur.cnf $cur.drat > $cur.k.out 2>&1
  grep -aq "^s UNSATISFIABLE" $cur.k.out || { echo "round $r: $cur not UNSAT"; exit 3; }
  "$DRAT_TRIM" $cur.cnf $cur.drat -c $cur.core > $cur.dt.out 2>&1
  grep -aq "s VERIFIED" $cur.dt.out || { echo "round $r: drat-trim did not verify $cur"; exit 4; }
  nxt=${PFX:-z}$r
  python3 "$C/shrink_core.py" W_$cur.json $cur.core W_$nxt.json
  python3 "$C/certify_z24.py" W_$nxt.json $nxt.cnf | tail -1
  rm -f -- "$PWD/$cur.drat"
  cur=$nxt
done
echo "last: $cur"
