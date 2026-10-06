#!/bin/bash
# usage: iter_shrink.sh START_NAME ROUNDS : repeat (kissat DRAT, drat-trim core, keep core vertices)
# KISSAT and DRAT_TRIM: the solver and the checker (environment variables; default: the names on the PATH)
KISSAT=${KISSAT:-kissat}; DRAT_TRIM=${DRAT_TRIM:-drat-trim}
cur=$1
for r in $(seq 1 $2); do
  $KISSAT --seed=$r $cur.cnf $cur.drat > $cur.k.out 2>&1
  grep -aq "^s UNSATISFIABLE" $cur.k.out || { echo "round $r: $cur not UNSAT"; exit 3; }
  $DRAT_TRIM $cur.cnf $cur.drat -c $cur.core > $cur.dt.out 2>&1
  grep -aq "s VERIFIED" $cur.dt.out || { echo "round $r: drat-trim did not verify $cur"; exit 4; }
  nxt=${PFX:-it}$r
  python3 shrink_core.py W_$cur.json $cur.core W_$nxt.json
  python3 certify4.py W_$nxt.json $nxt.cnf | tail -1
  rm -f "$cur.drat"
  cur=$nxt
done
echo "last: $cur"
