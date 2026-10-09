#!/bin/bash
# usage: shrink_cayley.sh IN.json ROUNDS NAME, in a scratch folder: ROUNDS rounds of clausal cores on the cycles of IN.json
# (kissat refutes the formula, drat-trim -c extracts the core, cayley.py reduce keeps the cycles that occur in it), then
# NAME.json.gz, NAME.cnf.gz and the proof NAME.drat of the last formula, checked by drat-trim. KISSAT and DRAT_TRIM name
# the binaries.
C=$(dirname "$(readlink -f "$0")")
cp "$1" cur.json
for r in $(seq 1 $2); do
  python3 "$C/cayley.py" cnf cur.json cur.cnf
  "$KISSAT" --time=3600 cur.cnf cur.drat > cur.k.out 2>&1
  grep -aq "^s UNSATISFIABLE" cur.k.out || { echo "round $r: not UNSAT"; exit 3; }
  "$DRAT_TRIM" cur.cnf cur.drat -c cur.core -t 20000 > cur.dt.out 2>&1
  grep -aq "s VERIFIED" cur.dt.out || { echo "round $r: drat-trim did not verify"; exit 4; }
  echo "round $r: $(python3 "$C/cayley.py" reduce cur.json cur.core cur.json)"
  rm -f -- "$PWD/cur.drat"
done
python3 "$C/cayley.py" final cur.json "$3"
gunzip -c "$3.cnf.gz" > "$3.cnf"
"$KISSAT" "$3.cnf" "$3.drat" > "$3.k.out" 2>&1
grep -a "^s " "$3.k.out"
"$DRAT_TRIM" "$3.cnf" "$3.drat" -t 20000 | tr '\r' '\n' | grep -a "^s "    # drat-trim ends its lines with \r
