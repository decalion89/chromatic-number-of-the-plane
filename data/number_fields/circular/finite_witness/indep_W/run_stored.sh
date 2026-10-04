#!/bin/bash
# Usage: run_stored.sh NAME   -- drat-trim (with LRAT output) and cake_lpr on the repository's stored formula and
# DRAT proof, decompressed into stored/ by validate_stored.py.  Solvers from DRAT_TRIM and CAKE_LPR (default: PATH).
set -u
DRAT=${DRAT_TRIM:-drat-trim}
CAKE=${CAKE_LPR:-cake_lpr}
N=$1
mkdir -p logs
nice -n 5 $DRAT stored/$N.stored.cnf stored/$N.stored.drat -L stored/$N.stored.lrat -t 40000 > logs/$N.stored.drat.log 2>&1
echo "drat-trim exit $?" >> logs/$N.stored.drat.log
nice -n 5 $CAKE stored/$N.stored.cnf stored/$N.stored.lrat > logs/$N.stored.cake.log 2>&1
echo "cake_lpr exit $?" >> logs/$N.stored.cake.log
echo DONE $N
