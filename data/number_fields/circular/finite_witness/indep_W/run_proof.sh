#!/bin/bash
# Usage: run_proof.sh NAME [SEED]   (expects cnf/NAME.ref.cnf; writes proofs/NAME.* and logs/NAME.*.log)
# kissat -> DRAT proof, drat-trim -> verification + LRAT, cake_lpr -> LRAT check.
# The solvers are taken from KISSAT, DRAT_TRIM and CAKE_LPR (default: on the PATH).
set -u
KISSAT=${KISSAT:-kissat}
DRAT=${DRAT_TRIM:-drat-trim}
CAKE=${CAKE_LPR:-cake_lpr}
N=$1
SEED=${2:-7}
mkdir -p proofs logs
CNF=cnf/$N.ref.cnf
PR=proofs/$N.drat
LR=proofs/$N.lrat
nice -n 5 $KISSAT --seed=$SEED $CNF $PR -f > logs/$N.kissat.log 2>&1
echo "kissat exit $?" >> logs/$N.kissat.log
nice -n 5 $DRAT $CNF $PR -L $LR -t 40000 > logs/$N.drat.log 2>&1
echo "drat-trim exit $?" >> logs/$N.drat.log
nice -n 5 $CAKE $CNF $LR > logs/$N.cake.log 2>&1
echo "cake_lpr exit $?" >> logs/$N.cake.log
echo "DONE $N"
