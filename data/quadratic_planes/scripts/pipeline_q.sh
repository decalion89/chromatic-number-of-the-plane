#!/bin/bash
# Full pipeline for one non-3-colourable instance: exact check, DRAT check, minimisation, certification.
# usage: pipeline_q.sh d tag   (files in SC/q47/unsat_<tag>/)
SC=${SC:-$(cd "$(dirname "$0")/.." && pwd)}
d=$1; t=$2; cd $SC/q47/unsat_$t || exit 1
python3 ../exactcheck.py $d grow3_${t}_state.json grow3_${t}_unsat_pts.npy $t.cnf > exact.log 2>&1
$SC/kissat/build/kissat --time=7200 $t.cnf /dev/shm/$t.drat > $t.kissat.log 2>&1
$SC/drat-trim/drat-trim $t.cnf /dev/shm/$t.drat -t 20000 > $t.drat-trim.log 2>&1; rm -f /dev/shm/$t.drat
python3 ../min3.py grow3_${t}_state.json grow3_${t}_unsat_pts.npy crit$d > min3.out 2>&1
python3 ../certify_q.py $d crit$d.json cert$d > cert$d.out 2>&1
echo "pipeline done" >> cert$d.out
