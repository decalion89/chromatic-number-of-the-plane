#!/bin/bash
# usage: scan_r.sh NAME "d:D d:D ..." -- like scan_q.sh but with grow3r.py (recolouring and soft rounds when stuck),
# maxn 40000 and 300 points per round. One line per run in scan_r_NAME.summary.
SC=${SC:-$(cd "$(dirname "$0")/.." && pwd)}
cd $SC/q47
name=$1
done_d=" 119 191 "
for pair in $2; do
  d=${pair%%:*}; D=${pair##*:}; t=q${d}d${D}r
  case "$done_d" in *" $d "*) continue;; esac
  rm -f grow3_${t}_state.json
  QD=$d SOFT=8 RECOL=4 timeout 1500 python3 grow3r.py $t 40000 900 $D 2 300 > grow3r_$t.out 2>&1
  if grep -q "NOT 3-COLOURABLE" grow3r_$t.out; then
    mkdir -p unsat_$t && cp grow3_${t}_state.json grow3_${t}_unsat_pts.npy unsat_$t/ && cp /dev/shm/g3_$t.cnf unsat_$t/$t.cnf
    echo "$(date -u +%H:%M) $t: $(grep 'NOT 3-COLOURABLE' grow3r_$t.out | tail -1)" >> scan_r_$name.summary
    bash pipeline_q.sh $d $t
    echo "$(date -u +%H:%M) $t pipeline: $(grep -E 'CERTIFIED|FAIL' unsat_$t/cert$d.out | tail -1)" >> scan_r_$name.summary
    done_d="$done_d$d "
  else
    echo "$(date -u +%H:%M) $t: $(tail -1 grow3r_$t.out | cut -c1-150)" >> scan_r_$name.summary
  fi
done
echo "$(date -u +%H:%M) scan finished" >> scan_r_$name.summary
