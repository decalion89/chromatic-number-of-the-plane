#!/bin/bash
# usage: scan_q.sh "d:D d:D ..."  -- grow for each pair in order; at the first UNSAT for a field d, run the
# certification pipeline and skip the remaining denominators of that d. One line per run in scan_q.summary.
SC=${SC:-$(cd "$(dirname "$0")/.." && pwd)}
cd $SC/q47
done_d=" "
for pair in $1; do
  d=${pair%%:*}; D=${pair##*:}; t=q${d}d${D}
  case "$done_d" in *" $d "*) continue;; esac
  rm -f grow3_${t}_state.json
  QD=$d timeout 2400 python3 grow3q.py $t 16000 900 $D 2 200 > grow3q_$t.out 2>&1
  if grep -q "NOT 3-COLOURABLE" grow3q_$t.out; then
    mkdir -p unsat_$t && cp grow3_${t}_state.json grow3_${t}_unsat_pts.npy unsat_$t/ && cp /dev/shm/g3_$t.cnf unsat_$t/$t.cnf
    echo "$(date -u +%H:%M) $t: $(grep 'NOT 3-COLOURABLE' grow3q_$t.out | tail -1)" >> scan_q.summary
    bash pipeline_q.sh $d $t
    echo "$(date -u +%H:%M) $t pipeline: $(grep -E 'CERTIFIED|FAIL' unsat_$t/cert$d.out | tail -1)" >> scan_q.summary
    done_d="$done_d$d "
  else
    echo "$(date -u +%H:%M) $t: $(tail -1 grow3q_$t.out | cut -c1-150)" >> scan_q.summary
  fi
done
echo "$(date -u +%H:%M) scan finished" >> scan_q.summary
