#!/bin/bash
# usage: scan_r5.sh NAME "d:D d:D ..." -- like scan_r2.sh (grow3r.py), but a graph with no 3-colouring is queued (shrinkq.txt)
# for shrink_worker.sh, so the scan goes on. RAD, TMO, MAXN, SOFT, RECOL from the environment.
SC=${SC:-$(cd "$(dirname "$0")/.." && pwd)}
cd "$(dirname "$0")"
name=$1
for pair in $2; do
  d=${pair%%:*}; D=${pair##*:}; t=q${d}d${D}r
  grep -q " $d " done_fields.txt 2>/dev/null && continue
  [ -n "$RESUME" ] || rm -f grow3_${t}_state.json
  QD=$d SOFT=${SOFT:-8} RECOL=${RECOL:-4} timeout ${TMO:-1500} python3 grow3r.py $t ${MAXN:-40000} 900 $D ${RAD:-2} 300 > grow3r_$t.out 2>&1
  if grep -q "NOT 3-COLOURABLE" grow3r_$t.out; then
    mkdir -p unsat_$t && cp grow3_${t}_state.json grow3_${t}_unsat_pts.npy unsat_$t/
    echo "$(date -u +%H:%M) $t: $(grep 'NOT 3-COLOURABLE' grow3r_$t.out | tail -1)" >> scan_r_$name.summary
    echo "$d $t" >> shrinkq.txt; echo " $d " >> done_fields.txt
  else
    echo "$(date -u +%H:%M) $t: $(tail -1 grow3r_$t.out | cut -c1-150)" >> scan_r_$name.summary
  fi
done
echo "$(date -u +%H:%M) scan finished" >> scan_r_$name.summary
