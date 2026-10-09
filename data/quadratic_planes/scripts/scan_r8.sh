#!/bin/bash
# usage: scan_r8.sh NAME "d:D:RAD:SOFT ..." -- growth (grow3r.py) with a radius and a number of soft rounds per field;
# skips a field whose grow3r_<tag>.out exists (claimed by another scan) or that done_Q14.txt lists, and lists every
# field it finishes there. A graph with no 3-colouring goes to shrinkq.txt for shrink_worker.sh. TMO, MAXN, RECOL
# from the environment.
SC=${SC:-$(cd "$(dirname "$0")/.." && pwd)}
cd "$(dirname "$0")"
name=$1
for spec in $2; do
  IFS=: read d D rad soft <<< "$spec"; t=q${d}d${D}r
  grep -q " $d:$D " done_Q14.txt 2>/dev/null && continue
  [ -e grow3r_$t.out ] && continue
  rm -f grow3_${t}_state.json
  QD=$d SOFT=$soft RECOL=${RECOL:-4} timeout ${TMO:-1500} python3 grow3r.py $t ${MAXN:-50000} 900 $D $rad 300 > grow3r_$t.out 2>&1
  if grep -q "NOT 3-COLOURABLE" grow3r_$t.out; then
    mkdir -p unsat_$t && cp grow3_${t}_state.json grow3_${t}_unsat_pts.npy unsat_$t/
    echo "$(date -u +%H:%M) $t: $(grep 'NOT 3-COLOURABLE' grow3r_$t.out | tail -1)" >> scan_r_$name.summary
    echo "$d $t" >> shrinkq.txt
  else
    echo "$(date -u +%H:%M) $t: $(tail -1 grow3r_$t.out | cut -c1-150)" >> scan_r_$name.summary
  fi
  echo " $d:$D " >> done_Q14.txt
done
echo "$(date -u +%H:%M) scan finished" >> scan_r_$name.summary
