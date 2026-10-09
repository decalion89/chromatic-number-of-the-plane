#!/bin/bash
# Shrink and certify the graphs queued in shrinkq.txt ("d tag" per line), one at a time, into incG_<tag>/; a line
# "<tag> certificate: ..." goes to shrinkq.log. Stops when the queue is done and scan NAME ($1) has finished.
SC=${SC:-$(cd "$(dirname "$0")/.." && pwd)}
cd "$(dirname "$0")"
while true; do
  next=$(while read d t; do grep -q "^$t " shrinkq.log 2>/dev/null || { echo "$d $t"; break; }; done < <(cat shrinkq.txt 2>/dev/null))
  if [ -z "$next" ]; then
    grep -q "scan finished" scan_r_$1.summary 2>/dev/null && break
    sleep 60; continue
  fi
  set -- $next "$1"; d=$1; t=$2; name=$3; set -- "$name"
  mkdir -p incG_$t && cd incG_$t
  python3 ../min3inc.py ../unsat_$t/grow3_${t}_state.json ../unsat_$t/grow3_${t}_unsat_pts.npy crit$d > min3inc.out 2>&1 \
    && python3 ../certify_q.py $d crit$d.json cert$d > cert$d.out 2>&1
  cd ..
  echo "$t certificate: $(grep -E 'CERTIFIED|Error|assert' incG_$t/cert$d.out 2>/dev/null | tail -1) $(tail -n 2 incG_$t/min3inc.out | head -1)" >> shrinkq.log
done
echo "worker finished $(date -u +%H:%M)" >> shrinkq.log
