#!/bin/bash
# reshrink.sh: min3multi.py on each grown graph (random deletion orders), then certify_q.py on the result
SC=${SC:-$(cd "$(dirname "$0")/.." && pwd)}
cd $SC/q47
for spec in "11 unsat_q11d30 grow3_q11d30" "191 unsat_q191d240 grow3_q191d240" "431 unsat_q431d600r grow3_q431d600r" \
            "179 unsat_q179d390r grow3_q179d390r" "239 unsat_q239d480r grow3_q239d480r" "59 unsat_q59d210 grow3_q59d210" \
            "119 unsat_q119d240 grow3_q119d240" "131 unsat_q131d390r grow3_q131d390r" "35 unsat_q35d390r grow3_q35d390r" \
            "71 unsat_q71d120 grow3_q71d120" "23 unsat_q23d120 grow3_q23d120" "359 unsat_q359d600r grow3_q359d600r" \
            "47 unsat9139 grow3_d240"; do
  set -- $spec; d=$1; dir=$2; pre=$3
  mkdir -p re$d; cd re$d
  python3 ../min3multi.py ../$dir/${pre}_state.json ../$dir/${pre}_unsat_pts.npy crit${d}m 400 300 > multi.out 2>&1
  python3 ../certify_q.py $d crit${d}m.json cert${d}m > cert${d}m.out 2>&1
  echo "$(date -u +%H:%M) d=$d: $(grep 'vertex-critical' multi.out | tail -1) / $(grep -E 'CERTIFIED|FAIL' cert${d}m.out | tail -1)" >> ../reshrink.summary
  cd ..
done
echo "$(date -u +%H:%M) reshrink finished" >> reshrink.summary
