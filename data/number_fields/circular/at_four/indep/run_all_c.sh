#!/bin/bash
# run the exhaustive C checker on all abelian groups of order <= 16 (and some of order 18, 20)
cd "$(dirname "$0")"
for g in 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 "2 2" "2 4" "2 2 2" "3 3" "2 6" "2 8" "4 4" "2 2 4" "2 2 2 2" 18 "3 6" 20 "2 10"; do
  f=out_c/Z_$(echo $g | tr ' ' '_').txt
  s=$(date +%s)
  timeout 1100 nice -n 10 ./at4check $g > $f 2>&1
  rc=$?
  e=$(date +%s)
  echo "$g rc=$rc $((e-s))s $(tail -1 $f)"
done
