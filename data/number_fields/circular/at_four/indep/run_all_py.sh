#!/bin/bash
cd "$(dirname "$0")"
for g in 2 3 4 5 6 7 8 9 10 "2 2" "2 4" "2 2 2" "3 3" 11 12 "2 6"; do
  f=out_py/Z_$(echo $g | tr ' ' '_').txt
  s=$(date +%s)
  timeout 1150 nice -n 10 python3 at4check_py.py $g > $f 2>&1
  rc=$?
  e=$(date +%s)
  echo "$g rc=$rc $((e-s))s $(tail -1 $f)"
done
