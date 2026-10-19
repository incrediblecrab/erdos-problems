#!/bin/bash
# Enumerate ALL k-APs of primes with last term <= B, for the Hardy-Littlewood
# count test.  B chosen per k so the count is large and the run is affordable.
set -e
cd "$(dirname "$0")/.."
for kb in "8 100000000" "9 100000000" "10 200000000" "11 1000000000" \
          "12 1000000000" "13 3000000000" "14 5000000000" "15 8000000000" \
          "16 8000000000" "17 20000000000" "18 30000000000" "19 150000000000"; do
  set -- $kb
  ./src/apsearch "$1" "$2" 10 2>>results/counts.log > "results/count_k$1.json"
  python3 -c "
import json;d=json.load(open('results/count_k$1.json'))
print(f\"k={d['k']:2d} B={d['B']:>13d} total={d['total']:>10d} min_last={d.get('min_last')} {d['seconds']}s\",flush=True)"
done
echo "COUNT BATCH DONE"
