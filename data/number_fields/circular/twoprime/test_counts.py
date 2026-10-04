"""Cross-check: the number of words (components) per level must equal the counts of levels3.py (snr.py) at
r = 0.3001611 (paper: 13, 33, 57, 81, 97, 109, 125, 133, 105, 89, 73, 49, 33, 25, 17, 9, 5, 5)."""
import sys, time
from gdyn import *
from collections import Counter
r = Fr(sys.argv[1]); K = int(sys.argv[2])
s = Fr(1, 2) - r
comps = level0(s)
for k in range(1, K + 1):
    t0 = time.time()
    comps = lift(comps, s)
    cnt = Counter(word_type(C.word) for C in comps)
    print(f"level {k}: {len(comps)} comps {dict(cnt)} ({time.time()-t0:.1f}s)", flush=True)
