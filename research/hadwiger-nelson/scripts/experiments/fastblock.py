import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle
from fractions import Fraction as Fr
from math import gcd
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
src = open("/tmp/hn/"
           "degreyfield.py").read()
exec(src[:src.index("t0 = time.time()\nrh = ")])
from hn.homcol import blocks_at, has_homomorphism, on_lattice
t0 = time.time()
with open("/tmp/hn/kneck.pkl", "rb") as fh:
    ws = pickle.load(fh)
ONE_PLUS = kadd(K1, Z6)


def flat2(x):
    out = []
    for part in x:
        out.extend(part)
    return tuple(out)


pts, B = [kz()], kz()
for w in ws:
    for off in (w, kmul(Z6, w), kmul(ONE_PLUS, w)):
        pts.append(kadd(B, off))
    B = kadd(B, kmul(ONE_PLUS, w))
uniq, seenq = [], set()
for q in pts:
    if q not in seenq:
        seenq.add(q)
        uniq.append(q)
E = [(i, j) for i in range(len(uniq)) for j in range(i + 1, len(uniq))
     if knorm2(ksub(uniq[j], uniq[i])) == K1]
rows = []
for a, b in E:
    rows.append([Fr(q) for q in flat2(ksub(uniq[b], uniq[a]))])
    rows.append([Fr(q) for q in flat2(ksub(uniq[a], uniq[b]))])
dn = 1
for r in rows:
    for q in r:
        dn = dn * q.denominator // gcd(dn, q.denominator)
ints = {tuple(int(q * dn) for q in r) for r in rows}
content = 0
for v in ints:
    for x in v:
        content = gcd(content, abs(x))
if content > 1:
    ints = {tuple(x // content for x in v) for v in ints}
ints = sorted(ints)
print(f"necklace: {len(uniq)} points, {len(E)} edges, {len(ints)} directions"
      f"  [{time.time()-t0:.0f}s]", flush=True)
for n in (2, 3, 4, 5):
    t1 = time.time()
    b = blocks_at(ints, n)
    print(f"  n = {n}: {'blocks' if b else 'a coset colouring exists'}  "
          f"[{time.time()-t1:.1f}s]", flush=True)
red = on_lattice(ints)
sizes = sorted({max(abs(x) for x in v) for v in red})
print(f"  largest lattice coordinate {sizes[-1]}  [{time.time()-t0:.0f}s]",
      flush=True)
