"""Re-check the necklace's blocking on the lattice its directions generate.

This is the one that matters.  Blocking can be bolted onto any graph by taking
a disjoint rotated copy -- it changes no colouring and adds a whole new orbit
of directions -- so blocking on its own is cheap.  What is not cheap is a
4-CRITICAL graph that blocks, because a critical graph has nothing spare to
bolt anything onto.  The necklace is that graph, so its verdict has to survive
the stricter test.

Two things are checked.  The directions are normalised by the GLOBAL content
rather than vector by vector (dividing each vector by its own gcd only makes
blocking harder to achieve, so the old result stands either way, but the
global version is the honest module).  And they are then rewritten in a
Z-basis of the lattice they generate, so that the search for phi runs over the
right group.
"""
import sys, time, pickle
from fractions import Fraction as Fr
from math import gcd
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
src = open("/tmp/claude-0/-home-user-darwin-50/"
           "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/"
           "degreyfield.py").read()
exec(src[:src.index("t0 = time.time()\nrh = ")])
from hn.homcol import has_homomorphism
exec(open("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-"
          "39c59179b415/scratchpad/hnfcheck.py").read()
     .split("print(\"re-checking")[0].split("t0 = time.time()")[1])

t0 = time.time()
with open("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-"
          "39c59179b415/scratchpad/kneck.pkl", "rb") as fh:
    ws = pickle.load(fh)
ONE_PLUS = kadd(K1, Z6)
assert knorm2(ONE_PLUS) == krat(3)
print(f"necklace of {len(ws)} rhombi loaded  [{time.time()-t0:.0f}s]",
      flush=True)


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
k = len(ws)
print(f"  {len(uniq)} points (3k+1 = {3*k+1}), {len(E)} edges "
      f"(5k+1 = {5*k+1})  [{time.time()-t0:.0f}s]", flush=True)

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
glob = sorted(ints)

per = set()
for v in glob:
    g = 0
    for x in v:
        g = gcd(g, abs(x))
    per.add(tuple(x // g for x in v) if g > 1 else v)
per = sorted(per)
print(f"  {len(glob)} directions under global content, {len(per)} under "
      f"per-vector gcd  [{time.time()-t0:.0f}s]", flush=True)

for name, VV in (("global content, Z^d", glob),
                 ("per-vector gcd, Z^d", per)):
    bl = [n for n in (2, 3, 4, 5) if has_homomorphism(VV, n)[0] is None]
    print(f"  {name}: blocks at {bl}  [{time.time()-t0:.0f}s]", flush=True)

red, rank = in_basis(glob, len(glob[0]))
bl = [n for n in (2, 3, 4, 5) if has_homomorphism(red, n)[0] is None]
print(f"  lattice rank {rank}, {len(red)} directions: blocks at {bl}  "
      f"[{time.time()-t0:.0f}s]", flush=True)
