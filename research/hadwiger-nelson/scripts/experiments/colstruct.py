"""Is the local-search colouring of a grown graph structured?  Recolour it with the tabu search,
put every vertex in module coordinates, and look for periods: difference vectors w with
c(x) = c(x + w) for every pair of vertices that differ by w (and many such pairs)."""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, json, subprocess, random, time
from fractions import Fraction as Fr
from collections import defaultdict
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
exec(open(HN_DIR + "/scripts/gate.py").read().split("def gate(g, label):")[0])
t0 = time.time()
d = json.load(open(sys.argv[1])); MODE = sys.argv[2] if len(sys.argv) > 2 else "plain"
F = Field(tuple(d["field_generators"]))
mk = lambda xy: Point(F.element([Fr(a, b) for a, b in xy[0]]), F.element([Fr(a, b) for a, b in xy[1]]))
V = [mk(xy) for xy in d["points"]]; U = [mk(xy) for xy in d["units"]]; A, Bi = d.get("A"), d.get("B")
key = lambda x, y: (round(x, 9), round(y, 9))
vkey = {key(p.fx, p.fy): i for i, p in enumerate(V)}
EA, EB = [], []
for i, p in enumerate(V):
    for u in U:
        j = vkey.get(key(p.fx + u.fx, p.fy + u.fy))
        if j is not None and j > i: EA.append(i); EB.append(j)
n = len(V)
ea, eb = EA, EB
if MODE == "apart": ea = [A if x == Bi else x for x in EA]; eb = [A if x == Bi else x for x in EB]
inp = f"{n} {len(ea)} 5 20000000 {int(sys.argv[3]) if len(sys.argv) > 3 else 1}\n" + "\n".join(f"{a} {b}" for a, b in zip(ea, eb)) + "\n" + "\n".join(["-1"] * n) + "\n"
out = subprocess.run([os.path.join(HN_DIR, "scripts", "tabucol")], input=inp, capture_output=True, text=True).stdout.split("\n")
col = [int(x) for x in out[1:1 + n]]
if MODE == "apart": col[Bi] = col[A]
print(f"n={n} m={len(EA)}; tabu: {out[0]}; proper: {all(col[a] != col[b] for a, b in zip(EA, EB))}   [{time.time()-t0:.0f}s]", flush=True)
# module coordinates relative to vertex 0
raw = []
for p in V:
    q = p - V[0]; raw.append(tuple(q.x.c) + tuple(q.y.c))
Ev = [tuple(u.x.c) + tuple(u.y.c) for u in U]
from math import gcd
den = 1
for v in Ev + raw:
    for x in v: den = den * Fr(x).denominator // gcd(den, Fr(x).denominator)
E = sorted({max(w, tuple(-x for x in w)) for w in (tuple(int(Fr(x) * den) for x in v) for v in Ev)})
B = echelon(E); r = len(B)
co = [tuple(coords(B, tuple(int(Fr(x) * den) for x in v))) for v in raw]
print(f"module rank {r}; coordinates done   [{time.time()-t0:.0f}s]", flush=True)
# periods: count agreement per difference vector over sampled pairs
agree = defaultdict(lambda: [0, 0])
rng = random.Random(0)
idx = list(range(n))
for _ in range(3_000_000):
    i, j = rng.randrange(n), rng.randrange(n)
    if i == j: continue
    w = tuple(a - b for a, b in zip(co[i], co[j]))
    if w[next(k for k, t in enumerate(w) if t)] < 0: w = tuple(-t for t in w); i, j = j, i
    e = agree[w]; e[0] += 1; e[1] += col[i] == col[j]
top = sorted(((c, s, w) for w, (c, s) in agree.items() if c >= 8), reverse=True)[:15]
print(f"difference vectors seen >= 8 times: {sum(1 for w,(c,s) in agree.items() if c >= 8)}; perfectly periodic among them: {sum(1 for w,(c,s) in agree.items() if c >= 8 and s == c)}")
for c, s, w in top: print(f"  seen {c:3d}, alike {s:3d}   w = {w}")
# residues: is colour a function of coords mod q?
for q in (2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 15):
    cls = defaultdict(set)
    for i in range(n): cls[tuple(x % q for x in co[i])].add(col[i])
    mono = sum(1 for s in cls.values() if len(s) == 1)
    print(f"  mod {q}M: {len(cls)} residue classes met, {mono} monochromatic")
