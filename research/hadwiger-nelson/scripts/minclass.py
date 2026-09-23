"""Shrink the carrier while keeping the whole forced class of three.

The one-pivot object -- three rotations about a, all the force of the class,
no group -- is 12097 points because the carrier is 3025.  Its five-colour test
is the only one of these that has a chance of finishing, and halving the
carrier halves it again.

The class is kept by keeping BOTH forcings, (a,b) and (a,c); the third, (b,c),
follows by transitivity.  Each is the cheap kind of instance -- "colour this
subgraph with four colours, with c(a)=0 and c(b)=1" -- so a deletion costs two
fast UNSAT proofs rather than one slow colourability question.
"""
import sys, time, random, json
from fractions import Fraction as Fr
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_Sa
from hn.geometry import Point, Rotation, rotation_joining, _rot60
from hn.field import Field
from hn.graph import build_graph
from pysat.solvers import Solver

SEED = int(sys.argv[1]) if len(sys.argv) > 1 else 1
ROUNDS = int(sys.argv[2]) if len(sys.argv) > 2 else 3
F = Field((3, 11, 23, 247))
rot60 = _rot60(F); Sa = build_Sa(F)
def orb(p, refl=False):
    out, q = [], p
    for _ in range(6):
        if q not in out: out.append(q)
        q = rot60(q)
    if refl:
        for q in list(out):
            r = Point(q.x, -q.y)
            if r not in out: out.append(r)
    return out
t0 = time.time()
seen, U = set(Sa), list(Sa)
for w in orb(Sa[199], True):
    rot = rotation_joining(Fr(1), F).about(w)
    for p in Sa:
        q = rot(p)
        if q not in seen: seen.add(q); U.append(q)
g = build_graph(U); n = g.n
A_, B_, C_ = 726, 1526, 2730
print(f"carrier n={n}; class {A_},{B_},{C_}   [{time.time()-t0:.0f}s]", flush=True)
K = 4
X = lambda v, c: 1 + v * K + c
SEL = lambda v: 1 + n * K + v
base = [[-SEL(v)] + [X(v, c) for c in range(K)] for v in range(n)]
for u, v in g.edges():
    for c in range(K):
        base.append([-X(u, c), -X(v, c)])
cnfs = []
for (x, y) in ((A_, B_), (A_, C_)):
    cnfs.append(list(base) + [[X(x, 0)], [X(y, 1)]])
def keeps(sub):
    a = sorted(sub)
    for cnf in cnfs:
        s = Solver(name="m22", bootstrap_with=cnf)
        ok = s.solve(assumptions=[SEL(v) for v in a]); s.delete()
        if ok: return False
    return True
# Start from the UNSAT cores, not from the whole carrier.  Each forcing
# instance hands back a subset that already refuses, and the union of the two
# refuses both -- a large first cut for two solver calls instead of 3025 times
# two.
def core_of(cnf, sub):
    s = Solver(name="m22", bootstrap_with=cnf)
    ok = s.solve(assumptions=[SEL(v) for v in sorted(sub)])
    c = s.get_core(); s.delete()
    if ok: return None
    return set(l - 1 - n * K for l in c) if c else set(sub)
cur = set(range(n))
assert keeps(cur), "the class is not forced in the full carrier"
keep3 = {A_, B_, C_}
for it in range(12):
    nxt = set()
    for cnf in cnfs:
        c = core_of(cnf, cur)
        if c is None: break
        nxt |= c
    else:
        nxt |= keep3
        if len(nxt) >= len(cur): break
        cur = nxt
        print(f"  cores: {len(cur)}   [{time.time()-t0:.0f}s]", flush=True)
        continue
    break
print(f"  after cores: {len(cur)}   [{time.time()-t0:.0f}s]", flush=True)
rng = random.Random(SEED)
for rd in range(ROUNDS):
    order = list(cur); rng.shuffle(order)
    dropped = 0
    for v in order:
        if v in (A_, B_, C_) or v not in cur: continue
        if keeps(cur - {v}):
            cur.discard(v); dropped += 1
    print(f"  round {rd}: dropped {dropped} -> carrier {len(cur)}"
          f"   [{time.time()-t0:.0f}s]", flush=True)
    if dropped == 0: break
V = [U[v] for v in sorted(cur)]
idx = {v: i for i, v in enumerate(sorted(cur))}
gv = build_graph(V)
print(f"  minimal carrier: n={gv.n} m={sum(len(a) for a in gv.adj)//2}"
      f"   [{time.time()-t0:.0f}s]", flush=True)
a2, b2, c2 = idx[A_], idx[B_], idx[C_]
u = V[b2] - V[a2]; v = V[c2] - V[a2]
P = u.x * v.x + u.y * v.y
Q = u.x * v.y - u.y * v.x
R = F.rational((Fr(64, 3) + Fr(256, 9) - 1) / 2)
D = F.sqrt(23) * F.rational(Fr(13, 18))
inv = F.rational(Fr(27, 16384))
rots = []
for sgn in (1, -1):
    co = (P * R + Q * D * F.rational(sgn)) * inv
    si = (Q * R - P * D * F.rational(sgn)) * inv
    if co * co + si * si == F.rational(1):
        r = Rotation(co, si)
        if (r.about(V[a2])(V[b2]) - V[c2]).norm2() == F.rational(1):
            rots.append(r)
rots.append(rotation_joining(Fr(64, 9), F))
seen2, W = set(V), list(V)
for r in rots:
    rot = r.about(V[a2])
    for p in V:
        z = rot(p)
        if z not in seen2: seen2.add(z); W.append(z)
gw = build_graph(W); nw = gw.n; mw = sum(len(x) for x in gw.adj) // 2
print(f"  one-pivot object: n={nw} m={mw} deg={2.0*mw/nw:.2f}"
      f"   [{time.time()-t0:.0f}s]", flush=True)
tri = gw.find_clique(3)
for KK in (4, 5, 6):
    XX = lambda v, c: 1 + v * KK + c
    cc = [[XX(v, c) for c in range(KK)] for v in range(nw)]
    for x, y in gw.edges():
        for c in range(KK):
            cc.append([-XX(x, c), -XX(y, c)])
    for i, x in enumerate(tri):
        cc.append([XX(x, i)])
        for c in range(KK):
            if c != i: cc.append([-XX(x, c)])
    t1 = time.time()
    sv = Solver(name="cd19", bootstrap_with=cc); ok = sv.solve(); sv.delete()
    print(f"  {KK}-colourable: {ok}   [{time.time()-t1:.0f}s]", flush=True)
    if ok:
        if KK >= 6: print("  *** chi = 6 ***", flush=True)
        break
    print(f"  *** refuses {KK} ***", flush=True)
    if KK == 5:
        json.dump({"field_generators": list(F.gens), "n": nw, "m": mw,
                   "points": [[[[c.numerator, c.denominator] for c in p.x.c],
                               [[c.numerator, c.denominator] for c in p.y.c]]
                              for p in gw.vertices]},
                  open("/home/user/darwin-50/research/hadwiger-nelson/data/six_candidate.json", "w"))
        print("  *** written data/six_candidate.json ***", flush=True)
