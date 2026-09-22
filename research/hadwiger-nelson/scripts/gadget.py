"""Search the space of LEMMAS, not the space of graphs.

Every search in this project, and every published one as far as it can tell,
grows a graph from a seed.  De Grey did not: he found a seven-point gadget
whose colour palette is capped -- the centre and its radius-2 ring take at
most two of four colours -- and the whole construction is that lemma plus a
way to make copies conflict.  He found it by insight, in 2018.  This project
measured that the cap is abundant at four colours and ABSENT at five, which
is where his engine runs out of fuel.

But that measurement ranged over his shapes: centres and rings.  Nobody has
enumerated small configurations and asked each one for its cap.  That is a
search over lemmas rather than over graphs, and it is cheap, because the
question "is S capped below k" needs exactly one solver call: assert a
proper k-colouring of the carrier AND that every one of the k colours
appears somewhere on S.  UNSAT means no colouring ever spreads k colours
across S -- a cap, in de Grey's sense, found rather than guessed.

The control comes first and decides whether any of it is believable: run at
k = 4 and the search must rediscover his gadget unaided.  If it cannot find
the lemma that is known to be there, its silence at k = 5 means nothing.
"""
import sys, time, math, pickle, itertools
from collections import defaultdict
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_G, build_Sa, build_Y
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

SC = ("/tmp/claude-0/-home-user-darwin-50/"
      "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/")
CAR = sys.argv[1] if len(sys.argv) > 1 else "Sa"
KC = int(sys.argv[2]) if len(sys.argv) > 2 else 4
LIM = int(sys.argv[3]) if len(sys.argv) > 3 else 4000
t0 = time.time()
P = ({"G": lambda k: build_G(k, as_graph=False), "Sa": build_Sa,
      "Y": build_Y}[CAR](K) if CAR in ("G", "Sa", "Y")
     else pickle.load(open(SC + CAR, "rb")))
b = IntBasis.covering(P)
rows = b.rows(P)
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, rows)))
n = len(P)
xy = [(float(t.x), float(t.y)) for t in P]
adj = defaultdict(set)
for a, c in E:
    adj[a].add(c)
    adj[c].add(a)
cls = [[1 + v * KC + c for c in range(KC)] for v in range(n)]
for a, c in E:
    for col in range(KC):
        cls.append([-(1 + a * KC + col), -(1 + c * KC + col)])
s = Solver(name="cd15", bootstrap_with=cls)
assert s.solve(), f"{CAR} is not {KC}-colourable"
print(f"{CAR}: {n} pts, {len(E)} edges, k = {KC}"
      f"  [{time.time()-t0:.0f}s]", flush=True)

SEL = n * KC + 1
nsel = [0]


def trivial(S):
    """A set with a common neighbour outside it is capped for free.

    If some w not in S is adjacent to every point of S then w needs a colour
    none of them has, so S shows at most k-1 colours whatever the geometry.
    A neighbourhood N(v) is the standard case -- v itself is that w -- and
    the first version of this search reported them as discoveries.  They are
    the definition of a proper colouring, not a lemma.
    """
    Ss = set(S)
    cands = set(adj[S[0]]) - Ss
    for v in S[1:]:
        cands &= adj[v]
        if not cands:
            return False
    return bool(cands)


def palette(S):
    """The largest number of colours any proper colouring puts on S."""
    S = sorted(set(S))
    lo, hi = 1, min(len(S), KC)
    while lo < hi:
        mid = (lo + hi + 1) // 2
        nsel[0] += 1
        sel = SEL + nsel[0]
        # at least mid colours appear: pick mid of them and demand each
        # shows up.  Over all choices that is expensive, so ask the
        # equivalent question with counting variables on the colours used.
        used = [SEL + 100000 + nsel[0] * KC + c for c in range(KC)]
        for c in range(KC):
            s.add_clause([-sel, -used[c]] + [1 + v * KC + c for v in S])
        from pysat.card import CardEnc, EncType
        from pysat.formula import IDPool
        pool = IDPool(start_from=max(used) + 1)
        for cl in CardEnc.atleast(lits=used, bound=mid, vpool=pool,
                                  encoding=EncType.seqcounter).clauses:
            s.add_clause([-sel] + list(cl))
        if s.solve(assumptions=[sel]):
            lo = mid
        else:
            hi = mid - 1
    return lo


# The families.  Rings are de Grey's own shape and must be included or the
# control cannot succeed; the rest are what nobody has swept.
d2 = defaultdict(list)
for v in range(n):
    for u in range(n):
        if u != v:
            r2 = round((xy[u][0] - xy[v][0]) ** 2
                       + (xy[u][1] - xy[v][1]) ** 2, 9)
            if r2 <= 16.0:
                d2[(v, r2)].append(u)
fams = []
for (v, r2), ring in d2.items():
    if 3 <= len(ring) <= 14:
        fams.append((f"ring r2={r2:.3f} about {v}", ring))
        fams.append((f"centre+ring r2={r2:.3f} about {v}", [v] + ring))
for v in range(n):
    if KC <= len(adj[v]) <= 16:
        fams.append((f"N({v})", sorted(adj[v])))
        fams.append((f"v{v}+N({v})", [v] + sorted(adj[v])))
for a, c in E:
    common = sorted(adj[a] & adj[c])
    if len(common) >= KC:
        fams.append((f"N({a}) cap N({c})", common))
print(f"{len(fams)} candidate configurations; testing {min(LIM,len(fams))}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
hits = []
skipped = 0
for i, (name, S) in enumerate(fams[:LIM]):
    S = sorted(set(S))
    if len(S) < KC or trivial(S):
        skipped += 1
        continue
    pal = palette(S)
    if pal < KC:
        hits.append((pal, len(S), name, S))
        print(f"   CAP {pal} of {KC}: |S|={len(S):3d}  {name}"
              f"  [{time.time()-t0:.0f}s]", flush=True)
    if i % 500 == 499:
        print(f"   {i+1} seen, {skipped} trivial or too small, "
              f"{len(hits)} caps  [{time.time()-t0:.0f}s]", flush=True)
print(f"\n{len(hits)} capped configurations out of {min(LIM,len(fams))}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
if hits:
    hits.sort()
    for pal, sz, name, S in hits[:12]:
        print(f"   cap {pal} of {KC}, |S|={sz:3d}  {name}", flush=True)
    pickle.dump([(nm, S) for _, _, nm, S in hits],
                open(SC + f"caps_{CAR}_{KC}.pkl", "wb"))
print("DONE", flush=True)
