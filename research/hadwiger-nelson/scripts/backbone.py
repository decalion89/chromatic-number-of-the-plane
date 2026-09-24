"""Backbone of the 5-colourings of a grown graph: is ANY pair already forced?

Any non-adjacent pair forced to the SAME colour in every 5-colouring proves
chi >= 6 (chain translated copies until the gap is >= 1/2, then close with
the spindle rotation whose chord at that radius is 1).  Any pair at distance
2 forced APART proves chi >= 6 through Exoo-Ismailescu.  The growth scripts
only ask about one chosen pair; this asks about all of them at once.

Filter: SAT colourings (random vertex order), then Kempe analysis -- a
same-coloured pair survives only if, for every other colour, both ends sit in
one Kempe component (else swapping one component splits it); a distance-2
pair with colours a != b survives only if both ends sit in one (a,b)
component (else a swap makes them alike).  Survivors go to an incremental SAT
query "some survivor breaks"; UNSAT means every survivor is forced.

usage: backbone.py <graph.json> [rounds]
"""
import sys, json, time, random
import numpy as np
from fractions import Fraction as Fr
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver
t0 = time.time(); K = 5
d = json.load(open(sys.argv[1])); ROUNDS = int(sys.argv[2]) if len(sys.argv) > 2 else 12
F = Field(tuple(d["field_generators"]))
V = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
G = build_graph(V); n = G.n; E = list(G.edges())
adj = [set() for _ in range(n)]
for a, b in E: adj[a].add(b); adj[b].add(a)
fx = np.array([p.fx for p in V]); fy = np.array([p.fy for p in V])
print(f"{sys.argv[1]}: n={n} m={len(E)}   [{time.time()-t0:.0f}s]", flush=True)
tri = next(((a, b, c) for a, b in E for c in sorted(adj[a] & adj[b])), None)
X = lambda v, c: 1 + v * K + c
base = [[X(v, c) for c in range(K)] for v in range(n)]
for a, b in E:
    for c in range(K): base.append([-X(a, c), -X(b, c)])
if tri:
    for v, c in zip(tri, range(3)): base.append([X(v, c)])
def decode(model):
    col = np.full(n, -1, dtype=np.int64)
    for v in range(n):
        for c in range(K):
            if model[X(v, c) - 1] > 0: col[v] = c; break
    return col
def random_colouring(seed):
    rng = random.Random(seed); perm = list(range(n)); rng.shuffle(perm); cp = [rng.sample(range(K), K) for _ in range(n)]
    # rename vertex v -> perm[v], colour c -> cp[v][c]; the pinned triangle keeps its colours (renamed consistently)
    Y = lambda v, c: 1 + perm[v] * K + cp[v][c]
    cl = [[Y(v, c) for c in range(K)] for v in range(n)]
    for a, b in E:
        for c in range(K): cl.append([-Y(a, c), -Y(b, c)])
    if tri:
        for v, c in zip(tri, range(3)): cl.append([Y(v, c)])
    rng.shuffle(cl)
    s = Solver(name="cd19", bootstrap_with=cl); assert s.solve(); m = s.get_model(); s.delete()
    col = np.full(n, -1, dtype=np.int64)
    for v in range(n):
        for c in range(K):
            if m[Y(v, c) - 1] > 0: col[v] = c; break
    return col
def kempe_labels(col):
    """lab[(a,b)][v] = component id of v in the subgraph induced by colours a, b (-1 if v has neither)."""
    lab = {}
    for a in range(K):
        for b in range(a + 1, K):
            L = np.full(n, -1, dtype=np.int64); cid = 0
            for s in range(n):
                if col[s] not in (a, b) or L[s] >= 0: continue
                L[s] = cid; st = [s]
                while st:
                    v = st.pop()
                    for w in adj[v]:
                        if L[w] < 0 and col[w] in (a, b): L[w] = cid; st.append(w)
                cid += 1
            lab[(a, b)] = lab[(b, a)] = L
    return lab
def kempe_swap(col, rng):
    a, b = rng.sample(range(K), 2); s = rng.randrange(n)
    if col[s] not in (a, b): return col
    col = col.copy(); st = [s]; seen = {s}
    while st:
        v = st.pop()
        for w in adj[v]:
            if w not in seen and col[w] in (a, b): seen.add(w); st.append(w)
    for v in seen: col[v] = b if col[v] == a else a
    return col
# ---- candidates from the first colouring
col = random_colouring(0)
d2 = lambda I, J: (fx[I] - fx[J]) ** 2 + (fy[I] - fy[J]) ** 2
SI, SJ = [], []
for c in range(K):
    cls = np.nonzero(col == c)[0]
    I, J = np.triu_indices(len(cls), 1); SI.append(cls[I]); SJ.append(cls[J])
SI = np.concatenate(SI); SJ = np.concatenate(SJ)      # same-coloured pairs (never adjacent)
AI, AJ = [], []
for a in range(n):
    dd = (fx - fx[a]) ** 2 + (fy - fy[a]) ** 2
    for b in np.nonzero(np.abs(dd - 4.0) < 1e-9)[0]:
        if b > a and V[a].dist2(V[b]) == 4: AI.append(a); AJ.append(b)
AI = np.array(AI, dtype=np.int64); AJ = np.array(AJ, dtype=np.int64)
print(f"  same candidates {len(SI)}; distance-2 pairs {len(AI)}   [{time.time()-t0:.0f}s]", flush=True)
def filt(col):
    global SI, SJ, AI, AJ
    keep = col[SI] == col[SJ]; SI, SJ = SI[keep], SJ[keep]
    lab = kempe_labels(col)
    for b in range(K):
        if len(SI) == 0: break
        a = col[SI]; m = a != b
        idx = np.nonzero(m)[0]
        if len(idx) == 0: continue
        # pairs coloured a (a != b): they survive only if same (a,b)-component
        la = np.array([lab[(int(x), b)][i] for x, i in zip(a[idx], SI[idx])]) if len(idx) < 50000 else None
        if la is None:
            la = np.empty(len(idx), dtype=np.int64); lb = np.empty(len(idx), dtype=np.int64)
            for aa in range(K):
                if aa == b: continue
                sel = a[idx] == aa; L = lab[(aa, b)]
                la[sel] = L[SI[idx][sel]]; lb[sel] = L[SJ[idx][sel]]
        else:
            lb = np.array([lab[(int(x), b)][j] for x, j in zip(a[idx], SJ[idx])])
        dead = np.zeros(len(SI), dtype=bool); dead[idx[la != lb]] = True
        SI, SJ = SI[~dead], SJ[~dead]
    if len(AI):
        ca, cb = col[AI], col[AJ]
        alike = ca == cb
        la = np.array([lab[(int(x), int(y))][i] if x != y else 0 for x, y, i in zip(ca, cb, AI)])
        lb = np.array([lab[(int(x), int(y))][j] if x != y else 1 for x, y, j in zip(ca, cb, AJ)])
        keep = (~alike) & (la == lb); AI, AJ = AI[keep], AJ[keep]
filt(col)
print(f"  after colouring 0 + Kempe: same {len(SI)}, apart {len(AI)}   [{time.time()-t0:.0f}s]", flush=True)
rng = random.Random(1)
for r in range(1, ROUNDS + 1):
    col = random_colouring(r); filt(col)
    for _ in range(4):
        for _ in range(200): col = kempe_swap(col, rng)
        filt(col)
    print(f"  round {r}: same {len(SI)}, apart {len(AI)}   [{time.time()-t0:.0f}s]", flush=True)
    if len(SI) == 0 and len(AI) == 0: break
# the survivors themselves, for gadget seeds (MODE=same on an apart survivor)
json.dump({"same": [[int(a), int(b)] for a, b in zip(SI, SJ)], "apart": [[int(a), int(b)] for a, b in zip(AI, AJ)]},
          open(sys.argv[1].replace(".json", "") + ".survivors.json", "w"))
# ---- exact phase: incremental "some survivor breaks"
def exact(kind, I, J):
    I, J = list(map(int, I)), list(map(int, J))
    s = Solver(name="cd19", bootstrap_with=base); top = n * K + 1; alive = set(range(len(I))); calls = 0
    sel = {}
    for p in range(len(I)):
        if kind == "same":                   # s_p -> colour sets of I[p], J[p] disjoint
            sp = top; top += 1; sel[p] = sp
            for c in range(K): s.add_clause([-sp, -X(I[p], c), -X(J[p], c)])
        else:                                 # s_p -> some colour shared
            sp = top; top += 1; sel[p] = sp; zs = []
            for c in range(K):
                z = top; top += 1; zs.append(z)
                s.add_clause([-z, X(I[p], c)]); s.add_clause([-z, X(J[p], c)])
            s.add_clause([-sp] + zs)
    while alive:
        act = top; top += 1; calls += 1
        s.add_clause([-act] + [sel[p] for p in alive])
        if not s.solve(assumptions=[act]): break
        col = decode(s.get_model())
        if kind == "same": alive = {p for p in alive if col[I[p]] == col[J[p]]}
        else: alive = {p for p in alive if col[I[p]] != col[J[p]]}
    s.delete()
    return [(I[p], J[p]) for p in sorted(alive)], calls
for kind, I, J in (("same", SI, SJ), ("apart", AI, AJ)):
    if len(I) == 0: continue
    if len(I) > 20000: print(f"  {kind}: {len(I)} survivors, too many for the exact phase", flush=True); continue
    forced, calls = exact(kind, I, J)
    print(f"  EXACT {kind}: {len(I)} survivors -> {len(forced)} FORCED after {calls} SAT calls   [{time.time()-t0:.0f}s]", flush=True)
    for a, b in forced[:20]:
        print(f"    FORCED {kind}: {a} {b}  d^2 = {V[a].dist2(V[b])}  (~{float((fx[a]-fx[b])**2+(fy[a]-fy[b])**2):.4f})", flush=True)
