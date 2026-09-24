"""Colouring-guided growth with local search first (lean memory).

Each iteration: repair the previous colouring on the enlarged graph with a C
tabu search (c/tabucol); only if that fails, ask the incremental CDCL solver
(seeded with the tabu search's best phases).  Then insert the R richest pool
points whose neighbourhoods show all five colours.  The pool stores only
(base vertex, unit index, neighbours) keyed by float coordinates; exact points
are built only for inserted vertices, and the unit directions are stored in
the checkpoint so a resume needs no exact rebuild.

MODE plain | apart (c(a) = c(b) imposed) | same (c(a) != c(b) imposed).
Any UNSAT is a claim until verify6 / verify_gadget re-check it exactly.

usage: grow_ls2.py <in.json> <out.json> [R] [near_weight]   (env MODE, SOLVER, BUDGET, LSIT)
"""
import sys, time, json, os, subprocess
from fractions import Fraction as Fr
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver
t0 = time.time(); K = 5
IN, OUT = sys.argv[1], sys.argv[2]
R = int(sys.argv[3]) if len(sys.argv) > 3 else 20
WN = float(sys.argv[4]) if len(sys.argv) > 4 else 1.0
MODE = os.environ.get("MODE", "plain"); SOLVER = os.environ.get("SOLVER", "cadical195")
BUDGET = int(os.environ.get("BUDGET", "200000000")); LSIT = int(os.environ.get("LSIT", "3000000"))
_here = os.path.dirname(os.path.abspath(__file__))
TABU = os.environ.get("TABU", os.path.join(_here, "tabucol"))
if not os.path.exists(TABU):          # build the tabu search from scripts/tabucol.c on first use
    subprocess.run(["gcc", "-O2", "-o", TABU, os.path.join(_here, "tabucol.c")], check=True)
d = json.load(open(IN))
F = Field(tuple(d["field_generators"]))
mk = lambda xy: Point(F.element([Fr(a, b) for a, b in xy[0]]), F.element([Fr(a, b) for a, b in xy[1]]))
V = [mk(xy) for xy in d["points"]]
A, Bi = d.get("A"), d.get("B")
key = lambda p: (round(p.fx, 9), round(p.fy, 9))
if "units" in d:
    U = [mk(xy) for xy in d["units"]]
else:
    G = build_graph(V); Ud = {}
    for a, b in G.edges():
        for u in (V[b] - V[a], V[a] - V[b]): Ud.setdefault(key(u), u)
    U = list(Ud.values())
Ufx = [u.fx for u in U]; Ufy = [u.fy for u in U]
Vfx = [p.fx for p in V]; Vfy = [p.fy for p in V]
vkey = {}
for i, p in enumerate(V): vkey.setdefault(key(p), i)
fkey = lambda x, y: (round(x, 9), round(y, 9))
adj = [set() for _ in V]; EA, EB = [], []
for i in range(len(V)):
    for k in range(len(U)):
        j = vkey.get(fkey(Vfx[i] + Ufx[k], Vfy[i] + Ufy[k]))
        if j is not None and j > i and j not in adj[i]:
            adj[i].add(j); adj[j].add(i); EA.append(i); EB.append(j)
pool = {}
def feed(i):
    for k in range(len(U)):
        kk = fkey(Vfx[i] + Ufx[k], Vfy[i] + Ufy[k])
        if kk in vkey: continue
        e = pool.get(kk)
        if e is None: pool[kk] = [i, k, {i}]
        else: e[2].add(i)
for i in range(len(V)): feed(i)
col = []
ser = lambda q: [[[t.numerator, t.denominator] for t in q.x.c], [[t.numerator, t.denominator] for t in q.y.c]]
def save(tag):
    json.dump({"field_generators": list(F.gens), "A": A, "B": Bi, "status": tag, "mode": MODE, "colouring": col[:len(V)] if len(col) == len(V) else None,
               "units": [ser(u) for u in U], "two_edges": TWO, "points": [ser(q) for q in V]}, open(OUT, "w"))
X = lambda v, c: 1 + v * K + c
s = Solver(name=SOLVER)
for v in range(len(V)): s.add_clause([X(v, c) for c in range(K)])
for a, b in zip(EA, EB):
    for c in range(K): s.add_clause([-X(a, c), -X(b, c)])
tri = next(((a, b, c) for a, b in zip(EA, EB) for c in sorted(adj[a] & adj[b])), None)
if tri:
    for k, v in enumerate(tri): s.add_clause([X(v, k)])
if MODE == "apart":
    for c in range(K): s.add_clause([-X(A, c), X(Bi, c)]); s.add_clause([X(A, c), -X(Bi, c)])
elif MODE == "same":
    for c in range(K): s.add_clause([-X(A, c), -X(Bi, c)])
TWO = [tuple(t) for t in d.get("two_edges", [])]
import numpy as _np
_tgt = None
def near_targets(kk):
    if _tgt is None or len(_tgt) == 0: return 0.0
    return float(_np.min(_np.hypot(_tgt[:, 0] - kk[0], _tgt[:, 1] - kk[1])))
if TWO:
    near = near_targets
elif A is not None:
    ax, ay, bx, by = Vfx[A], Vfy[A], Vfx[Bi], Vfy[Bi]
    near = lambda kk: min(((kk[0] - ax) ** 2 + (kk[1] - ay) ** 2) ** 0.5, ((kk[0] - bx) ** 2 + (kk[1] - by) ** 2) ** 0.5)
if not TWO and A is None:
    cx = sum(Vfx) / len(V); cy = sum(Vfy) / len(V)
    near = lambda kk: ((kk[0] - cx) ** 2 + (kk[1] - cy) ** 2) ** 0.5
TABU2 = os.environ.get("TABU2", os.path.join(_here, "tabu2"))
if not os.path.exists(TABU2):
    subprocess.run(["gcc", "-O2", "-o", TABU2, os.path.join(_here, "tabu2.c")], check=True)
IT2 = int(os.environ.get("IT2", "300000"))
def tabucol(n, init, seed):
    """two-phase tabu search: a proper colouring (MODE constraint included), then as few alike
    skeleton two-edges as it can find without breaking a unit edge.  Returns (proper?, #alike, colouring)."""
    ea, eb = EA, EB
    if MODE == "apart":
        ea = [A if x == Bi else x for x in EA]; eb = [A if x == Bi else x for x in EB]
    elif MODE == "same":
        ea = EA + [A]; eb = EB + [Bi]
    inp = (f"{n} {len(ea)} {len(TWO)} {K} {LSIT} {IT2} {seed}\n" + "\n".join(f"{a} {b}" for a, b in zip(ea, eb)) + "\n"
           + "\n".join(f"{i} {j}" for i, j in TWO) + "\n" + "\n".join(map(str, init)) + "\n")
    out = subprocess.run([TABU2], input=inp, capture_output=True, text=True).stdout.split("\n")
    ok = out[0].startswith("OK"); f = int(out[0].split()[1]); col = [int(x) for x in out[1:1 + n]]
    if MODE == "apart": col[Bi] = col[A]
    return ok, f, col
def proper(col):
    if any(col[a] == col[b] for a, b in zip(EA, EB)): return False
    if MODE == "apart" and col[A] != col[Bi]: return False
    if MODE == "same" and col[A] == col[Bi]: return False
    return True
print(f"{IN}: n={len(V)} m={len(EA)}, {len(U)} directions, pool {len(pool)}; mode {MODE}, solver {SOLVER}"
      + (f", pair {A},{Bi} at d^2 = {V[A].dist2(V[Bi])}" if A is not None else "") + f"   [{time.time()-t0:.0f}s]", flush=True)
last = 0; lsfail = 0
if d.get("colouring") and len(d["colouring"]) == len(V): col = list(d["colouring"])   # resume from the saved colouring
for it in range(1, 10 ** 7):
    n = len(V)
    init = col + [-1] * (n - len(col))
    ok, info, lcol = tabucol(n, init, it)
    if ok and proper(lcol):
        r = True; col = lcol; conf = 0; how = f"LSw alike {info}"
    else:
        lsfail += 1
        try: s.set_phases([X(v, c) if lcol[v] == c else -X(v, c) for v in range(n) for c in range(K)])
        except Exception: pass
        save("pre_cdcl")                      # the exact instance CDCL is about to decide
        if os.environ.get("FALLBACK") == "kissat":
            # external kissat on the exact clause set (hard edges + MODE constraint + pinned triangle)
            cnf = OUT + ".cnf"; KB = os.environ.get("KISSAT", "kissat")
            cls = [[X(v, c) for c in range(K)] for v in range(n)] + [[-X(a, c), -X(b, c)] for a, b in zip(EA, EB) for c in range(K)]
            if MODE == "apart": cls += [[-X(A, c), X(Bi, c)] for c in range(K)] + [[X(A, c), -X(Bi, c)] for c in range(K)]
            elif MODE == "same": cls += [[-X(A, c), -X(Bi, c)] for c in range(K)]
            if tri: cls += [[X(v, k)] for k, v in enumerate(tri)]
            with open(cnf, "w") as fh:
                fh.write(f"p cnf {n * K} {len(cls)}\n"); fh.write("".join(" ".join(map(str, c)) + " 0\n" for c in cls))
            kout = subprocess.run([KB, f"--time={int(os.environ.get('KTIME', '7200'))}", cnf], capture_output=True, text=True).stdout
            status = next((l for l in kout.split("\n") if l.startswith("s ")), "s UNKNOWN")
            conf = -1
            if "UNSATISFIABLE" in status: r = False
            elif "SATISFIABLE" in status:
                pos = {int(x) for l in kout.split("\n") if l.startswith("v") for x in l.split()[1:] if int(x) > 0}
                col = [next((c for c in range(K) if X(v, c) in pos), 0) for v in range(n)]
                r = True if proper(col) else None
            else: r = None
            how = f"KISSAT[{status}]"
        else:
            s.conf_budget(BUDGET)
            r = s.solve_limited(); st = s.accum_stats(); conf = st.get("conflicts", 0) - last; last = st.get("conflicts", 0)
        if not os.environ.get("FALLBACK") == "kissat": how = f"CDCL[LS best {info}]"
        if r is True and not os.environ.get("FALLBACK") == "kissat":
            mdl = s.get_model(); col = [next(c for c in range(K) if mdl[X(v, c) - 1] > 0) for v in range(n)]
    if r is not True:
        save(("UNSAT_" + MODE) if r is False else "hard")
        print(f"  iter {it}: n={n} m={len(EA)}: {'UNSAT (' + MODE + ') -- a claim until verified' if r is False else 'beyond budget'}   [{time.time()-t0:.0f}s]", flush=True)
        break
    if TWO:
        alike = [(i, j) for i, j in TWO if col[i] == col[j]]
        _tgt = _np.array([[Vfx[v], Vfy[v]] for p_ in alike for v in p_]) if alike else None
        if it % 10 == 1: print(f"    {len(alike)} of {len(TWO)} skeleton two-edges alike in this colouring", flush=True)
    rb = sorted(((0.25 * len(e[2]) - WN * near(kk), kk) for kk, e in pool.items()
                 if len(e[2]) >= 5 and len({col[j] for j in e[2]}) == K), reverse=True)
    if it % 10 == 1 or not rb or how.startswith(("CDCL", "KISSAT")):
        print(f"  iter {it}: n={n} m={len(EA)} {how} conflicts {conf}; {len(rb)} rainbow; LS fails so far {lsfail}   [{time.time()-t0:.0f}s]", flush=True)
    if it % 10 == 0: save("checkpoint")
    pick = [kk for _, kk in rb[:R]]
    if len(pick) < R:
        # no (or few) rainbow points: take conflicting TIGHT pairs -- two adjacent pool points that
        # each see four colours and miss the same one, so the colouring cannot extend to both
        tight = {}
        for kk, e in pool.items():
            if len(e[2]) >= 4:
                cs = {col[j] for j in e[2]}
                if len(cs) == K - 1: tight[kk] = (K * (K - 1) // 2 - sum(cs), len(e[2]))
        pairs = []
        for kk, (fc, dg) in tight.items():
            for k2 in range(len(U)):
                k3 = fkey(kk[0] + Ufx[k2], kk[1] + Ufy[k2])
                t3 = tight.get(k3)
                if t3 and t3[0] == fc and kk < k3: pairs.append((0.25 * (dg + t3[1]) - WN * (near(kk) + near(k3)) / 2, kk, k3))
        pairs.sort(reverse=True); have = set(pick)
        for _, k1, k2 in pairs:
            if len(pick) >= R: break
            for kx in (k1, k2):
                if kx not in have: pick.append(kx); have.add(kx)
        if it % 10 == 1 or not rb: print(f"    {len(tight)} tight pool points, {len(pairs)} conflicting tight pairs", flush=True)
    if not pick:
        # densify: the richest pool points, nearest first among equals
        pick = [kk for _, kk in sorted(((len(e[2]) - WN * near(kk), kk) for kk, e in pool.items()), reverse=True)[:R]]
        if it % 10 == 1: print(f"    no rainbow and no conflicting tight pair: densifying with the {len(pick)} richest pool points", flush=True)
    if not pick:
        save("pool_exhausted"); print("  pool exhausted", flush=True); break
    for kk in pick:
        i, k, nb = pool.pop(kk)
        p = V[i] + U[k]; v = len(V)
        vkey[key(p)] = v; vkey.setdefault(kk, v)
        V.append(p); Vfx.append(p.fx); Vfy.append(p.fy); adj.append(set())
        s.add_clause([X(v, c) for c in range(K)])
        for k2 in range(len(U)):
            j = vkey.get(fkey(p.fx + Ufx[k2], p.fy + Ufy[k2]))
            if j is not None and j != v and j not in adj[v]:
                adj[v].add(j); adj[j].add(v); EA.append(v); EB.append(j)
                for c in range(K): s.add_clause([-X(v, c), -X(j, c)])
    for i in range(n, len(V)): feed(i)
