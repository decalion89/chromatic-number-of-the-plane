"""f(k): the size of the smallest unit-distance configuration carrying a forced
monochromatic pair under k colours.

f(3) = 4: the unit rhombus forces its two far tips, and nothing smaller can --
three points cannot force anything with three colours available.

f(4) has never been tabulated as far as this project has found. Measure it by
shrinking, from a ball known to carry forcing, using the solver's UNSAT core
over vertex selectors and then greedy deletion.

The ratio is what matters. It says what scale a k=5 forcing configuration
would sit at, and therefore whether a search that reaches ten thousand
vertices is close or is missing orders of magnitude.
"""
import sys, time
from fractions import Fraction as Fr
sys.path.insert(0,'/home/user/darwin-50/research/hadwiger-nelson')
from pysat.formula import CNF
from pysat.solvers import Solver
from hn.degrey import build_G
from hn.spindle import local_ball

G = build_G()
centre = min(range(G.n), key=lambda v: G.vertices[v].fx**2 + G.vertices[v].fy**2)

def forced_pair_cnf(g, k, pivot, targets, active):
    """CNF whose UNSAT under {a_v : v in active} means the (pivot,targets)
    disjunction is forced in the induced subgraph on `active`."""
    n = g.n
    x = lambda v,c: 1 + v*k + c
    a = lambda v: 1 + n*k + v
    sel = 2 + n*k + n
    cnf = CNF()
    for v in range(n):
        cnf.append([-a(v)] + [x(v,c) for c in range(k)])
    for u,v in g.edges():
        for c in range(k):
            cnf.append([-a(u), -a(v), -x(u,c), -x(v,c)])
    for q in targets:
        for c in range(k):
            cnf.append([-sel, -x(pivot,c), -x(q,c)])
    return cnf, a, sel

def minimise(g, k, pivot, targets, verbose=True):
    cnf, a, sel = forced_pair_cnf(g, k, pivot, targets, None)
    s = Solver(name='cd19', bootstrap_with=cnf)
    keep = set(range(g.n))
    try:
        assume = [sel] + [a(v) for v in keep]
        if s.solve(assumptions=assume) is not False:
            return None
        core = set(l - 1 - g.n*k for l in (s.get_core() or []) if l != sel and l > g.n*k)
        core |= {pivot} | set(targets)
        if verbose: print(f"    first core: {len(core)}", flush=True)
        for _ in range(30):
            res = s.solve(assumptions=[sel] + [a(v) for v in core])
            if res is not False: break
            new = set(l - 1 - g.n*k for l in (s.get_core() or []) if l != sel and l > g.n*k)
            new |= {pivot} | set(targets)
            if not new or len(new) >= len(core): break
            core = new
            if verbose: print(f"    shrunk to {len(core)}", flush=True)
        order = sorted(core, key=lambda v: len(g.adj[v]))
        for v in order:
            if v in (pivot,) or v in targets or v not in core: continue
            trial = core - {v}
            if s.solve(assumptions=[sel] + [a(t) for t in trial]) is False:
                core = trial
        return core
    finally:
        s.delete()

print("Measuring f(4): smallest configuration with a forced pair at k=4", flush=True)
for r in (2.5, 3.0, 3.5):
    ball, bp = local_ball(G, centre, r)
    p = ball.vertices[bp]
    groups = {}
    for j in range(ball.n):
        if j == bp: continue
        d2 = p.dist2(ball.vertices[j])
        if d2.is_rational() and Fr(1,4) <= d2.c[0] <= Fr(40):
            groups.setdefault(d2.c[0], []).append(j)
    for val in sorted(groups):
        t = time.time()
        core = minimise(ball, 4, bp, groups[val], verbose=False)
        if core:
            print(f"  radius {r} (n={ball.n}), d2={val}: FORCED, minimal core "
                  f"{len(core)} vertices  [{time.time()-t:.0f}s]", flush=True)
            sub = ball.induced(sorted(core))
            print(f"     induced subgraph: {sub}", flush=True)
            break
    else:
        continue
    break
