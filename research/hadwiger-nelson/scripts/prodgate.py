"""The product gate: periodic colourings through NON-cyclic quotients.

A twisted colouring psi + t*floor(phi/L) is periodic with quotient Z/5L x Z/5 --
not cyclic, so the Z/n gates never saw it.  Here chi = (phi_1 mod n_1, ...,
phi_k mod n_k) : M -> A = Z/n_1 x ... x Z/n_k is sampled (optionally with one
coordinate an admissible psi mod 5), and ANY proper 5-colouring G of
Cay(A, chi(U)) is searched for, not just the staircase ones:
  apart : G(0) = G(chi(2e)) for some direction e  -> refutes a 2e gadget along e
  pair  : G(0) != G(chi(5e)) for some direction e -> refutes forcing 5e along e
(vertex-transitivity lets the witness sit at 0).  Refuted directions are
dropped and the same quotient is asked again for the rest.

usage: prodgate.py <module.json> <apart|pair> <samples> <n_1,n_2,...> [adm]
"""
import sys, json, random, itertools, time
from fractions import Fraction as Fr
exec(open(__import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "gate.py")).read().split("def gate(g, label):")[0])
from pysat.solvers import Solver
ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
name, MODE, NS = sys.argv[1], sys.argv[2], int(sys.argv[3])
MODS = [int(x) for x in sys.argv[4].split(",")]; ADM = len(sys.argv) > 5 and sys.argv[5] == "adm"
d = json.load(open(f"{ROOT}/data/{name}"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
E = edge_vectors(build_graph(P))
B = echelon(E); C = [coords(B, v) for v in E]; r = len(B)
adm = [psi for psi in itertools.product(range(5), repeat=r) if all(sum(a * b for a, b in zip(psi, c)) % 5 for c in C)] if ADM else None
print(f"{name}: {len(C)} directions, rank {r}; quotient Z/" + " x Z/".join(map(str, MODS))
      + (f"; first factor an admissible psi ({len(adm)})" if ADM else "") + f"; mode {MODE}", flush=True)
elts = list(itertools.product(*[range(n) for n in MODS])); idx = {v: i for i, v in enumerate(elts)}; N = len(elts)
add = lambda a, b: tuple((x + y) % n for x, y, n in zip(a, b, MODS))
X = lambda v, c: 1 + v * 5 + c
rng = random.Random(7); t0 = time.time()
alive = set(range(len(C))); tried = loopless = colourable = 0
for smp in range(NS):
    phis = [[rng.randrange(n) for _ in range(r)] for n in MODS]
    if ADM: phis[0] = list(rng.choice(adm))           # needs MODS[0] == 5
    img = [tuple(sum(a * b for a, b in zip(ph, c)) % n for ph, n in zip(phis, MODS)) for c in C]
    tried += 1
    zero = tuple(0 for _ in MODS)
    if zero in img: continue
    loopless += 1
    S = set(img) | {tuple((-x) % n for x, n in zip(v, MODS)) for v in img}
    cl = [[X(v, c) for c in range(5)] for v in range(N)]
    for i, v in enumerate(elts):
        for s_ in S:
            j = idx[add(v, s_)]
            if i < j:
                for c in range(5): cl.append([-X(i, c), -X(j, c)])
    cl.append([X(0, 0)])
    sol = Solver(name="cd19", bootstrap_with=cl)
    if not sol.solve(): sol.delete(); continue
    colourable += 1
    top = N * 5 + 1
    while alive and MODE != "plain":
        k2 = 2 if MODE == "apart" else 5
        targets = {k: idx[tuple((k2 * x) % n for x, n in zip(img[k], MODS))] for k in alive}
        sel = {}
        for k, w in targets.items():
            z = top; top += 1; sel[k] = z
            if MODE == "apart":     # z -> G(0) = G(w): G(0) = 0 is pinned, so z -> X(w, 0)
                sol.add_clause([-z, X(w, 0)])
            else:                   # z -> G(w) != 0
                sol.add_clause([-z, -X(w, 0)])
        act = top; top += 1
        sol.add_clause([-act] + list(sel.values()))
        if not sol.solve(assumptions=[act]): break
        m = sol.get_model()
        col = [next(c for c in range(5) if m[X(v, c) - 1] > 0) for v in range(N)]
        hit = {k for k, w in targets.items() if (col[w] == col[0]) == (MODE == "apart")}
        alive -= hit
        print(f"  sample {smp}: chi images {sorted(set(img))[:6]}..., refutes {len(hit)} directions; {len(alive)} left   [{time.time()-t0:.0f}s]", flush=True)
    sol.delete()
    if not alive: break
print(f"  {tried} sampled, {loopless} loopless, {colourable} 5-colourable; directions unrefuted: {len(alive)} of {len(C)}   [{time.time()-t0:.0f}s]", flush=True)
