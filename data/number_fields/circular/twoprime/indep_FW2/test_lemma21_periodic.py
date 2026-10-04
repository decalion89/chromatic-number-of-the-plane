"""Test the mechanism of Lemmas 21 and 22 on random *periodic* homomorphisms of Cay(Gamma, S) to K_{p/q},
p/q = chi_c = 1/kappa(S) in (2,4), Gamma = Z^r x Z/m.

A homomorphism of the finite quotient Cay(Gamma/Lambda, S) (Lambda = N1 Z x ... a full-rank sublattice of the free
part) lifts to a periodic homomorphism c of Cay(Gamma, S); for it the invariant mean is the plain average over the
quotient.  For random such c (SAT with random assumptions) we check:
 (1) a(s) = mean_x l(x,s) lies in [q, p-q] and a(-s) = p - a(s);
 (2) s -> a(s)/p mod 1 extends to a character of Gamma (Lemma 21);
 (3) T = {s : a(s) = q} carries a positive integer relation (the key step of Lemma 22);
 (4) the closed walk with n_t steps t (from the relation), started at every vertex of the period, is tight at every
     step (Lemma 21), so c has a tight cycle.
"""
import random, sys, json, itertools, math
from fractions import Fraction as Fr
from pysat.solvers import Solver
from kappa_exact import Group, symmetrize, half, kappa, positive_relation, check_relation
from sat_tight import Enc, encode_hom, decode, verify_hom


def torus(G, S, Ns):
    r, m = G.r, G.m
    V = [tuple(x) + (y,) for x in itertools.product(*[range(n) for n in Ns]) for y in range(m)]
    idx = {v: i for i, v in enumerate(V)}

    def red(v):
        return tuple(a % n for a, n in zip(v[:r], Ns)) + (v[r] % m,)

    E = set()
    for v in V:
        for s in S:
            w = red(tuple(a + b for a, b in zip(v[:r], s[:r])) + (v[r] + s[r],))
            assert w != v, "quotient too small: a generator vanishes"
            a, b = idx[v], idx[w]
            E.add((min(a, b), max(a, b)))
    return V, idx, red, sorted(E)


def is_character(G, S, vals):
    """vals: s -> Fraction in [0,1). Is there a character xi with xi(s) = vals[s] mod 1 for all s?"""
    r, m = G.r, G.m
    free = [s for s in S if any(s[:r])]
    # choose r linearly independent free parts
    basis = None
    for comb in itertools.combinations(free, r):
        M = [[Fr(v) for v in s[:r]] for s in comb]
        det = M[0][0] if r == 1 else M[0][0] * M[1][1] - M[0][1] * M[1][0]
        if det != 0:
            basis = comb
            break
    assert basis is not None
    M = [[s[i] for i in range(r)] for s in basis]
    det = M[0][0] if r == 1 else M[0][0] * M[1][1] - M[0][1] * M[1][0]
    D = abs(det)
    for b in range(m):
        for js in itertools.product(range(D), repeat=r):
            # solve alpha . x_s = vals[s] - b y_s/m + j_s for s in basis
            rhs = [vals[s] - Fr(b * s[r], m) + j for s, j in zip(basis, js)]
            if r == 1:
                alpha = (rhs[0] / M[0][0],)
            else:
                a11, a12, a21, a22 = M[0][0], M[0][1], M[1][0], M[1][1]
                alpha = ((rhs[0] * a22 - a12 * rhs[1]) / det, (a11 * rhs[1] - a21 * rhs[0]) / det)
            ok = True
            for s in S:
                v = sum(a * x for a, x in zip(alpha, s[:r])) + Fr(b * s[r], m)
                if (v - vals[s]).denominator != 1:
                    ok = False
                    break
            if ok:
                return alpha, b
    return None


def run_example(name, r, m, Sh, Ns, n_homs=30, seed=0):
    rng = random.Random(seed)
    G = Group(r, m)
    S = symmetrize(G, Sh)
    kap, opt = kappa(G, S)
    chi = 1 / kap
    p, q = chi.numerator, chi.denominator
    V, idx, red, E = torus(G, S, Ns)
    N = len(V)
    enc = Enc()
    x = encode_hom(N, E, p, q, enc)
    out = dict(name=name, chi_c=str(chi), torus=Ns, V=N, homs=0, non_character_homs=0, ok=0, fail=[])
    with Solver(name='cadical153', bootstrap_with=enc.cls) as sol:
        if not sol.solve():
            out['note'] = 'quotient has no hom to K_{p/q}'
            return out
        attempts = 0
        while out['homs'] < n_homs and attempts < 20 * n_homs:
            attempts += 1
            k = rng.randint(1, max(1, N // 6))
            assum = []
            for v in rng.sample(range(N), k):
                assum.append(x[v][rng.randrange(p)])
            if not sol.solve(assumptions=assum):
                continue
            model = set(l for l in sol.get_model() if l > 0)
            col = decode(model, x, N, p)
            assert verify_hom(col, E, p, q)
            out['homs'] += 1

            def c(v):
                return col[idx[red(v)]]

            def ell(v, s):
                w = tuple(a + b for a, b in zip(v[:r], s[:r])) + (v[r] + s[r],)
                d = (c(w) - c(v)) % p
                assert q <= d <= p - q
                return d
            a = {}
            const = True
            for s in S:
                vals = [ell(v, s) for v in V]
                if len(set(vals)) > 1:
                    const = False
                a[s] = Fr(sum(vals), N)
            if not const:
                out['non_character_homs'] += 1
            ok1 = all(q <= a[s] <= p - q and a[G.neg(s)] == p - a[s] for s in S)
            vals = {s: (a[s] / p) % 1 for s in S}
            ok2 = is_character(G, S, vals) is not None
            T = [s for s in S if a[s] == q]
            rel = positive_relation(G, T)
            ok3 = rel is not None and check_relation(G, rel)
            ok4 = False
            if ok3:
                steps = [t for t, n in rel.items() for _ in range(n)]
                rng.shuffle(steps)
                ok4 = True
                for v0 in V:
                    v = v0
                    for t in steps:
                        if ell(v, t) != q:
                            ok4 = False
                            break
                        v = tuple(a_ + b_ for a_, b_ in zip(v[:r], t[:r])) + (v[r] + t[r],)
                    if not ok4:
                        break
            if ok1 and ok2 and ok3 and ok4:
                out['ok'] += 1
            else:
                out['fail'].append(dict(ok=[ok1, ok2, ok3, ok4], a={str(s): str(v) for s, v in a.items()}))
    return out


EX = [
    ("Z, D={3,4,9,12} (7/2)", 1, 1, [(3, 0), (4, 0), (9, 0), (12, 0)], [14]),
    ("Z, D={3,4,9,12} (7/2)", 1, 1, [(3, 0), (4, 0), (9, 0), (12, 0)], [21]),
    ("Z, D={3,4,9,12} (7/2)", 1, 1, [(3, 0), (4, 0), (9, 0), (12, 0)], [28]),
    ("Z, D={1,3,8} (11/4)", 1, 1, [(1, 0), (3, 0), (8, 0)], [22]),
    ("Z, D={1,3,8} (11/4)", 1, 1, [(1, 0), (3, 0), (8, 0)], [33]),
    ("Z x Z/2, S={(5,1),(3,0),(2,1)} (16/5)", 1, 2, [(5, 1), (3, 0), (2, 1)], [16]),
    ("Z x Z/2, S={(5,1),(3,0),(2,1)} (16/5)", 1, 2, [(5, 1), (3, 0), (2, 1)], [32]),
    ("Z^2, S={(-1,2),(-1,1),(0,1),(3,0)} (10/3)", 2, 1, [(-1, 2, 0), (-1, 1, 0), (0, 1, 0), (3, 0, 0)], [10, 10]),
    ("Z^2, S={(-1,2),(-1,1),(0,1),(3,0)} (10/3)", 2, 1, [(-1, 2, 0), (-1, 1, 0), (0, 1, 0), (3, 0, 0)], [20, 10]),
    ("Z^2, S={(1,0),(0,1),(1,2)}", 2, 1, [(1, 0, 0), (0, 1, 0), (1, 2, 0)], [10, 10]),
    ("Z x Z/4, S={(4,1),(2,2),(1,2)} (12/5)", 1, 4, [(4, 1), (2, 2), (1, 2)], [24]),
]

if __name__ == '__main__':
    res = []
    for (name, r, m, Sh, Ns) in EX:
        o = run_example(name, r, m, Sh, Ns, n_homs=int(sys.argv[1]) if len(sys.argv) > 1 else 25)
        print(json.dumps({k: v for k, v in o.items() if k != 'fail'}), ('FAILS: ' + str(o['fail'][:2])) if o['fail'] else '', flush=True)
        res.append(o)
    json.dump(res, open('lemma21_periodic_results.json', 'w'), indent=1)
