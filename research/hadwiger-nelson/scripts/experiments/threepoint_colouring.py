"""A variant of scripts/threepoint.py with --colour k: the three-point programme for a proper k-colouring,
kept as the record of a negative result (notes/local_colourings.md, section 14; research log, 25 September).

For a k-colouring c, recentre at a random point t, turn by a random isometry and relabel the colours so that
the origin gets colour 1 (the other k - 1 labels in random order).  g(w) = P(c(w) = c(0)), z(a, b) =
P(c(0) = c(a) = c(b)), and sum_w g(w) = sum_C |C|^2 / n >= n / k.  Besides M1 and M0 (whose empty-set row now
has the exact corner k - 1), the moment matrix of the colour indicators has a component on the standard
representation of the symmetric group of the other k - 1 colours:
    M2 = [(k - 2) p_ABB(a, b) - p_ABC(a, b)] >= 0,   p_ABB = g(b - a) - z,   p_ABC = 1 - g(a) - g(b) - g(b - a) + 2 z,
and p_ABC >= 0 on every triangle.  A maximum below n / k proves chi >= k + 1.

Result: the new constraints change the bound by about 0.01 at most (G_17: 64.556 -> 64.545; F_23^2:
108.547 -> 108.542, and 107.04 with the localizing matrices either way).  A real 5-colouring of F_13^2 satisfies every
constraint (checked).  Not maintained, and scripts/threepoint_verify.py does not read its --colour output.

The original docstring follows.

Three-point bound for the independence number of a finite unit-distance plane
(notes/local_colourings.md, section 14).

The plane is F_q^2 with Q(x, y) = x^2 + y^2 ('std', Moorhouse's table) or x^2 - n y^2 with n a non-square
('inert', the anisotropic plane G_q of section 4); two points are adjacent when Q of their difference is 1.
The graph is vertex-transitive, so chi >= q^2 / alpha, and alpha < q^2 / 5 gives chi >= 6.

**Variables.**  Recentre an independent set S at a random point of S and turn it by a random isometry
fixing that point.  g_c is the probability that a given point w with Q(w) = c lies in the moved set (c runs
over the circles and, when Q is isotropic, the nonzero isotropic vectors), and z(a, b) the probability that
a and b both do; z depends only on the isometry class of the triangle {0, a, b} (with its vertices
permuted), and vanishes when a side is a unit vector.  |S| = sum_w g(w), with g(0) = 1.

**Constraints.**
- M1 = [z(a, b)]_{a, b} (z(a, a) = g(a), z(0, b) = g(b), z(0, 0) = 1) and M0 = [g(b - a) - z(a, b)]_{a, b != 0}
  are positive semidefinite: they average x_0 v v^T and (1 - x_0) v v^T, v the indicator of the moved set.
- 0 <= z <= g(side) for each side, and 1 - g(a) - g(b) + z(a, b) >= 0 (with the vertices permuted).
- The Delsarte inequalities sum_c g_c P_c(xi) + 1 >= 0 for every character xi.
Optional (--tri, --trilocal, --edge, --pentagon): conditional triangle inequalities, and the localizing
matrices of 1 - x_0 - x_e - x_f for a unit triangle {0, e, f}, of 1 - x_0 - x_e for a unit vector e and of
2 - sum_t x_t for a unit 5-cycle.  With a target K0 (--test K0) these matrices and that of 1 - x_0 get their
empty-set rows, whose corner k/delta - |T| (delta = |S|/n) is replaced by k n/K0 - |T|: valid for every S
with |S| >= K0.  A bound below K0 then proves alpha < K0.

**Symmetry.**  Both matrices commute with the rotation group SO(Q), cyclic of order N = q -/+ 1.  In the
basis of the omega^k-eigenvectors of a rotation, made real by a reflection (J = reflection o conjugation),
they split into N/2 + 1 real blocks of size about q.  The localizing matrices split by the characters of
the Klein group of the edge and of the symmetric group of the triangle.

**Solvers.**  The default ('dsdp') solves the primal with DSDP (through cvxopt), which is accurate in the
primal but not in the dual, and then finds the dual on the null spaces of the blocks at the optimum: each
psd multiplier is U W U^T, with U a basis of the null space of its block and W >= 0 small, and with the
multipliers of the active linear constraints this is a small semidefinite programme (Clarabel).  dual_lp,
used for three of the stored certificates, restricts W to nonnegative combinations of rank-one directions
and solves a linear programme instead.  'clarabel' and 'cvxopt' solve the whole programme at once (Clarabel
needs about 2 GB per thousand triangle variables; cvxopt is slow beyond q = 37).
scripts/threepoint_verify.py checks the saved dual rigorously.

usage: python3 scripts/threepoint.py q std|inert [--solver dsdp|clarabel|cvxopt] [--tri] [--trilocal]
                                                  [--edge] [--test K0] [--save certificate.npz]
"""
import sys, time, json, argparse, itertools
import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
from scipy.optimize import linprog


# ------------------------------------------------------------------------------------------ the plane
def nonsq(q):
    sq = {t * t % q for t in range(1, q)}
    return next(n for n in range(2, q) if n not in sq)


def form(q, kind):
    """n0 with Q(x, y) = x^2 + n0 y^2 (mod q)"""
    return 1 if kind == 'std' else (-nonsq(q)) % q


def orthogonal_group(q, n0):
    """all 2 x 2 matrices over F_q preserving Q(x, y) = x^2 + n0 y^2"""
    Qf = lambda x, y: (x * x + n0 * y * y) % q
    V1 = [(x, y) for x in range(q) for y in range(q) if Qf(x, y) == Qf(1, 0)]
    V2 = [(x, y) for x in range(q) for y in range(q) if Qf(x, y) == Qf(0, 1)]
    return [np.array([[a, b], [c, d]]) for (a, c) in V1 for (b, d) in V2 if (a * b + n0 * c * d) % q == 0]


def setup(q, kind):
    """(n0, Q, a generator sigma of SO(Q), a reflection tau, N = |SO(Q)|)"""
    n0 = form(q, kind)
    Qf = lambda x, y: (x * x + n0 * y * y) % q
    mats = orthogonal_group(q, n0)
    det = lambda R: (R[0, 0] * R[1, 1] - R[0, 1] * R[1, 0]) % q
    rots = [R for R in mats if det(R) == 1]
    refl = [R for R in mats if det(R) != 1]
    N = len(rots)

    def order(R):
        M = R.copy(); k = 1
        while not (M[0, 0] == 1 and M[1, 1] == 1 and M[0, 1] == 0 and M[1, 0] == 0):
            M = (M @ R) % q; k += 1
        return k
    sig = next(R for R in rots if order(R) == N)
    return n0, Qf, sig, refl[0], N


def character_table(q, kind):
    """(class sizes, distinct rows P[xi, c] = sum_{w in class c} cos(2 pi xi.w / q), class array);
    class 0 = {0}, class r = {Q(w) = r} for r != 0, class q = nonzero isotropic vectors"""
    n0 = form(q, kind)
    X, Y = np.meshgrid(np.arange(q), np.arange(q), indexing='ij')
    Qv = (X * X + n0 * Y * Y) % q
    cls = np.where(Qv == 0, q, Qv); cls[0, 0] = 0
    size = np.bincount(cls.ravel(), minlength=q + 1)
    W = np.stack([X.ravel(), Y.ravel()], 1); C = cls.ravel()
    rows = {}
    for xi in range(1, q * q):
        a, b = divmod(xi, q)
        ph = np.cos(2 * np.pi * ((a * W[:, 0] + b * W[:, 1]) % q) / q)
        r = np.bincount(C, weights=ph, minlength=q + 1)
        rows.setdefault(tuple(np.round(r, 6)), r)
    return size, np.array(list(rows.values())), cls


# ------------------------------------------------------------------------------------------ the programme
def build(q, kind, verbose=True, K0=None, tri=False, trilocal=False, edge=False, pentagon=False, colour=None):
    """blocks (size, {(i, j): ({var: coef}, const)}), linear rows ({var: coef}, rhs) meaning coef.v <= rhs.
    colour = k: the programme for a k-colouring (see the module docstring); K0 must then be None."""
    t0 = time.time()
    n = q * q
    assert colour is None or (K0 is None and colour >= 3)
    delta0 = None if K0 is None else K0 / n
    if colour is not None:
        delta0 = 1.0 / colour                         # corners k * colour - |T|: exact for a colouring
    n0, Qf, sig, tau, N = setup(q, kind)
    X, Y = np.divmod(np.arange(n), q)
    idx = lambda x, y: (x % q) * q + (y % q)
    Qv = Qf(X, Y)
    cls = np.where(Qv == 0, q, Qv); cls[0] = 0
    apply = lambda R, x, y: ((R[0, 0] * x + R[0, 1] * y) % q, (R[1, 0] * x + R[1, 1] * y) % q)
    s_img = idx(*apply(sig, X, Y))
    t_img = idx(*apply(tau, X, Y))
    neg = idx(-X, -Y)
    powimg = np.zeros((N, n), int); powimg[0] = np.arange(n)
    for j in range(1, N):
        powimg[j] = s_img[powimg[j - 1]]
    # rotation orbits; tau-stable ones get the half-shift t_a/2, pairs swapped by tau are re-based (tau p_a = p_a')
    orb = -np.ones(n, int); reps = []
    for v in range(1, n):
        if orb[v] < 0:
            orb[powimg[:, v]] = len(reps); reps.append(v)
    m = len(reps)
    partner = [-1] * m; shift = np.zeros(m); done = [False] * m
    for a in range(m):
        if done[a]:
            continue
        tp = t_img[reps[a]]; b = orb[tp]
        if b == a:
            shift[a] = int(np.nonzero(powimg[:, reps[a]] == tp)[0][0]) / 2; done[a] = True
        else:
            reps[b] = int(tp); partner[a] = b; partner[b] = a; done[a] = done[b] = True
    # triangle classes: orbits of ordered pairs (a, b) under rotation, reflection, swap and recentring at a
    P = np.arange(n * n); A_, B_ = np.divmod(P, n)
    sub = lambda u, v: idx(X[u] - X[v], Y[u] - Y[v])
    gens = [s_img[A_] * n + s_img[B_], t_img[A_] * n + t_img[B_], B_ * n + A_, neg[A_] * n + sub(B_, A_)]
    G = sp.coo_matrix((np.ones(4 * n * n, np.int8), (np.tile(P, 4), np.concatenate(gens))), shape=(n * n, n * n))
    ncomp, lab = connected_components(G, directed=True, connection='weak')
    del G, gens
    unit = (cls[A_] == 1) | (cls[B_] == 1) | (cls[sub(B_, A_)] == 1)
    degen = (A_ == 0) | (B_ == 0) | (A_ == B_)
    tri_ok = ~unit & ~degen
    good = np.unique(lab[tri_ok])
    assert len(np.intersect1d(good, np.unique(lab[unit & ~degen]))) == 0
    zid = -np.ones(ncomp, int); zid[good] = np.arange(len(good)); nz = len(good)
    free = [c for c in range(2, q + 1) if np.any(cls == c)]
    gcol = {c: nz + i for i, c in enumerate(free)}
    nv = nz + len(free)
    if verbose:
        print(f"{kind}{q}: n={n}, N={N}, {m} rotation orbits ({sum(p >= 0 for p in partner)} paired), "
              f"{nz} triangle variables, {len(free)} circle variables [{time.time()-t0:.0f}s]", flush=True)

    def Mrow(x, ys, which):
        """M1(x, y) (which = 1), M0(x, y) (which = 0) or, for a colouring with k = colour colours, the standard
        block M2(x, y) = (k-1) g(y-x) + g(x) + g(y) - 1 - k z(x, y), M2(x, x) = (k-2)(1 - g(x)) (which = 2), for
        y in ys, as ([(var array, coef array), ...], const array)"""
        ys = np.asarray(ys); L = len(ys)
        slots = []; cst = np.zeros(L)
        new_slot = lambda: (slots.append((-np.ones(L, int), np.zeros(L))) or slots[-1])
        pair = x * n + ys
        if which == 1:
            vi, co = new_slot()
            deg = (ys == 0) | (ys == x) | (x == 0)
            cdeg = np.where(ys == 0, cls[x], np.where(x == 0, cls[ys], cls[x]))
            for i in np.nonzero(deg)[0]:
                cc = cdeg[i]
                if cc == 0: cst[i] = 1.0
                elif cc != 1: vi[i] = gcol[cc]; co[i] = 1.0
            nd = ~deg & ~unit[pair]
            vi[nd] = zid[lab[pair[nd]]]; co[nd] = 1.0
        elif which == 0:
            vi, co = new_slot(); vi2, co2 = new_slot()
            for i in np.nonzero(ys == x)[0]:
                cst[i] = 1.0
                if cls[x] != 1: vi[i] = gcol[cls[x]]; co[i] = -1.0
            cdiff = cls[sub(ys, np.full(L, x))]
            for i in np.nonzero(ys != x)[0]:
                if cdiff[i] != 1: vi[i] = gcol[cdiff[i]]; co[i] = 1.0
                if not unit[pair[i]]: vi2[i] = zid[lab[pair[i]]]; co2[i] = -1.0
        else:
            kc = colour
            vd, cd = new_slot(); vx, cx = new_slot(); vy, cy = new_slot(); vz, cz = new_slot()
            for i in np.nonzero(ys == x)[0]:
                cst[i] = kc - 2.0
                if cls[x] != 1: vx[i] = gcol[cls[x]]; cx[i] = -(kc - 2.0)
            cdiff = cls[sub(ys, np.full(L, x))]
            for i in np.nonzero(ys != x)[0]:
                cst[i] = -1.0
                if cdiff[i] != 1: vd[i] = gcol[cdiff[i]]; cd[i] = kc - 1.0
                if cls[x] != 1: vx[i] = gcol[cls[x]]; cx[i] = 1.0
                if cls[ys[i]] != 1: vy[i] = gcol[cls[ys[i]]]; cy[i] = 1.0
                if not unit[pair[i]]: vz[i] = zid[lab[pair[i]]]; cz[i] = -float(kc)
        return slots, cst

    omega = np.exp(2j * np.pi / N)
    blocks = []
    for which in ((1, 0) if colour is None else (1, 0, 2)):
        rowdata = {(a, b): Mrow(reps[a], powimg[:, reps[b]], which) for a in range(m) for b in range(m)}
        for k in range(N // 2 + 1):
            wl = omega ** (k * np.arange(N))
            H = {}
            for (a, b), (slots, cst) in rowdata.items():
                ph = omega ** (k * (shift[a] - shift[b]))
                acc = {}
                for arr_v, arr_c in slots:
                    msk = arr_v >= 0
                    for v_, c_ in zip(arr_v[msk], (arr_c * wl)[msk]):
                        acc[v_] = acc.get(v_, 0) + c_ * ph
                H[a, b] = (acc, ph * np.sum(cst * wl))
            cols = []; seen = set()
            for a in range(m):
                if a in seen:
                    continue
                b = partner[a]
                if b < 0:
                    cols.append({a: 1.0}); seen.add(a)
                else:
                    cols.append({a: 1 / np.sqrt(2), b: 1 / np.sqrt(2)})
                    cols.append({a: 1j / np.sqrt(2), b: -1j / np.sqrt(2)})
                    seen.update((a, b))
            with_origin = (which == 1 and k == 0)
            with_empty = (which == 0 and k == 0 and delta0 is not None)
            off = 1 if (with_origin or with_empty) else 0
            ent = {}
            for i, ci in enumerate(cols):
                for j in range(i, len(cols)):
                    cj = cols[j]; acc = {}; cst = 0j
                    for a, ta in ci.items():
                        for b, tb in cj.items():
                            w = np.conj(ta) * tb
                            hv, hc = H[a, b]
                            cst += w * hc
                            for v_, c_ in hv.items():
                                acc[v_] = acc.get(v_, 0) + w * c_
                    ent[i + off, j + off] = (acc, cst)
            if with_origin:
                ent[0, 0] = ({}, 1.0)
                for j, cj in enumerate(cols):
                    acc = {}; cst = 0j
                    for b, tb in cj.items():
                        slots0, c0 = Mrow(0, np.array([reps[b]]), 1)
                        vi, co = slots0[0]
                        if vi[0] >= 0: acc[vi[0]] = acc.get(vi[0], 0) + tb * np.sqrt(N) * co[0]
                        cst += tb * np.sqrt(N) * c0[0]
                    ent[0, j + 1] = (acc, cst)
            if with_empty:
                ent[0, 0] = ({}, 1.0 / delta0 - 1.0)
                for j, cj in enumerate(cols):
                    acc = {}; cst = 0j
                    for b, tb in cj.items():
                        c = cls[reps[b]]
                        cst += tb * np.sqrt(N)
                        if c != 1: acc[gcol[c]] = acc.get(gcol[c], 0) - tb * np.sqrt(N)
                    ent[0, j + 1] = (acc, cst)
            mx = max([abs(np.imag(c)) for (acc, cst) in ent.values() for c in list(acc.values()) + [cst]] + [0])
            assert mx < 1e-8, (which, k, mx)
            blocks.append((len(cols) + off, ent))
    nrot = len(blocks)
    if verbose:
        print(f"{nrot} rotation blocks of sizes {sorted(set(s for s, _ in blocks))} [{time.time()-t0:.0f}s]", flush=True)
    # linear constraints on triangle classes
    reps_t = {}
    for p in np.nonzero(tri_ok)[0]:
        reps_t.setdefault(lab[p], p)
    lin = []
    for L, p in reps_t.items():
        a, b = divmod(int(p), n)
        cs = [cls[a], cls[b], cls[sub(b, a)]]
        zi = zid[L]
        for c in cs:
            lin.append(({zi: 1.0, gcol[c]: -1.0}, 0.0))
        for i in range(3):
            for j in range(i + 1, 3):
                d = {zi: -1.0}
                d[gcol[cs[i]]] = d.get(gcol[cs[i]], 0) + 1.0
                d[gcol[cs[j]]] = d.get(gcol[cs[j]], 0) + 1.0
                lin.append((d, 1.0))
    if colour is not None:
        # three distinct colours on {0, a, b}: 1 - g(a) - g(b) - g(b - a) + 2 z(a, b) >= 0
        for L, p in reps_t.items():
            a, b = divmod(int(p), n)
            d = {zid[L]: -2.0}
            for c in (cls[a], cls[b], cls[sub(b, a)]):
                d[gcol[c]] = d.get(gcol[c], 0) + 1.0
            lin.append((d, 1.0))
        # points w, w' one unit step apart get different colours: g(w) + g(w') <= 1 (triangles {0, w, w'} with
        # a unit side)
        Uc = np.nonzero(cls == 1)[0]
        seenp = set()
        for w in range(1, n):
            if cls[w] == 1: continue
            for u in Uc:
                w2 = int(idx(X[w] + X[u], Y[w] + Y[u]))
                if w2 == 0 or cls[w2] == 1: continue
                key = tuple(sorted((int(cls[w]), int(cls[w2]))))
                if key in seenp: continue
                seenp.add(key)
                d = {}
                for c in key: d[gcol[c]] = d.get(gcol[c], 0) + 1.0
                lin.append((d, 1.0))
    size_c, Pm, _ = character_table(q, kind)
    U_ = np.nonzero(cls == 1)[0]; Us_ = set(U_.tolist())
    if tri:
        tri0 = [(u, v) for u in U_ for v in U_ if u < v and sub(v, u) in Us_]
        T = np.array([[x, idx(X[x] + X[u], Y[x] + Y[u]), idx(X[x] + X[v], Y[x] + Y[v])] for x in range(n) for (u, v) in tri0])
        T = np.unique(np.sort(T, 1), axis=0) if len(T) else T
        seen = set()
        for t3 in T:                                  # given 0 in S: sum_{t in T} g(t) <= 1
            if 0 in t3:
                continue
            d = {}
            for t in t3:
                if cls[t] != 1: d[gcol[cls[t]]] = d.get(gcol[cls[t]], 0) + 1.0
            key = tuple(sorted(d.items()))
            if d and key not in seen:
                seen.add(key); lin.append((d, 1.0))
        for c in free:                                # given 0, a in S: sum_{t in T} z(a, t) <= g(a)
            a = int(np.nonzero(cls == c)[0][0])
            for t3 in T:
                if 0 in t3 or a in t3:
                    continue
                d = {}
                for t in t3:
                    p = a * n + t
                    if not unit[p]: d[zid[lab[p]]] = d.get(zid[lab[p]], 0) + 1.0
                key = (c,) + tuple(sorted(d.items()))
                if d and key not in seen:
                    seen.add(key); d[gcol[c]] = d.get(gcol[c], 0) - 1.0; lin.append((d, 0.0))
    # localizing matrices: (1 - sum_{t in T} x_t) (1, x)(1, x)^T >= 0 for a clique T containing 0
    def localizing(T, group, chars, corner, k=1):
        """localizing matrix of k - sum_{t in T} x_t (valid when T has no independent set of k + 1 points);
        group: list of index maps preserving T; chars: list of lists of coefficient functions, one per
        irreducible character, giving for each group element its matrix coefficient (several functions for a
        2-dimensional character, one per copy)"""
        pts = [p for p in range(n) if p not in T]
        subs = [np.arange(n) if t == 0 else idx(X - int(X[t]), Y - int(Y[t])) for t in T]

        def L(xa, xb):
            d = {}; cst = 0.0
            if xa == xb:
                cst = float(k)
                for es in subs:
                    c = cls[es[xa]]
                    if c != 1: d[gcol[c]] = d.get(gcol[c], 0) - 1.0
                return d, cst
            c = cls[sub(xb, xa)]
            if c != 1: d[gcol[c]] = d.get(gcol[c], 0) + float(k)
            for es in subs:
                p_ = int(es[xa]) * n + int(es[xb])
                if not unit[p_]: d[int(zid[lab[p_]])] = d.get(int(zid[lab[p_]]), 0) - 1.0
            return d, cst
        seen_p = set(); orbs = []
        for p in pts:
            if p not in seen_p:
                seen_p.update(int(h[p]) for h in group); orbs.append(p)
        out = []; allvecs = []
        for ci, fns in enumerate(chars):
            vecs = []
            for p in orbs:
                for fn in fns:
                    coefs = {}
                    for gi, h in enumerate(group):
                        coefs[int(h[p])] = coefs.get(int(h[p]), 0) + fn(gi)
                    coefs = {k_: v_ for k_, v_ in coefs.items() if abs(v_) > 1e-12}
                    if coefs:
                        vecs.append(coefs)
            empty = (ci == 0 and delta0 is not None)
            off = 1 if empty else 0
            ent = {}
            for i in range(len(vecs)):
                for j in range(i, len(vecs)):
                    acc = {}; cst = 0.0
                    for xa, ca in vecs[i].items():
                        for xb, cb in vecs[j].items():
                            d_, c_ = L(xa, xb)
                            cst += ca * cb * c_
                            for k_, v_ in d_.items(): acc[k_] = acc.get(k_, 0) + ca * cb * v_
                    ent[i + off, j + off] = (acc, cst)
            if empty:
                ent[0, 0] = ({}, k / delta0 - corner)
                for j in range(len(vecs)):
                    acc = {}; cst = 0.0
                    for xb, cb in vecs[j].items():
                        d_, c_ = L(xb, xb)
                        cst += cb * c_
                        for k_, v_ in d_.items(): acc[k_] = acc.get(k_, 0) + cb * v_
                    ent[0, j + 1] = (acc, cst)
            out.append((len(vecs) + off, ent)); allvecs.append(vecs)
        loc_vectors.extend(allvecs)
        return out
    loc_vectors = []                                  # the basis vectors of each localizing block
    mats = orthogonal_group(q, n0)
    det = lambda R: (R[0, 0] * R[1, 1] - R[0, 1] * R[1, 0]) % q
    if edge:
        e = int(U_[0]); ex_, ey_ = int(X[e]), int(Y[e])
        te = next(R for R in mats if det(R) != 1 and apply(R, ex_, ey_) == (ex_, ey_))
        h1 = idx(*apply(te, X, Y)); h2 = idx(ex_ - X, ey_ - Y)
        group = [np.arange(n), h1, h2, h2[h1]]
        table = [(1, 1, 1, 1), (1, 1, -1, -1), (1, -1, 1, -1), (1, -1, -1, 1)]
        chars = [[(lambda gi, ch=ch: ch[gi])] for ch in table]
        blocks += localizing([0, e], group, chars, 2.0)
    if trilocal:
        e = int(U_[0]); f = next(int(u) for u in U_ if sub(u, e) in Us_)
        T3 = [0, e, f]
        ex_, ey_, fx_, fy_ = int(X[e]), int(Y[e]), int(X[f]), int(Y[f])
        dinv = pow((ex_ * fy_ - ey_ * fx_) % q, q - 2, q)
        Binv = np.array([[fy_, -fx_], [-ey_, ex_]]) * dinv % q
        group, perms = [], []
        for perm in itertools.permutations(range(3)):
            im = [T3[i] for i in perm]
            t0x, t0y = int(X[im[0]]), int(Y[im[0]])
            Mimg = np.array([[(int(X[im[1]]) - t0x) % q, (int(X[im[2]]) - t0x) % q],
                             [(int(Y[im[1]]) - t0y) % q, (int(Y[im[2]]) - t0y) % q]])
            Rl = (Mimg @ Binv) % q
            for (x_, y_) in [(1, 0), (0, 1), (1, 1)]:
                assert Qf(*apply(Rl, x_, y_)) == Qf(x_, y_)
            group.append(idx(Rl[0, 0] * X + Rl[0, 1] * Y + t0x, Rl[1, 0] * X + Rl[1, 1] * Y + t0y))
            perms.append(perm)
        sgn = [1 if sum(p[i] > p[j] for i in range(3) for j in range(i + 1, 3)) % 2 == 0 else -1 for p in perms]
        # standard representation on {x in R^3 : sum x = 0}, with w_i = 3 e_i - (1, 1, 1): the G-maps
        # phi_i(w_0) = sum_h <w_0, h w_i> / <w_0, w_0> e_{h p}, i = 0, 1, span the multiplicity space;
        # <w_0, h w_i> / <w_0, w_0> is 1 if h sends vertex i to vertex 0, else -1/2
        r0 = [1.0 if p[0] == 0 else -0.5 for p in perms]
        r1 = [1.0 if p[1] == 0 else -0.5 for p in perms]
        chars = [[lambda gi: 1.0], [lambda gi: float(sgn[gi])], [lambda gi: r0[gi], lambda gi: r1[gi]]]
        blocks += localizing(T3, group, chars, 3.0)
    if pentagon:
        # an induced unit pentagon T through 0: 2 - sum_{t in T} x_t >= 0; symmetry: the isometries preserving
        # T that fix it pointwise or reverse it (a reflection at most is used: characters +-1)
        Pts = None
        for u1 in U_:
            for u2 in U_:
                p2 = idx(X[u1] + X[u2], Y[u1] + Y[u2])
                if p2 == 0 or cls[p2] == 1: continue
                for u3 in U_:
                    p3 = idx(X[p2] + X[u3], Y[p2] + Y[u3])
                    if p3 in (0, u1) or cls[p3] == 1 or cls[sub(p3, u1)] == 1: continue
                    for p4 in U_:                  # p4 adjacent to 0, to p3; not to u1, p2
                        if p4 in (u1, p2, p3) or cls[sub(p4, p3)] != 1 or cls[sub(p4, u1)] == 1 or cls[sub(p4, p2)] == 1:
                            continue
                        Pts = [0, int(u1), int(p2), int(p3), int(p4)]; break
                    if Pts: break
                if Pts: break
            if Pts: break
        assert Pts is not None, "no induced unit pentagon"
        # isometries mapping T to itself: try the reflection reversing the cycle 0 -> 0, u1 <-> p4, p2 <-> p3
        group = [np.arange(n)]
        for Rm in mats:
            if det(Rm) == 1: continue
            if all(idx(*apply(Rm, int(X[a_]), int(Y[a_]))) == b_ for a_, b_ in ((Pts[1], Pts[4]), (Pts[2], Pts[3]))):
                group.append(idx(*apply(Rm, X, Y))); break
        chars = [[lambda gi: 1.0]] if len(group) == 1 else [[lambda gi: 1.0], [lambda gi: (1.0, -1.0)[gi]]]
        blocks += localizing(Pts, group, chars, 5.0, k=2)
    if len(blocks) > nrot and verbose:
        print(f"localizing blocks of sizes {[b[0] for b in blocks[nrot:]]} [{time.time()-t0:.0f}s]", flush=True)
    meta = dict(q=q, kind=kind, N=N, m=m, reps=[int(r) for r in reps], partner=[int(p) for p in partner], shift=shift.tolist(),
                sig=sig.tolist(), tau=tau.tolist(), n0=int(n0), K0=K0, nrot=nrot,
                options=dict(tri=tri, trilocal=trilocal, edge=edge, pentagon=pentagon, colour=colour))
    return dict(nv=nv, nz=nz, free=free, gcol=gcol, blocks=blocks, lin=lin, size_c=size_c, Pm=Pm, n=n,
                meta=meta, lab=lab, zid=zid, unit=unit, cls=cls, t0=t0, localizing_vectors=loc_vectors)


# ------------------------------------------------------------------------------------------ solvers
def conic_form(D):
    """rows of A v + s = b with s in (nonnegative cone) x (psd cones, svec with sqrt 2 off the diagonal)"""
    nv = D['nv']; rows, cols, vals, b = [], [], [], []
    r = 0
    for i in range(nv):
        rows.append(r + i); cols.append(i); vals.append(-1.0)
    b += [0.0] * nv; r += nv
    for d, rhs in D['lin']:
        for i, c in d.items():
            rows.append(r); cols.append(int(i)); vals.append(float(c))
        b.append(float(rhs)); r += 1
    for xi in range(D['Pm'].shape[0]):
        for c in D['free']:
            rows.append(r); cols.append(D['gcol'][c]); vals.append(-float(D['Pm'][xi, c]))
        b.append(float(D['Pm'][xi, 0])); r += 1
    nonneg = r
    s2 = np.sqrt(2)
    for size, ent in D['blocks']:
        for j in range(size):
            for i in range(j + 1):
                acc, cst = ent.get((i, j), ({}, 0.0))
                sc = 1.0 if i == j else s2
                for vi, c in acc.items():
                    rows.append(r); cols.append(int(vi)); vals.append(-sc * float(np.real(c)))
                b.append(sc * float(np.real(cst))); r += 1
    A = sp.csc_matrix((vals, (rows, cols)), shape=(r, nv))
    c = np.zeros(nv)
    for cc in D['free']:
        c[D['gcol'][cc]] = -float(D['size_c'][cc])
    return A, np.array(b), c, nonneg


def solve_clarabel(D, verbose=False, eps=1e-8):
    import clarabel
    A, b, c, nonneg = conic_form(D)
    cones = [clarabel.NonnegativeConeT(nonneg)] + [clarabel.PSDTriangleConeT(s) for s, _ in D['blocks']]
    st = clarabel.DefaultSettings()
    st.verbose = verbose; st.tol_gap_abs = eps; st.tol_gap_rel = eps; st.tol_feas = eps; st.max_iter = 400
    sol = clarabel.DefaultSolver(sp.csc_matrix((D['nv'], D['nv'])), c, A, b, cones, st).solve()
    return 1 - sol.obj_val, np.array(sol.z), str(sol.status)


def solve_cvxopt(D, verbose=False):
    """cvxopt's interior-point method: little memory (one dense Schur matrix of size nv)"""
    from cvxopt import matrix, spmatrix, solvers
    A, b, c, nonneg = conic_form(D)
    nv = D['nv']
    Al = A[:nonneg].tocoo()
    Gl = spmatrix(Al.data.tolist(), Al.row.tolist(), Al.col.tolist(), (nonneg, nv))
    Gs, hs = [], []
    for size, ent in D['blocks']:
        I_, J_, V_ = [], [], []
        H0 = np.zeros((size, size))
        for (i, j), (acc, cst) in ent.items():
            for (a, bb) in (((i, j), (j, i)) if i != j else ((i, j),)):
                H0[a, bb] = float(np.real(cst))
                for vi, cc in acc.items():
                    I_.append(int(bb * size + a)); J_.append(int(vi)); V_.append(-float(np.real(cc)))
        Gs.append(spmatrix(V_, I_, J_, (size * size, nv))); hs.append(matrix(H0))
    solvers.options.update(show_progress=verbose, abstol=1e-9, reltol=1e-9, feastol=1e-9, maxiters=200)
    sol = solvers.sdp(matrix(c), Gl=Gl, hl=matrix(b[:nonneg]), Gs=Gs, hs=hs)
    parts = [np.array(sol['zl']).ravel()]
    for Zm, (size, _) in zip(sol['zs'], D['blocks']):
        Z = np.array(Zm); Z = (Z + Z.T) / 2
        parts.append(np.array([Z[i, j] if i == j else np.sqrt(2) * Z[i, j] for j in range(size) for i in range(j + 1)]))
    return 1 - sol['primal objective'], np.concatenate(parts), sol['status']


def block_matrices(D):
    """for each block: (sparse (size^2, nv) coefficient matrix, constant vector), row-major vec"""
    out = []
    for size, ent in D['blocks']:
        rows, cols, vals = [], [], []
        c0 = np.zeros(size * size)
        for (i, j), (acc, cst) in ent.items():
            for (a, b) in (((i, j), (j, i)) if i != j else ((i, j),)):
                c0[a * size + b] = np.real(cst)
                for vi, c in acc.items():
                    rows.append(a * size + b); cols.append(int(vi)); vals.append(float(np.real(c)))
        out.append((sp.csr_matrix((vals, (rows, cols)), shape=(size * size, D['nv'])), c0))
    return out


def solve_primal_dsdp(D):
    """DSDP (through cvxopt): accurate in the primal, not in the dual"""
    from cvxopt import matrix, spmatrix, solvers
    A, b, c, nonneg = conic_form(D)
    nv = D['nv']
    Al = A[:nonneg].tocoo()
    Gl = spmatrix(Al.data.tolist(), Al.row.tolist(), Al.col.tolist(), (nonneg, nv))
    Gs, hs = [], []
    for size, ent in D['blocks']:
        I_, J_, V_ = [], [], []
        H0 = np.zeros((size, size))
        for (i, j), (acc, cst) in ent.items():
            for (a, bb) in (((i, j), (j, i)) if i != j else ((i, j),)):
                H0[a, bb] = float(np.real(cst))
                for vi, cc in acc.items():
                    I_.append(int(bb * size + a)); J_.append(int(vi)); V_.append(-float(np.real(cc)))
        Gs.append(spmatrix(V_, I_, J_, (size * size, nv))); hs.append(matrix(H0))
    solvers.options.update(show_progress=False, dsdp_Monitor=0, dsdp_GapTolerance=1e-9)
    sol = solvers.sdp(matrix(c), Gl=Gl, hl=matrix(b[:nonneg]), Gs=Gs, hs=hs, solver='dsdp')
    Xs = [np.array(Zm) for Zm in sol['zs']]
    return np.array(sol['x']).ravel(), 1 - sol['primal objective'], sol['status'], Xs


def dual_lp(D, v, tol=1e-5, max_dirs=6000, verbose=True, Xs=None, active_tol=1e-5):
    """a dual solution from a primal optimum v: the psd multipliers are nonnegative combinations of u u^T over
    the null vectors u of each block at v (and their sums and differences, while the total number of
    directions stays below max_dirs), and the linear multipliers are free; one linear programme"""
    nv = D['nv']
    F = np.zeros(nv)
    for c in D['free']: F[D['gcol'][c]] = float(D['size_c'][c])
    mats = block_matrices(D)
    nulls = []
    for bi, ((Bk, c0), (size, _)) in enumerate(zip(mats, D['blocks'])):
        Mv = (Bk @ v + c0).reshape(size, size); Mv = (Mv + Mv.T) / 2
        lam, V = np.linalg.eigh(Mv)
        N0 = V[:, lam < tol * max(1.0, np.abs(lam).max())]
        if Xs is not None and N0.shape[1] > 0:
            # the solver's dual matrix, compressed to the null space: its range gives the directions
            X = (Xs[bi] + Xs[bi].T) / 2
            Xn = N0.T @ X @ N0
            mu, Wv = np.linalg.eigh((Xn + Xn.T) / 2)
            keep = mu > 1e-7 * max(1e-300, mu.max())
            N0 = N0 @ Wv[:, keep] if keep.any() else N0[:, :0]
        nulls.append(N0)
    npairs = sum(U.shape[1] ** 2 for U in nulls)
    use_pairs = npairs <= max_dirs
    consts, coefs, owner = [], [], []
    for bi, ((Bk, c0), U) in enumerate(zip(mats, nulls)):
        d = U.shape[1]
        dirs = [U[:, i] for i in range(d)]
        if use_pairs:
            for i in range(d):
                for j in range(i + 1, d):
                    dirs += [(U[:, i] + U[:, j]) / np.sqrt(2), (U[:, i] - U[:, j]) / np.sqrt(2)]
        if not dirs:
            continue
        W = np.array([np.outer(u, u).ravel() for u in dirs])             # directions x size^2
        consts.append(W @ c0); coefs.append(np.asarray((Bk.T @ W.T)).T)   # directions x nv
        owner += [(bi, u) for u in dirs]
    Cl = np.concatenate(consts) if consts else np.zeros(0)
    Gl = np.vstack(coefs) if coefs else np.zeros((0, nv))
    nlam = len(owner)
    # complementary slackness: only the linear constraints tight at v can carry a multiplier
    active = [l for l, (d, rhs) in enumerate(D['lin'])
              if rhs - sum(c * v[i] for i, c in d.items()) < active_tol * (1 + abs(rhs))]
    lin = [D['lin'][l] for l in active]; nlin = len(lin); Pm = D['Pm']; nch = Pm.shape[0]
    rows, cols, vals = [], [], []
    for l, (d, rhs) in enumerate(lin):
        for i, c in d.items(): rows.append(int(i)); cols.append(l); vals.append(-float(c))
    Ay = sp.csr_matrix((vals, (rows, cols)), shape=(nv, nlin))
    Aw = np.zeros((nv, nch))
    for c in D['free']: Aw[D['gcol'][c]] = Pm[:, c]
    Aub = sp.hstack([sp.csr_matrix(Gl.T), Ay, sp.csr_matrix(Aw), -sp.identity(nv)]).tocsr()
    obj = np.concatenate([Cl, [float(r) for _, r in lin], Pm[:, 0], np.ones(nv)])
    t0 = time.time()
    res = linprog(obj, A_ub=Aub, b_ub=-F, bounds=[(0, None)] * len(obj), method='highs')
    if verbose:
        print(f"dual LP: null spaces {[U.shape[1] for U in nulls]}, {nlam} directions "
              f"({'with' if use_pairs else 'without'} pairs), {nlin} of {len(D['lin'])} linear constraints active, "
              f"status {res.status}, bound {1 + res.fun:.6f} [{time.time()-t0:.0f}s]", flush=True)
    lamv = res.x[:nlam]; w = res.x[nlam + nlin:nlam + nlin + nch]
    y = np.zeros(len(D['lin'])); y[active] = res.x[nlam:nlam + nlin]
    Z = [np.zeros((s, s)) for s, _ in D['blocks']]
    for lv, (bi, u) in zip(lamv, owner):
        if lv > 0: Z[bi] += lv * np.outer(u, u)
    parts = [np.zeros(nv), y, w]
    for (s, _), Zk in zip(D['blocks'], Z):
        parts.append(np.array([Zk[i, j] if i == j else np.sqrt(2) * Zk[i, j] for j in range(s) for i in range(j + 1)]))
    return 1 + res.fun, np.concatenate(parts)


def dual_nullsdp(D, v, tol=1e-5, verbose=True, active_tol=1e-5):
    """a dual solution from a primal optimum v: each psd multiplier is U_k W_k U_k^T with U_k a basis of the
    null space of block k at v and W_k >= 0 a small psd matrix; with the multipliers y, w of the (active)
    linear constraints this is a small semidefinite programme, solved by Clarabel"""
    import clarabel
    nv = D['nv']
    F = np.zeros(nv)
    for c in D['free']: F[D['gcol'][c]] = float(D['size_c'][c])
    mats = block_matrices(D)
    Ucols, Gcols, consts, dims, owners = [], [], [], [], []
    for bi, ((Bk, c0), (size, _)) in enumerate(zip(mats, D['blocks'])):
        Mv = (Bk @ v + c0).reshape(size, size); Mv = (Mv + Mv.T) / 2
        lam, V = np.linalg.eigh(Mv)
        U = V[:, lam < tol * max(1.0, np.abs(lam).max())]
        d = U.shape[1]
        if d == 0:
            continue
        # svec basis of d x d symmetric matrices: E_ii, (E_ij + E_ji)/sqrt 2; tr(W G) = <svec W, svec G>
        basis = []
        for j in range(d):
            for i in range(j + 1):
                if i == j: basis.append(np.outer(U[:, i], U[:, i]))
                else: basis.append((np.outer(U[:, i], U[:, j]) + np.outer(U[:, j], U[:, i])) / np.sqrt(2))
        Wm = np.array([b_.ravel() for b_ in basis])                        # svec-dim x size^2
        consts.append(Wm @ c0); Gcols.append(np.asarray((Bk.T @ Wm.T)).T)  # svec-dim x nv
        dims.append(d); owners.append((bi, U))
    active = [l for l, (dd, rhs) in enumerate(D['lin'])
              if rhs - sum(c * v[i] for i, c in dd.items()) < active_tol * (1 + abs(rhs))]
    lin = [D['lin'][l] for l in active]; nlin = len(lin); Pm = D['Pm']; nch = Pm.shape[0]
    nW = sum(d * (d + 1) // 2 for d in dims)
    Gw = np.vstack(Gcols) if Gcols else np.zeros((0, nv))
    rows, cols, vals = [], [], []
    for l, (dd, rhs) in enumerate(lin):
        for i, c in dd.items(): rows.append(int(i)); cols.append(l); vals.append(-float(c))
    Ay = sp.csr_matrix((vals, (rows, cols)), shape=(nv, nlin))
    Aw = np.zeros((nv, nch))
    for c in D['free']: Aw[D['gcol'][c]] = Pm[:, c]
    nx = nW + nlin + nch + nv
    # rows: (a) coef_i - t_i <= -F_i  ->  [Gw^T, Ay, Aw, -I] x + s = -F, s >= 0
    #       (b) y, w, t >= 0          ->  -x_{y,w,t} + s = 0, s >= 0
    #       (c) W_k psd               ->  -svec(W_k) + s = 0, s in psd
    Aa = sp.hstack([sp.csr_matrix(Gw.T), Ay, sp.csr_matrix(Aw), -sp.identity(nv)])
    Ab = sp.hstack([sp.csr_matrix((nlin + nch + nv, nW)), -sp.identity(nlin + nch + nv)])
    Ac = sp.hstack([-sp.identity(nW), sp.csr_matrix((nW, nlin + nch + nv))])
    A = sp.vstack([Aa, Ab, Ac]).tocsc()
    b = np.concatenate([-F, np.zeros(nlin + nch + nv), np.zeros(nW)])
    cvec = np.concatenate([np.concatenate(consts) if consts else np.zeros(0),
                           [float(r) for _, r in lin], Pm[:, 0], np.ones(nv)])
    cones = [clarabel.NonnegativeConeT(nv + nlin + nch + nv)] + [clarabel.PSDTriangleConeT(d) for d in dims]
    st = clarabel.DefaultSettings(); st.verbose = False
    st.tol_gap_abs = 1e-9; st.tol_gap_rel = 1e-9; st.tol_feas = 1e-9; st.max_iter = 300
    t0 = time.time()
    sol = clarabel.DefaultSolver(sp.csc_matrix((nx, nx)), cvec, A, b, cones, st).solve()
    x = np.array(sol.x)
    if verbose:
        print(f"null-space dual: null spaces {dims}, {nlin} of {len(D['lin'])} linear constraints active, "
              f"{sol.status}, bound {1 + sol.obj_val:.6f} [{time.time()-t0:.0f}s]", flush=True)
    Z = [np.zeros((s, s)) for s, _ in D['blocks']]
    off = 0
    for d, (bi, U) in zip(dims, owners):
        sv = x[off:off + d * (d + 1) // 2]; off += d * (d + 1) // 2
        Wk = np.zeros((d, d)); t = 0
        for j in range(d):
            for i in range(j + 1):
                if i == j: Wk[i, i] = sv[t]
                else: Wk[i, j] = Wk[j, i] = sv[t] / np.sqrt(2)
                t += 1
        lam_, V_ = np.linalg.eigh(Wk)
        Wk = (V_ * np.maximum(lam_, 0)) @ V_.T
        Z[bi] = U @ Wk @ U.T
    y = np.zeros(len(D['lin'])); y[active] = np.maximum(x[nW:nW + nlin], 0)
    w = np.maximum(x[nW + nlin:nW + nlin + nch], 0)
    parts = [np.zeros(nv), y, w]
    for (s, _), Zk in zip(D['blocks'], Z):
        parts.append(np.array([Zk[i, j] if i == j else np.sqrt(2) * Zk[i, j] for j in range(s) for i in range(j + 1)]))
    return 1 + sol.obj_val, np.concatenate(parts)


def solve_dsdp_lp(D, verbose=False, tol=1e-5):
    """primal by DSDP, dual by a small semidefinite programme on the null spaces of the blocks"""
    v, val, st, Xs = solve_primal_dsdp(D)
    if verbose: print(f"DSDP: value {val:.6f} ({st})", flush=True)
    best = None
    for t in (tol, tol * 10, tol * 100):
        bound, z = dual_nullsdp(D, v, tol=t, verbose=verbose)
        if best is None or bound < best[0]:
            best = (bound, z)
        if bound < val + 1e-3 * max(1.0, val):
            break
    bound, z = best
    return bound, z, f"DSDP {st}; dual on the null spaces"


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('q', type=int); ap.add_argument('kind', choices=['std', 'inert'])
    ap.add_argument('--solver', default='dsdp', choices=['dsdp', 'clarabel', 'cvxopt'])
    ap.add_argument('--tri', action='store_true'); ap.add_argument('--trilocal', action='store_true')
    ap.add_argument('--edge', action='store_true'); ap.add_argument('--test', type=int, default=None)
    ap.add_argument('--pentagon', action='store_true')
    ap.add_argument('--colour', type=int, default=None)
    ap.add_argument('--save', default=None); ap.add_argument('--verbose', action='store_true')
    a = ap.parse_args()
    D = build(a.q, a.kind, K0=a.test, tri=a.tri, trilocal=a.trilocal, edge=a.edge, pentagon=a.pentagon,
              colour=a.colour)
    solve = {'dsdp': solve_dsdp_lp, 'clarabel': solve_clarabel, 'cvxopt': solve_cvxopt}[a.solver]
    val, z, status = solve(D, verbose=a.verbose)
    n = a.q * a.q
    if a.colour:
        verdict = "so no proper colouring" if val < n / a.colour else "no conclusion"
        print(f"{a.kind}{a.q}: {a.colour}-colouring programme, E|colour class| <= {val:.4f} against n/{a.colour} = "
              f"{n / a.colour:.2f}: {verdict}; {status} [{time.time() - D['t0']:.0f}s]", flush=True)
    else:
        print(f"{a.kind}{a.q}: three-point bound alpha <= {val:.4f} ({val / n:.4f} n); n/5 = {n / 5:.1f}, "
              f"n/6 = {n / 6:.1f}; {status} [{time.time() - D['t0']:.0f}s]", flush=True)
    if a.save:
        np.savez_compressed(a.save, z=z, meta=json.dumps(D['meta']), value=val, status=status)
        print(f"dual solution written to {a.save}")
