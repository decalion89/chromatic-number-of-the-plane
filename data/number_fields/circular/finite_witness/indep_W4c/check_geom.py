#!/usr/bin/env python3
"""Independent referee checker for claims 1-5 (and the unit-distance count).

Usage:  python3 check_geom.py WITNESS.json.gz [--no-unit] [--json OUT]

All decisions are made with exact integer / Fraction arithmetic.
Floating point is used only for an optional cross-check that is reported
but never decides anything.
"""
import gzip
import json
import math
import os
import sys
import time
from collections import deque
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from biquad import Field, int_sq_norm_8, vec_to_int8  # noqa: E402

FIELDS = {
    "Q(sqrt3, sqrt11)": (3, 11, "(1, sqrt3, sqrt11, sqrt33)"),
    "Q(sqrt2, sqrt3)": (2, 3, "(1, sqrt2, sqrt3, sqrt6)"),
}


def load(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        return json.load(fh)


def is_int(x):
    return isinstance(x, int) and not isinstance(x, bool)


# --------------------------------------------------------------------------
# Construction of the expected generator sets (claims 1 and 2)
# --------------------------------------------------------------------------
def F_(q):
    return Fraction(q)


def elt(c0=0, c1=0, c2=0, c3=0):
    return (F_(c0), F_(c1), F_(c2), F_(c3))


def cpow(F, z, k):
    """z^k for a complex number z = (re, im) of norm 1 (k may be negative;
    the inverse of a norm-1 element is its conjugate, which we verify)."""
    one = (elt(1), elt(0))
    if k < 0:
        assert F.sq_norm(z) == elt(1), "inverse via conjugate needs norm 1"
        z = (z[0], F.neg(z[1]))
        k = -k
    r = one
    for _ in range(k):
        r = F.cmul(r, z)
    return r


def expected_q3_11(F, D):
    """The 27 vectors R60^j RA^k RG^l (1,0), j,k,l in {-1,0,1}, as integer
    8-tuples over D, plus float cross-check data."""
    # rotations as complex numbers cos + i sin; basis (1, sqrt3, sqrt11, sqrt33)
    R60 = (elt(Fraction(1, 2)), elt(0, Fraction(1, 2)))          # cos 1/2, sin sqrt3/2
    RA = (elt(Fraction(5, 6)), elt(0, 0, Fraction(1, 6)))         # cos 5/6, sin sqrt11/6
    RG = (elt(Fraction(11, 14)), elt(0, Fraction(5, 14)))         # cos 11/14, sin 5sqrt3/14
    info = {}
    for name, R in (("R60", R60), ("RA", RA), ("RG", RG)):
        n = F.sq_norm(R)
        info[name + "_norm_is_1"] = (n == elt(1))
        if n != elt(1):
            raise AssertionError(name + " is not a rotation")
    # float cross-check of the angles
    info["float_angles_deg"] = {
        "R60": math.degrees(math.atan2(F.to_float(R60[1]), F.to_float(R60[0]))),
        "RA": math.degrees(math.atan2(F.to_float(RA[1]), F.to_float(RA[0]))),
        "RG": math.degrees(math.atan2(F.to_float(RG[1]), F.to_float(RG[0]))),
    }
    e1 = (elt(1), elt(0))  # the vector (1, 0)
    vecs = {}
    for j in (-1, 0, 1):
        for k in (-1, 0, 1):
            for l in (-1, 0, 1):
                # apply RG^l first, then RA^k, then R60^j (they commute anyway)
                v = F.cmul(cpow(F, RG, l), e1)
                v = F.cmul(cpow(F, RA, k), v)
                v = F.cmul(cpow(F, R60, j), v)
                vecs[(j, k, l)] = vec_to_int8(v, D)
    return vecs, info


def expected_q2_3(F, D):
    """The 120 vectors zeta^j w^l, j in Z/24, l in {-2..2}."""
    # basis (1, sqrt2, sqrt3, sqrt6)
    zeta = (elt(0, Fraction(1, 4), 0, Fraction(1, 4)),       # (sqrt6 + sqrt2)/4
            elt(0, Fraction(-1, 4), 0, Fraction(1, 4)))      # (sqrt6 - sqrt2)/4
    w = (elt(Fraction(1, 3)), elt(0, Fraction(2, 3)))        # (1/3, 2 sqrt2/3)
    info = {}
    info["zeta_norm_is_1"] = F.sq_norm(zeta) == elt(1)
    info["w_norm_is_1"] = F.sq_norm(w) == elt(1)
    if not (info["zeta_norm_is_1"] and info["w_norm_is_1"]):
        raise AssertionError("zeta or w is not of norm 1")
    one = (elt(1), elt(0))
    # zeta^2 must be e^{i pi/6} = (sqrt3/2, 1/2); with Re zeta > 0 this pins zeta = e^{i pi/12}
    z2 = F.cmul(zeta, zeta)
    info["zeta^2 == (sqrt3/2, 1/2)"] = (z2 == (elt(0, 0, Fraction(1, 2)), elt(Fraction(1, 2))))
    powers = [one]
    for _ in range(24):
        powers.append(F.cmul(powers[-1], zeta))
    info["zeta^12 == -1"] = powers[12] == (elt(-1), elt(0))
    info["zeta^24 == 1"] = powers[24] == one
    info["zeta primitive 24th root"] = all(powers[k] != one for k in range(1, 24))
    info["float zeta angle deg"] = math.degrees(math.atan2(F.to_float(zeta[1]), F.to_float(zeta[0])))
    info["float w"] = (F.to_float(w[0]), F.to_float(w[1]))
    # w = (1 + 2 sqrt(-2))/3 : real part 1/3, imaginary part 2 sqrt2/3
    vecs = {}
    for j in range(24):
        for l in (-2, -1, 0, 1, 2):
            v = F.cmul(powers[j], cpow(F, w, l))
            vecs[(j, l)] = vec_to_int8(v, D)
    return vecs, info


def compare_up_to_sign(gens, expected_vecs, expect_closed_under_neg):
    """Compare the generator list with the constructed vectors, up to sign.
    Returns (ok, details)."""
    neg = lambda t: tuple(-x for x in t)
    det = {}
    gpairs = set(frozenset((tuple(g), neg(tuple(g)))) for g in gens)
    det["num_generators"] = len(gens)
    det["num_distinct_generator_pairs"] = len(gpairs)
    ev = list(expected_vecs.values())
    det["num_constructed"] = len(ev)
    det["num_distinct_constructed_vectors"] = len(set(ev))
    epairs = set(frozenset((v, neg(v))) for v in ev)
    det["num_distinct_constructed_pairs"] = len(epairs)
    if expect_closed_under_neg:
        det["constructed_closed_under_negation"] = all(neg(v) in set(ev) for v in ev)
    det["generator_pairs_minus_constructed"] = len(gpairs - epairs)
    det["constructed_pairs_minus_generator"] = len(epairs - gpairs)
    ok = (len(gpairs) == len(gens) and gpairs == epairs and len(epairs) == len(gens))
    if expect_closed_under_neg:
        ok = ok and det["constructed_closed_under_negation"] and len(set(ev)) == 2 * len(gens)
    else:
        ok = ok and len(epairs) == len(ev)
    det["ok"] = ok
    return ok, det


# --------------------------------------------------------------------------
# The main check
# --------------------------------------------------------------------------
def check_all(data, do_unit=True, do_construction=True, log=print):
    res = {}      # name -> bool
    stats = {}
    t0 = time.time()

    # ---- format ----
    field = data.get("field")
    if field not in FIELDS:
        res["field_known"] = False
        return res, stats
    A, B, basis = FIELDS[field]
    res["basis_string_matches"] = (data.get("basis") == basis)
    F = Field(A, B)
    D = data["denominator"]
    res["denominator_positive_int"] = is_int(D) and D > 0
    gens = data["generators"]
    pts = data["points"]
    edges = data["edges"]
    col = data["colouring"]
    cycles = data["cycles"]
    fv = data["fixed_vertex"]
    n = len(pts)
    stats["n_points"] = n
    stats["n_edges_listed"] = len(edges)
    stats["n_generators"] = len(gens)
    stats["n_cycles"] = len(cycles)

    res["generators_are_int8"] = all(isinstance(g, list) and len(g) == 8 and all(is_int(x) for x in g) for g in gens)
    res["points_are_int8"] = all(isinstance(p, list) and len(p) == 8 and all(is_int(x) for x in p) for p in pts)
    P = [tuple(p) for p in pts]
    res["points_distinct"] = (len(set(P)) == n)
    stats["n_distinct_points"] = len(set(P))

    # ---- claims 1/2: generator construction ----
    if do_construction:
        if field == "Q(sqrt3, sqrt11)":
            ev, info = expected_q3_11(F, D)
            ok, det = compare_up_to_sign(gens, ev, expect_closed_under_neg=False)
        else:
            ev, info = expected_q2_3(F, D)
            ok, det = compare_up_to_sign(gens, ev, expect_closed_under_neg=True)
        stats["construction_info"] = info
        stats["construction_compare"] = det
        res["claim12_generators_match_construction"] = ok and all(v for k, v in info.items() if isinstance(v, bool))

    # ---- claim 3a: generators have length exactly 1 ----
    unit = (D * D, 0, 0, 0)
    res["claim3_generators_unit_length"] = all(int_sq_norm_8(F, g) == unit for g in gens)
    neg = lambda t: tuple(-x for x in t)
    Gpm = set()
    for g in gens:
        Gpm.add(tuple(g))
        Gpm.add(neg(tuple(g)))
    res["claim3_no_zero_generator"] = (0,) * 8 not in Gpm
    stats["n_signed_generators"] = len(Gpm)

    # ---- claim 3b: listed edges == Cayley edges among the points ----
    idx = {p: i for i, p in enumerate(P)}
    cay = set()
    for i, p in enumerate(P):
        for g in Gpm:
            q = tuple(a + b for a, b in zip(p, g))
            j = idx.get(q)
            if j is not None and j != i:
                cay.add((min(i, j), max(i, j)))
    stats["n_cayley_pairs_by_generator_search"] = len(cay)
    edges_ok_format = all(isinstance(e, list) and len(e) == 2 and all(is_int(x) for x in e)
                          and 0 <= e[0] < n and 0 <= e[1] < n and e[0] != e[1] for e in edges)
    res["edges_well_formed"] = edges_ok_format
    E = set()
    dup = 0
    for e in edges:
        if not (isinstance(e, list) and len(e) == 2):
            continue
        k = (min(e), max(e))
        if k in E:
            dup += 1
        E.add(k)
    stats["n_duplicate_edges"] = dup
    res["edges_no_duplicates"] = (dup == 0)
    not_cay = [e for e in E if e not in cay]
    missing = [e for e in cay if e not in E]
    stats["n_listed_edges_not_cayley"] = len(not_cay)
    stats["n_cayley_pairs_not_listed"] = len(missing)
    if not_cay:
        stats["example_listed_not_cayley"] = sorted(not_cay)[:5]
    if missing:
        stats["example_cayley_not_listed"] = sorted(missing)[:5]
    # direct check of every listed edge's difference (independent of the hash search)
    bad_diff = 0
    for (i, j) in E:
        if 0 <= i < n and 0 <= j < n:
            d = tuple(b - a for a, b in zip(P[i], P[j]))
            if d not in Gpm:
                bad_diff += 1
        else:
            bad_diff += 1
    stats["n_listed_edges_with_bad_difference"] = bad_diff
    res["claim3_listed_edges_are_generator_differences"] = (bad_diff == 0 and len(not_cay) == 0)
    res["claim3_every_cayley_pair_is_listed"] = (len(missing) == 0)

    # ---- claim 3c: connectivity from the fixed vertex; fixed vertex = origin ----
    res["fixed_vertex_in_range"] = is_int(fv) and 0 <= fv < n
    if res["fixed_vertex_in_range"]:
        res["claim3_fixed_vertex_is_origin"] = (P[fv] == (0,) * 8)
        adj = [[] for _ in range(n)]
        for (i, j) in E:
            if 0 <= i < n and 0 <= j < n:
                adj[i].append(j)
                adj[j].append(i)
        seen = [False] * n
        seen[fv] = True
        dq = deque([fv])
        while dq:
            u = dq.popleft()
            for v in adj[u]:
                if not seen[v]:
                    seen[v] = True
                    dq.append(v)
        stats["n_reachable_from_fixed_vertex"] = sum(seen)
        res["claim3_all_points_reachable"] = all(seen)
        degs = [len(a) for a in adj]
        stats["degree_min_max"] = (min(degs), max(degs))
    else:
        res["claim3_fixed_vertex_is_origin"] = False
        res["claim3_all_points_reachable"] = False

    # ---- claim 4: proper 4-colouring ----
    res["colouring_length_n"] = (len(col) == n)
    res["colouring_values_in_0123"] = all(is_int(c) and 0 <= c <= 3 for c in col)
    bad_col = 0
    if res["colouring_length_n"]:
        for (i, j) in E:
            if 0 <= i < n and 0 <= j < n and col[i] == col[j]:
                bad_col += 1
    stats["n_monochromatic_edges"] = bad_col
    res["claim4_proper_4_colouring"] = (res["colouring_length_n"] and res["colouring_values_in_0123"] and bad_col == 0)
    if res["colouring_length_n"]:
        stats["colour_class_sizes"] = [sum(1 for c in col if c == k) for k in range(4)]

    # ---- claim 5: cycles ----
    bad = {"format": 0, "empty": 0, "len_not_div4": 0, "repeated_vertex": 0, "non_edge_step": 0}
    lens = {}
    seen_cyc = set()
    dup_cyc = 0
    rev_pairs = 0
    for cyc in cycles:
        if not (isinstance(cyc, list) and all(is_int(v) and 0 <= v < n for v in cyc)):
            bad["format"] += 1
            continue
        m = len(cyc)
        lens[m] = lens.get(m, 0) + 1
        if m == 0:
            bad["empty"] += 1
            continue
        if m % 4 != 0:
            bad["len_not_div4"] += 1
        if len(set(cyc)) != m:
            bad["repeated_vertex"] += 1
        for t in range(m):
            u, v = cyc[t], cyc[(t + 1) % m]
            if (min(u, v), max(u, v)) not in E:
                bad["non_edge_step"] += 1
                break
        key = tuple(cyc)
        if key in seen_cyc:
            dup_cyc += 1
        seen_cyc.add(key)
    # count cycles whose reversal (as a cyclic sequence) is also listed
    canon = set()
    for cyc in cycles:
        if isinstance(cyc, list) and cyc:
            m = len(cyc)
            rots = [tuple(cyc[s:] + cyc[:s]) for s in range(m)]
            canon.add(min(rots))
    for cyc in cycles:
        if isinstance(cyc, list) and cyc:
            r = cyc[::-1]
            m = len(r)
            rr = min(tuple(r[s:] + r[:s]) for s in range(m))
            if rr in canon:
                rev_pairs += 1
    stats["cycle_length_histogram"] = dict(sorted(lens.items()))
    stats["cycle_defects"] = bad
    stats["n_identical_duplicate_cycles"] = dup_cyc
    stats["n_cycles_whose_reverse_is_also_listed"] = rev_pairs
    stats["n_distinct_cycles_up_to_rotation"] = len(canon)
    res["claim5_cycles_valid"] = all(v == 0 for v in bad.values())

    # ---- unit-distance pairs (exact, all pairs) ----
    if do_unit:
        t1 = time.time()
        unit_pairs = []
        gen_diff_pairs = 0
        for i in range(n):
            pi = P[i]
            a0, a1, a2, a3, b0, b1, b2, b3 = pi
            for j in range(i + 1, n):
                q = P[j]
                d = (q[0] - a0, q[1] - a1, q[2] - a2, q[3] - a3,
                     q[4] - b0, q[5] - b1, q[6] - b2, q[7] - b3)
                if d in Gpm:
                    gen_diff_pairs += 1
                if int_sq_norm_8(F, d) == unit:
                    unit_pairs.append((i, j))
        stats["allpairs_n_pairs_with_generator_difference"] = gen_diff_pairs
        res["allpairs_generator_difference_count_equals_edges"] = (gen_diff_pairs == len(E) == len(cay))
        U = set(unit_pairs)
        stats["n_unit_distance_pairs"] = len(U)
        nonedge_unit = sorted(U - E)
        stats["n_unit_distance_pairs_not_edges"] = len(nonedge_unit)
        stats["n_edges_not_unit_distance"] = len(E - U)
        res["all_edges_unit_distance"] = (len(E - U) == 0)
        # classify the extra unit vectors by their difference
        diffs = {}
        for (i, j) in nonedge_unit:
            d = tuple(b - a for a, b in zip(P[i], P[j]))
            key = max(d, neg(d))
            diffs[key] = diffs.get(key, 0) + 1
        stats["n_distinct_extra_unit_differences_up_to_sign"] = len(diffs)
        stats["extra_unit_differences_up_to_sign"] = sorted(([list(k), c] for k, c in diffs.items()), key=lambda x: -x[1])
        stats["unit_check_seconds"] = round(time.time() - t1, 1)
        # floating-point cross-check (does not decide anything)
        try:
            import numpy as np
            sa, sb = math.sqrt(A), math.sqrt(B)
            basisf = np.array([1.0, sa, sb, sa * sb])
            Pa = np.array(P, dtype=float)
            X = Pa[:, :4] @ basisf / D
            Y = Pa[:, 4:] @ basisf / D
            cnt = 0
            for i in range(n):
                dd = np.hypot(X[i + 1:] - X[i], Y[i + 1:] - Y[i])
                cnt += int(np.sum(np.abs(dd - 1.0) < 1e-9))
            stats["float_crosscheck_unit_pairs"] = cnt
        except Exception as ex:  # pragma: no cover
            stats["float_crosscheck_unit_pairs"] = "failed: %r" % (ex,)

    stats["seconds"] = round(time.time() - t0, 1)
    return res, stats


def main():
    args = sys.argv[1:]
    path = args[0]
    do_unit = "--no-unit" not in args
    out = None
    if "--json" in args:
        out = args[args.index("--json") + 1]
    data = load(path)
    res, stats = check_all(data, do_unit=do_unit)
    print("file:", path)
    print("field:", data.get("field"), "basis:", data.get("basis"), "D:", data.get("denominator"))
    for k, v in res.items():
        print("  CHECK %-55s %s" % (k, "PASS" if v else "FAIL"))
    for k, v in stats.items():
        if k == "extra_unit_differences_up_to_sign":
            print("  STAT %s: (%d entries, first 10) %s" % (k, len(v), v[:10]))
        else:
            print("  STAT %s: %s" % (k, v))
    allok = all(res.values())
    print("OVERALL:", "PASS" if allok else "FAIL")
    if out:
        with open(out, "w") as fh:
            json.dump({"res": res, "stats": stats, "overall": allok}, fh, indent=1, default=str)
    sys.exit(0 if allok else 1)


if __name__ == "__main__":
    main()
