"""Identify each component of S^theta(2,1) (output of window21_components.py) with the type points it contains, and
compare the index vectors with the stored ones in cert_K2M1.txt (29 lines) / cert_K2M1_seven.txt."""
import os
import sys
from fractions import Fraction as Fr
sys.path.insert(0, ".")
from window21_components import KEYS, FORMS

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

def load(fn):
    comps = []
    for ln in open(fn):
        if not ln.startswith("COMP"):
            continue
        head, verts = ln.split("|")
        t = head.split()
        a, b = int(t[1]), int(t[2])
        vec = tuple(int(x) for x in t[3:])
        poly = [tuple(Fr(u) for u in v.split(",")) for v in verts.split()]
        comps.append((a, b, vec, poly))
    return comps

def floors(c):
    out = []
    for k in KEYS:
        g0, g1 = FORMS[k]
        v = (g0 * c[0] + g1 * c[1]) / 325
        out.append(v.numerator // v.denominator)
    return tuple(out)

def least_margin(c):
    m = None
    for k in KEYS:
        g0, g1 = FORMS[k]
        v = (g0 * c[0] + g1 * c[1]) / 325
        fl = v.numerator // v.denominator
        mm = min(v - fl, fl + 1 - v)
        m = mm if m is None else min(m, mm)
    return m

# type points modulo 325 Z[i], represented in [0,325)^2
A7 = {(2, 3), (2, 4), (3, 2), (3, 5), (4, 2), (4, 5), (5, 3), (5, 4)}
TP = [("c", (Fr(325, 2), Fr(325, 2)))]
TP += [(f"q:({al}+{be}i)/3", (Fr(325 * al, 3), Fr(325 * be, 3))) for al in (1, 2) for be in (1, 2)]
TP += [(f"7:({a}+{b}i)/7", (Fr(325 * a, 7), Fr(325 * b, 7))) for a in range(7) for b in range(7) if ((3 * a) % 7, (3 * b) % 7) in A7]
assert len(TP) == 13
tp_vec = {name: floors(pt) for name, pt in TP}
for name, pt in TP:
    print(f"type point {name}: least margin {least_margin(pt)}")

theta_file = sys.argv[1]
comps = load(theta_file)
print(f"{theta_file}: {len(comps)} components")
stored = {}
for ln in open(f"{D}/cert_K2M1.txt"):
    if ln.startswith("component"):
        head, rest = ln.split("|", 1)
        stored[tuple(int(x) for x in head.split()[2:])] = rest.strip().split()[0]
matched = {}
for a, b, vec, poly in comps:
    names = [n for n, v in tp_vec.items() if v == vec]
    kind = "point" if len(poly) == 1 else ("segment" if len(poly) == 2 else f"{len(poly)}-gon")
    st = stored.get(vec, "NOT IN cert_K2M1.txt")
    print(f"  cell ({a},{b}) {kind:8s} index vector = type point {names if names else 'NONE'}; stored label: {st}")
    for n in names:
        matched[n] = matched.get(n, 0) + 1
print("type points matched:", len(matched), "of 13;", "every component matched:", all(any(v == vec for v in tp_vec.values()) for _, _, vec, _ in comps))
print("all component vectors among the 29 stored:", all(vec in stored for _, _, vec, _ in comps))
