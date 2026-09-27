"""scripts/g13/g13_audit.py: every clause of the five formulas of alpha(G_13) <= 36 (E37_A6, E37_A7, E37_A9, E37_A11
and E37_B, as scripts/g13/enum_cert.py writes them) holds in the intended models, checked from the DIMACS text alone;
and the audit rejects a formula with one wrong clause, for the right reason. The formulas are written and parsed
once; each wrong formula is a copy of a list of clauses with one change."""
import os
import subprocess
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
G13 = os.path.join(ROOT, "scripts", "g13")
sys.path.insert(0, G13)
import g13                  # noqa: E402
import g13_audit as au      # noqa: E402
import g13cnf as gc         # noqa: E402

ARGS = {f"E37_A{c}": ["--L", "0", "--rosette", str(c)] for c in (6, 7, 9, 11)}
ARGS["E37_B"] = ["--norosette"]
FOUND_A = {"clauses": 7803, "edge clauses": 1183, "domination clauses": 169, "unit clauses": 15, "totalizer nodes": 168}


@pytest.fixture(scope="module")
def texts(tmp_path_factory):
    """the five formulas, as enum_cert.py writes them"""
    d = tmp_path_factory.mktemp("g13_audit")
    out = {}
    for name, args in ARGS.items():
        path = str(d / f"{name}.cnf")
        subprocess.run([sys.executable, os.path.join(G13, "enum_cert.py"), path, "37"] + args, check=True,
                       capture_output=True)
        with open(path) as fh:
            out[name] = fh.read()
    return out


@pytest.fixture(scope="module")
def parsed(texts):
    """name -> (number of variables, clauses)"""
    return {name: au.parse(text) for name, text in texts.items()}


def test_the_audit_rebuilds_g13_on_its_own():
    edges = {frozenset((u, v)) for u in au.VERTICES for v in au.NEIGHBOURS[u]}
    assert len(au.VERTICES) == 169 and all(len(n) == 14 for n in au.NEIGHBOURS)
    assert edges == {frozenset(e) for e in g13.EDGES} and len(edges) == 1183
    assert {c: sorted(p) for c, p in au.CIRCLES.items()} == {c: sorted(g13.INDEX[z] for z in g13.circle(c))
                                                             for c in range(1, 13)}
    assert all(len(p) == 14 for p in au.CIRCLES.values())
    assert au.INDEPENDENT == [1, 6, 7, 9, 11] and au.ROSETTE_CIRCLES == [6, 7, 9, 11]
    assert len(au.MAPS) == 4732 and set(au.MAPS) == {tuple(m) for m in g13.automorphisms()}
    # every map is a linear map followed by a translation; these are automorphisms, and they send the circle C_c
    # around a point p onto the circle C_c around the image of p
    around = {pc: s for s, pc in au.AROUND.items()}
    for g in au.LINEAR + au.TRANSLATIONS:
        assert all(frozenset((g[u], g[v])) in edges for u, v in map(tuple, edges))
        assert all(frozenset(g[x] for x in s) == around[g[p], c] for (p, c), s in around.items())


def test_the_audit_reads_nothing_but_dimacs():
    assert au.parse("p cnf 3 2\n1 -2 0\n3 0\n") == (3, [(1, -2), (3,)])
    for bad in ("p cnf 3 3\n1 -2 0\n3 0\n", "p cnf 2 2\n1 -2 0\n3 0\n", "p cnf 3 2\n1 -2\n3 0\n",
                "p cnf 3 2\n1 0 -2 0\n3 0\n", "c a comment\np cnf 3 2\n1 -2 0\n3 0\n"):
        with pytest.raises(au.AuditError):
            au.parse(bad)


def test_every_clause_of_the_five_formulas_holds_in_the_intended_models(texts):
    for c in (6, 7, 9, 11):
        assert au.audit_formula_a(c, texts[f"E37_A{c}"]) == FOUND_A
    order = gc.lex_order()[:25]
    moved = sum(m[v] != v for m in g13.automorphisms() for v in order)    # 4731 * 25 - 25 * 27: each point is
    assert moved == 117600                                                # fixed by 27 maps other than the identity
    assert au.audit_formula_b(texts["E37_B"]) == {
        "clauses": 351802, "edge clauses": 1183, "domination clauses": 169, "no-rosette clauses": 676,
        "totalizer nodes": 168, "chains": 4731, "comparisons": moved, "maps": 4731, "prefix": tuple(order)}


# the clauses of formula E37_B, in the order enum_cert.py writes them: the edges, the domination clauses, the
# totalizer (its root unit last), the chains, the no-rosette clauses

def first(clauses, kind):
    return next(k for k, c in enumerate(clauses) if kind(c))


def replaced(clauses, k, clause):
    return clauses[:k] + [clause] + clauses[k + 1:]


def root_of(clauses):
    """the variable of the unit clause o_37 of the root of the totalizer: the last counter"""
    return next(c[0] for c in clauses if len(c) == 1)


def first_chain(clauses):
    """(start, end, e_0): clauses[start:end] is the first chain, x_v0 | -x_w0, x_v0 | e_0, -x_w0 | e_0, then
    -e_0 | x_v1 | -x_w1, -e_0 | x_v1 | e_1, -e_0 | -x_w1 | e_1, and so on; e_0 is the variable after the counters"""
    root = root_of(clauses)
    start = clauses.index((root,)) + 1
    end = next(k for k in range(start + 1, len(clauses)) if max(map(abs, clauses[k])) <= 169)
    return start, end, root + 1


def rosette_of(p, c):
    """the vertices of the point p with its circle N = c"""
    return [g13.INDEX[z] for z in [p] + [g13.add(p, q) for q in g13.circle(c)]]


def non_edge(clauses):
    k = first(clauses, lambda c: len(c) == 2 and c[0] < 0 and c[1] < 0)
    u = -clauses[k][0] - 1
    w = next(v for v in range(169) if v != u and v not in g13.ADJ[u])
    return replaced(clauses, k, (-(u + 1), -(w + 1)))


def domination_swapped(clauses):
    """a neighbour in a domination clause swapped for a vertex that is not one"""
    k = first(clauses, lambda c: len(c) == 15 and c[0] > 0)
    w = next(v for v in range(1, 170) if v not in clauses[k])
    return replaced(clauses, k, clauses[k][:-1] + (w,))


def rosette_point_moved(clauses):
    k = first(clauses, lambda c: len(c) == 15 and c[0] < 0)
    w = next(v for v in range(1, 170) if -v not in clauses[k])
    return replaced(clauses, k, clauses[k][:5] + (-w,) + clauses[k][6:])


def rosette_on_circle_2(clauses):
    """a clause against a point with its circle N = 2, which has a unit distance inside"""
    k = first(clauses, lambda c: len(c) == 15 and c[0] < 0)
    return replaced(clauses, k, tuple(-(v + 1) for v in rosette_of((2, 3), 2)))


def rosette_on_circle_1(clauses):
    """a clause against a point with its circle N = 1: independent, but the neighbours of the point"""
    k = first(clauses, lambda c: len(c) == 15 and c[0] < 0)
    return replaced(clauses, k, tuple(-(v + 1) for v in rosette_of((2, 3), 1)))


def rosette_clause_dropped(clauses):
    """one clause against a rosette left out: the maps no longer send the ones that are left onto themselves"""
    k = first(clauses, lambda c: len(c) == 15 and c[0] < 0)
    return clauses[:k] + clauses[k + 1:]


def totalizer_clause_dropped(clauses):
    root = root_of(clauses)
    ks = [k for k, c in enumerate(clauses) if len(c) == 3 and 169 < max(map(abs, c)) <= root]
    k = ks[len(ks) // 2]
    return clauses[:k] + clauses[k + 1:]


def totalizer_leaf_changed(clauses):
    """a leaf in a clause of the totalizer changed to the next vertex"""
    root = root_of(clauses)
    k = first(clauses, lambda c: 169 < max(map(abs, c)) <= root and any(0 < l <= 169 for l in c))
    return replaced(clauses, k, tuple(l % 169 + 1 if 0 < l <= 169 else l for l in clauses[k]))


def threshold_36(clauses):
    """the unit clause o_36 of the root instead of o_37"""
    k = clauses.index((root_of(clauses),))
    return replaced(clauses, k, (clauses[k][0] - 1,))


def chain_image_changed(clauses):
    """the vertex that the first chain compares with its second position, changed in both clauses that have it"""
    start, _, e0 = first_chain(clauses)
    (_, x, y), (_, _, e1) = clauses[start + 3], clauses[start + 4]      # -e_0 | x_v1 | -x_w1, -e_0 | x_v1 | e_1
    w = next(u for u in (-y % 169 + 1, -y % 169 + 2) if u != x)
    out = list(clauses)
    out[start + 3], out[start + 5] = (-e0, x, -w), (-e0, -w, e1)
    return out


def chain_skips_a_position(clauses):
    """the second pair of the first chain left out: its three clauses dropped, the chain continued from e_0"""
    start, _, e0 = first_chain(clauses)
    return clauses[:start + 3] + [(-e0,) + c[1:] for c in clauses[start + 6:start + 9]] + clauses[start + 9:]


def chain_order_swapped(clauses):
    """the second and the third pair of the first chain compared in the other order"""
    start, _, _ = first_chain(clauses)
    out = list(clauses)
    one, two = clauses[start + 3], clauses[start + 6]                   # -e_0 | x_v1 | -x_w1, -e_1 | x_v2 | -x_w2
    for k, (_, x, y) in ((start + 3, two), (start + 6, one)):
        e = clauses[k][0]
        out[k], out[k + 1], out[k + 2] = (e, x, y), (e, x, clauses[k + 1][2]), (e, y, clauses[k + 2][2])
    return out


def chain_clause_dropped(clauses):
    start, _, _ = first_chain(clauses)
    return clauses[:start + 7] + clauses[start + 8:]


def stray_comparison(clauses):
    """x_0 | -x_1 once more: a chain of one pair, whose maps (0 -> 1) all move other positions of the prefix"""
    return clauses + [(1, -2)]


def stray_unit(clauses):
    return clauses + [(5,)]


WRONG_B = [(non_edge, "not adjacent"), (domination_swapped, "not the closed neighbourhood"),
           (rosette_point_moved, "not a point with its whole circle"), (rosette_on_circle_2, "not independent"),
           (rosette_on_circle_1, "not independent"), (rosette_clause_dropped, "no-rosette clauses onto themselves"),
           (totalizer_clause_dropped, "the clauses of counter"), (totalizer_leaf_changed, "node"),
           (threshold_36, "not o_37 of the root"), (chain_image_changed, "no map of the group"),
           (chain_skips_a_position, "skips position"), (chain_order_swapped, "one common order"),
           (chain_clause_dropped, "has 1 clauses that define it"), (stray_comparison, "skips position"),
           (stray_unit, "a unit clause in formula E37_B")]


@pytest.mark.parametrize("mutation, reason", WRONG_B, ids=[m.__name__ for m, _ in WRONG_B])
def test_the_audit_rejects_a_wrong_formula_b(parsed, mutation, reason):
    nv, clauses = parsed["E37_B"]
    wrong = mutation(clauses)
    assert wrong != clauses
    with pytest.raises(au.AuditError, match=reason):
        au.audit_clauses_b(nv, wrong)


def test_the_audit_needs_neither_every_chain_nor_every_circle(parsed):
    """formula E37_B without its first chain and without the clauses against the rosettes of the circle N = 11: a
    weaker formula, still symmetric"""
    nv, clauses = parsed["E37_B"]
    start, end, _ = first_chain(clauses)
    circle11 = {frozenset(-(v + 1) for v in rosette_of(p, 11)) for p in g13.POINTS}
    weaker = [c for c in clauses[:start] + clauses[end:] if frozenset(c) not in circle11]
    found = au.audit_clauses_b(nv, weaker)
    assert (found["chains"], found["maps"], found["no-rosette clauses"]) == (4730, 4730, 3 * 169)


def with_units(clauses, c):
    """the clauses of a formula E37_A_c other than its unit clauses on vertices, then x_0 and x_q for q in C_c"""
    rest = [cl for cl in clauses if not (len(cl) == 1 and cl[0] <= 169)]
    return rest + [(v + 1,) for v in rosette_of((0, 0), c)]


def test_the_audit_rejects_a_formula_a_with_the_wrong_units(parsed):
    nv, clauses = parsed["E37_A6"]
    assert au.audit_clauses_a(7, nv, with_units(clauses, 7)) == FOUND_A      # formula E37_A7 in another order
    with pytest.raises(au.AuditError, match="not on 0 and the points of the circle N = 6"):
        au.audit_clauses_a(6, nv, with_units(clauses, 7))
    moved = with_units(clauses, 6)
    moved[-1] = (next(v for v in range(1, 170) if (v,) not in moved),)
    with pytest.raises(au.AuditError, match="not on 0 and the points of the circle N = 6"):
        au.audit_clauses_a(6, nv, moved)
    for c in (1, 2):            # {0} + C_c is not independent: 0 is adjacent to C_1, and C_2 has a unit distance
        with pytest.raises(au.AuditError, match="not a rosette"):
            au.audit_clauses_a(c, nv, with_units(clauses, c))
    with pytest.raises(au.AuditError, match="a first comparison clause in formula E37_A_6"):
        au.audit_clauses_a(6, nv, clauses + [(1, -2)])


@pytest.mark.parametrize("leaves, threshold, reason", [(range(169, 0, -1), 37, None),
                                                        (range(1, 169), 37, "every leaf"),
                                                        (range(1, 170), 36, "not o_1, ..., o_37")])
def test_the_totalizer_may_have_any_tree_but_must_count_all_169_vertices_to_37(parsed, leaves, threshold, reason):
    """formula E37_A6 with the totalizer of g13cnf.py over other leaves, or with another threshold"""
    nv, clauses = parsed["E37_A6"]
    cnf = gc.CNF()
    cnf.top = 169
    gc.totalizer_atleast(cnf, list(leaves), threshold)
    other = [c for c in clauses if max(map(abs, c)) <= 169] + [tuple(c) for c in cnf.clauses]
    if reason is None:
        assert au.audit_clauses_a(6, cnf.top, other) == FOUND_A
    else:
        with pytest.raises(au.AuditError, match=reason):
            au.audit_clauses_a(6, cnf.top, other)
