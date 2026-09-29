"""The 852-vertex Moser-spindle-free 5-chromatic unit-distance graph (notes/flat852.md): exact checks, no solver."""
import copy
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import verify_flat852 as vf  # noqa: E402

D = os.path.join(ROOT, "data", "flat852")


def load(name):
    with open(os.path.join(D, name)) as fh:
        return json.load(fh)


def colours(name):
    c = load(name)
    return c.get("colours", c.get("colouring")) if isinstance(c, dict) else c


def test_cyclotomic_polynomial():
    assert vf.PHI == [1, -1, 0, 1, -1, 0, 1, 0, -1, 1, 0, -1, 1]      # Phi_21, lowest degree first
    z = [0, 1] + [0] * 10
    p = [1] + [0] * 11
    for _ in range(21):
        p = vf.mul(p, z)
    assert p == [1] + [0] * 11                                           # z^21 = 1
    assert vf.mul(z, vf.conj(z)) == [1] + [0] * 11                       # z conj(z) = 1


def test_directions():
    g = load("core852a.json")
    U = [tuple(u) for u in g["U"]]
    assert vf.check_directions(U)
    omega7 = (4, -3, 1, 6, 1, -5, 1, 2, -7, 4, 1, -2)                   # 7 omega in the power basis of zeta21
    assert omega7 in set(U) and tuple(vf.conj(list(omega7))) in set(U)
    # 7^12 omega^12 - 13 * 7^6 omega^6 * 7^6 + 7 * 7^12 = 0, i.e. 7 w^12 - 13 w^6 + 7 = 0 for w = omega
    w6 = [1] + [0] * 11
    for _ in range(6):
        w6 = vf.mul(w6, list(omega7))
    w12 = vf.mul(w6, w6)
    lhs = [7 * a - 13 * 7 ** 6 * b + 7 * 7 ** 12 * (1 if i == 0 else 0) for i, (a, b) in enumerate(zip(w12, w6))]
    assert lhs == [0] * 12


def test_core852_is_a_unit_distance_graph():
    ok, msg = vf.check_graph(load("core852a.json"))
    assert ok, msg
    g = load("core852a.json")
    assert len(g["vertices"]) == 852 and len(g["edges"]) == 4487
    assert {j for _, _, j in g["edges"]} == set(range(126))              # every one of the 126 directions is used


def test_core852_five_colouring():
    ok, msg = vf.check_colouring(load("core852a.json"), colours("core852.5colouring.json"))
    assert ok and "5 colours" in msg, msg


def test_stored_cnf_is_the_formula_of_the_graph():
    ok, msg = vf.check_stored_cnf(load("core852a.json"), os.path.join(D, "core852a.cnf"))
    assert ok, msg


def test_certificate_logs():
    with open(os.path.join(D, "logs", "core852a.drat-trim.log")) as fh:
        assert "s VERIFIED" in fh.read()
    with open(os.path.join(D, "independent-check", "mine852.drat-trim.log")) as fh:
        assert "s VERIFIED" in fh.read()


def test_a_wrong_point_is_caught():
    g = copy.deepcopy(load("core852a.json"))
    g["vertices"][5]["exact7"][0] += 1
    assert not vf.check_graph(g)[0]
    g = load("core852a.json")
    bad = list(colours("core852.5colouring.json"))
    a, b, _ = g["edges"][0]
    bad[b] = bad[a]
    assert not vf.check_colouring(g, bad)[0]


def test_g1023():
    g = load("g1023.json")
    ok, msg = vf.check_graph(g)
    assert ok, msg
    assert vf.check_colouring(g, colours("g1023.5colouring.json"))[0]
    assert vf.check_stored_cnf(g, os.path.join(D, "g1023.cnf"))[0]
