"""Four colours for every d = 11 (mod 12) (notes/four_colours_11_mod_12.md).

Checked here: the exact facts behind the structure lemma for S_{5^k} (structure_lemma.py) and behind its crude form,
the one proved in lean/FourColours.lean (crude_check.py); that the configurations of
Theorems 1a and 1b consist of unit vectors and satisfy the vector relations used in the proofs; the 2-adic and 3-adic
valuations of Steps 2 and 3 of the proof of Theorem 1b for every d = 11 (mod 24) below 20000; the large-N limit
criterion (limit.py) on small cases; and the stored exact certificates (check_w.py), which must describe exactly the
configurations of family23.py / family11.py."""
import glob
import gzip
import json
import os
import subprocess
import sys
from fractions import Fraction as Fr

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WDIR = os.path.join(ROOT, "data", "quadratic_planes", "winding")
FDIR = os.path.join(WDIR, "family")
sys.path.insert(0, WDIR)
sys.path.insert(0, FDIR)
from check_w import check  # noqa: E402
import family11  # noqa: E402
import family23  # noqa: E402


def v(p, x):
    s = 0
    while x % p == 0:
        s += 1
        x //= p
    return s


def test_structure_lemma_exact_facts():
    out = subprocess.run([sys.executable, os.path.join(FDIR, "structure_lemma.py"), "4", "2"],
                         capture_output=True, text=True, timeout=600)
    assert out.returncode == 0, out.stdout + out.stderr
    assert "ALL CHECKS PASSED" in out.stdout


def test_crude_structure_lemma_facts():
    """The exact facts behind the crude form of the structure lemma used by lean/FourColours.lean."""
    out = subprocess.run([sys.executable, os.path.join(FDIR, "crude_check.py")], capture_output=True, text=True,
                         check=True).stdout
    assert "ALL PLAN CHECKS PASSED" in out


@pytest.mark.parametrize("d,k", [(23, 2), (47, 2), (119, 3)])
def test_family23_configuration(d, k):
    L, U = family23.config(d, k)
    assert len(U) == 6 * (2 * k + 1)
    for a, b, c, e in U:
        assert a * a + d * b * b + c * c + d * e * e == L * L and a * b + c * e == 0


@pytest.mark.parametrize("d,k", [(11, 2), (35, 2), (59, 2), (107, 3)])
def test_family11_configuration(d, k):
    n = family11.default_n(d)
    assert n % 4 == 3 and (n - 1) % 3 ** (v(3, d + 1) + 1) == 0
    L, U = family11.config(d, k, n)
    assert len(U) == 8 * (2 * k + 1)
    for a, b, c, e in U:
        assert a * a + d * b * b + c * c + d * e * e == L * L and a * b + c * e == 0


def test_vector_relations():
    """a(u1 + conj u1) + b = 0 and (n^2 + d) u_n = (n - 1)(n + d) + n(1 + d) u1, as (X, Y) with v = X + i Y sqrt d."""
    for d in range(11, 2000, 12):
        for n in (1, 3, 19, 55, 163):
            def u(m, sg=1):
                t = Fr(m * m + d)
                return (Fr(m * m - d) / t, Fr(2 * m * sg) / t)
            a, b = Fr(d + 1, 4), Fr(d - 1, 2)
            u1, u1c, un = u(1), u(1, -1), u(n)
            assert a * (u1[0] + u1c[0]) + b == 0 and a * (u1[1] + u1c[1]) == 0
            assert (n * n + d) * un[0] == (n - 1) * (n + d) + n * (1 + d) * u1[0]
            assert (n * n + d) * un[1] == n * (1 + d) * u1[1]


def test_valuations_theorem_1b():
    for d in range(11, 20000, 24):
        s = v(3, d + 1)
        assert s >= 1 and v(3, (d + 1) // 4) == s and (d - 1) // 2 % 3 != 0
        n = family11.default_n(d)
        e = v(3, n - 1)
        assert e >= s + 1 and n % 3 != 0
        # Step 2 (at 2)
        assert v(2, n * n + d) == 2 and v(2, n * (1 + d)) == 2 and v(2, (n - 1) * (n + d)) == 2
        # Step 3 (at 3)
        assert v(3, n * n + d) == s and v(3, n * (1 + d)) == s and v(3, (n - 1) * (n + d)) == e + s
        assert e + s - 1 >= 2 * s >= 2


def test_limit_criterion_small_cases():
    from limit import solutions, ONE, ntype, tconj, utype
    assert not solutions([ONE, ntype(23, 1), tconj(ntype(23, 1))])
    assert solutions([ONE, ntype(11, 1), tconj(ntype(11, 1))])
    for d in (11, 35, 59, 107):
        n = family11.default_n(d)
        assert not solutions([ONE, utype(d, 1, 1), utype(d, 1, -1), utype(d, n, 1)])


CERTS = sorted(glob.glob(os.path.join(FDIR, "cert_*.json.gz")))
# CI checks the small certificates (seconds each); the others take minutes and are checked locally.
FAST = {"cert_f23_23_k2.json.gz", "cert_f23_71_k2.json.gz", "cert_f23_95_k2.json.gz", "cert_f23_119_k3.json.gz",
        "cert_f11_11_k2.json.gz"}


def _check_certificate(path):
    msg = check(path)
    assert msg.startswith("VERIFIED")
    C = json.load(gzip.open(path))
    d, D = C["d"], C["D"]
    name = os.path.basename(path)
    k = int(name.split("_k")[1].split(".")[0])
    if name.startswith("cert_f23_"):
        L, U = family23.config(d, k)
    else:
        L, U = family11.config(d, k, family11.default_n(d))
    assert L == D and sorted(map(tuple, U)) == sorted(map(tuple, C["units"]))


def test_fast_certificates_present():
    assert FAST <= {os.path.basename(p) for p in CERTS}


@pytest.mark.parametrize("path", [p for p in CERTS if os.path.basename(p) in FAST],
                         ids=lambda p: os.path.basename(p))
def test_certificates(path):
    _check_certificate(path)


@pytest.mark.slow
@pytest.mark.parametrize("path", [p for p in CERTS if os.path.basename(p) not in FAST],
                         ids=lambda p: os.path.basename(p))
def test_certificates_slow(path):
    _check_certificate(path)


def test_sharpness_decision():
    """The referee's exact decision for d = 23 (mod 24) (family/fc_case23.py, given the shape of S_N): the bound
    d < 21N/5 is sharp for N = 5, ..., 3125 and not for N = 15625 (22 exceptions, 65 639 <= d <= 66 143)."""
    out = subprocess.run([sys.executable, "-B", "fc_case23.py"], cwd=FDIR, capture_output=True, text=True,
                         check=True, timeout=600).stdout
    assert out == open(os.path.join(FDIR, "fc_case23.out")).read()
    assert "N=3125: 21N/5 = 13125; first feasible d (d = 23 mod 24) in the window: 13127" in out
    assert "N=15625: 21N/5 = 65625; first feasible d (d = 23 mod 24) in the window: 66167" in out
    assert "[65639, 65663, 65687, 65711, 65735]... (22 values)" in out


def _indep_lines():
    with open(os.path.join(WDIR, "check_w_indep.out")) as fh:
        return {line.split(" ", 1)[0]: line for line in fh}


@pytest.mark.parametrize("name", ["family/cert_f23_23_k2.json.gz", "family/cert_f11_11_k2.json.gz",
                                  "cert_11_30.json.gz"])
def test_second_checker(name):
    """check_w_indep.py, written separately from check_w.py, accepts the certificate with the stored output."""
    out = subprocess.run([sys.executable, "-B", "check_w_indep.py", name], cwd=WDIR, capture_output=True,
                         text=True, check=True, timeout=600).stdout
    assert out == _indep_lines()[name]
    assert "'rankU': 4" in out and "'saturated': True" in out


@pytest.mark.slow
def test_second_checker_all():
    """check_w_indep.py on all 26 stored certificates reproduces check_w_indep.out (about a minute)."""
    names = sorted(_indep_lines())
    out = subprocess.run([sys.executable, "-B", "check_w_indep.py", *names], cwd=WDIR, capture_output=True,
                         text=True, check=True, timeout=3600).stdout
    assert sorted(out.splitlines(keepends=True)) == sorted(_indep_lines().values()) and len(names) == 26
