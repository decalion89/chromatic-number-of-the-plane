"""Theorems E and F (notes/circular_planes.md §6, papers/three-colours Section 10): the computer-assisted steps and
the finite facts at 7, rerun from the stored programs and certificates (data/number_fields/circular/twoprime/ and the
referees' subfolders indep_E/ and indep_F/)."""
import gzip, os, re, shutil, subprocess, sys

import pytest

P = os.path.join(os.path.dirname(__file__), "..", "data", "number_fields", "circular", "twoprime")


def run(cwd, *args, timeout=900, check=True):
    return subprocess.run([sys.executable, "-B", *args], cwd=cwd, capture_output=True, text=True, check=check,
                          timeout=timeout)


def stage(tmp_path, names, trees=()):
    """Copy programs and certificates into tmp_path and decompress the stored trees there."""
    for n in names:
        shutil.copy(os.path.join(P, n), tmp_path / n)
    for t in trees:
        with gzip.open(os.path.join(P, t + ".gz"), "rb") as src, open(tmp_path / t, "wb") as dst:
            shutil.copyfileobj(src, dst)
    return tmp_path


def window65(tmp_path):
    return stage(tmp_path, ["verify_tree.py", "verify_prop7.py", "prop7_certificates.txt"], ["tree_r7_25.txt"])


def test_window_65_certificates(tmp_path):
    """The window G(1,1) modulo 65 (Proposition E1 of the note): completeness tree and the 13 index vectors."""
    d = window65(tmp_path)
    out = run(d, "verify_tree.py").stdout
    assert "tree verified: 4225 cells, 15620 branch nodes, 4364 Farkas-closed branches, 13 leaves" in out
    out = run(d, "verify_prop7.py").stdout
    assert "5 main components" in out and "8 extra components with kappa = 2/7 exactly" in out


def corrupt(path, pattern, repl):
    s = path.read_text()
    t = re.sub(pattern, repl, s, count=1, flags=re.M)
    assert t != s
    path.write_text(t)


@pytest.mark.parametrize("what", ["multiplier", "dropped branch", "certificate"])
def test_window_65_rejects_corrupted_copies(tmp_path, what):
    d = window65(tmp_path)
    if what == "multiplier":
        corrupt(d / "tree_r7_25.txt", r"^(D -?\d+ \S+ \S+ )(\d+)/(\d+)", lambda m: f"{m[1]}{int(m[2]) + 1}/{m[3]}")
        prog = "verify_tree.py"
    elif what == "dropped branch":
        corrupt(d / "tree_r7_25.txt", r"^D [^\n]*\n", "")
        prog = "verify_tree.py"
    else:
        corrupt(d / "prop7_certificates.txt", r"( ; cert -?\d+ -?\d+ \S+ \S+ )(\d+)/(\d+)",
                lambda m: f"{m[1]}{int(m[2]) + 1}/{m[3]}")
        prog = "verify_prop7.py"
    r = run(d, prog, check=False)
    assert r.returncode != 0 and "verified" not in r.stdout


@pytest.mark.slow
def test_window_325_certificates(tmp_path):
    """The window G(2,1) modulo 325 (Proposition F1 of the note): 29 index vectors, 16 of them excluded at 1/4
    (about 45 seconds)."""
    d = stage(tmp_path, ["verify_window.py", "cert_K2M1.txt"], ["tree_K2M1.txt"])
    out = run(d, "verify_window.py", "2", "1").stdout
    assert "tree verified (105625 cells, 550474 branch nodes, 199920 Farkas-closed branches, 29 leaves)" in out
    assert "labels {'7': 8, 'extra': 16, 'Q': 4, 'C': 1}" in out and "EXTRA components = 1/4" in out


@pytest.mark.parametrize("prog", ["seven_local", "seven_patterns"])
def test_facts_at_seven(prog):
    """A'_7: eight elements, one mu_8-orbit, no coset of a subgroup; the digit patterns of the 7-adic characters."""
    with open(os.path.join(P, prog + ".txt")) as fh:
        assert run(P, prog + ".py").stdout == fh.read()


def test_facts_at_seven_referee():
    out = run(os.path.join(P, "indep_F"), "f49_and_digits.py").stdout
    with open(os.path.join(P, "indep_F", "f49_and_digits.txt")) as fh:
        assert out == fh.read()
    for line in ["(F1) single mu8-orbit: True ; 0 not in A'_7: True", "(F2) cosets of non-zero subgroups inside "
                 "A'_7: none", "the 8 shifts are distinct and exhaust Z/8: True", "2 consecutive digits suffice: True",
                 "7-type points consistent with the characters: True"]:
        assert line in out


@pytest.mark.parametrize("prog", ["b4_model_check", "b5_f3_check"])
def test_lemmas_at_seven_on_models(prog):
    """Lemma F6 on the finite models F_49[pi]/(pi^e), e <= 3, and Lemma F7 for f = 3."""
    d = os.path.join(P, "indep_F")
    with open(os.path.join(d, prog + ".txt")) as fh:
        assert run(d, prog + ".py").stdout == fh.read()


def test_window_65_referee_vertices():
    """The referee's vertex enumeration: 13 components at 7/25 and 2/7, 5 above 2/7."""
    d = os.path.join(P, "indep_E")
    counts = {}
    for a, b in [(7, 25), (2, 7), (201, 700), (3, 10), (1, 3)]:
        head = run(d, "ref_vertex.py", str(a), str(b)).stdout.splitlines()[0]
        counts[(a, b)] = int(head.rsplit("components:", 1)[1])
    assert counts == {(7, 25): 13, (2, 7): 13, (201, 700): 5, (3, 10): 5, (1, 3): 5}


def test_window_65_referee_clip_and_kappa(tmp_path):
    """The referee's clipping at 7/25 and the exact largest least margin of each of the 13 components."""
    for n in ["ref_clip.py", "ref_kappa.py", "ref_rot.py"]:
        shutil.copy(os.path.join(P, "indep_E", n), tmp_path / n)
    run(tmp_path, "ref_clip.py", "7", "25", "clip.txt")
    out = run(tmp_path, "ref_kappa.py", "clip.txt").stdout
    vals = re.findall(r"kappa primal = (\S+), dual = (\S+), equal: True", out)
    assert sorted(p for p, q in vals if p == q) == sorted(["2/7"] * 8 + ["1/3"] * 4 + ["1/2"])


def test_torsion_characters():
    """One prime: orders 41 and 76 beat 2/7; two primes: only types c and q (orders <= 400)."""
    d = os.path.join(P, "indep_E")
    one = run(d, "ref_torsion.py", "one", "400").stdout
    two = run(d, "ref_torsion.py", "two", "400").stdout
    assert "exact orders with an orbit of kappa > 2/7: [2, 3, 41, 76]" in one
    assert "M=41: max kappa 12/41" in one and "M=76: max kappa 11/38" in one
    assert "exact orders with an orbit of kappa > 2/7: [2, 3]" in two


def test_padic_points_killed_by_sigma():
    out = run(os.path.join(P, "indep_E"), "ref_ck.py").stdout
    assert out.count("min margin of c_k over G_N = 3/10") == 5
    assert "k=1: min margin of c_k over G_N = 3/10 ; max over m of min margin of c_k+Nm over G(k,1) = 41/250" in out
