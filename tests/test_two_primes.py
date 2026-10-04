"""Theorems E and F (notes/circular_planes.md §6, papers/three-colours Section 10): the computer-assisted steps and
the finite facts at 7, rerun from the stored programs and certificates (data/number_fields/circular/twoprime/ and the
referees' subfolders indep_E/, indep_F/, indep_S10/, indep_FW/ and indep_FW2/)."""
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


def test_window_325_seven_certificates(tmp_path):
    """The 8 index vectors of family 7 in the window G(2,1) have largest least margin exactly 2/7, so for theta > 2/7
    that window leaves only the types c and q (Theorem E again, and Theorem F from Proposition F1 alone)."""
    d = stage(tmp_path, ["verify_seven_K2M1.py", "cert_K2M1.txt", "cert_K2M1_seven.txt"])
    assert "largest least margin exactly 2/7" in run(d, "verify_seven_K2M1.py").stdout
    corrupt(d / "cert_K2M1_seven.txt", r"(\| kappa = 2/7 ; cert -?\d+ -?\d+ \S+ \S+ )(\d+)/(\d+)",
            lambda m: f"{m[1]}{int(m[2]) + 1}/{m[3]}")
    r = run(d, "verify_seven_K2M1.py", check=False)
    assert r.returncode != 0 and "verified" not in r.stdout


def test_section10_referee(tmp_path):
    """The referee of Section 10 of the paper: the facts at 7 and the type points modulo 325, the points c_k, torsion,
    an archimedean point of Q(i), the family-7 vectors of G(2,1) restricted to G(1,1), and its own clipping of the
    window G(1,1) at 2/7 (the same 13 index vectors as the certificates)."""
    d = os.path.join(P, "indep_S10")
    for n in os.listdir(d):
        shutil.copy(os.path.join(d, n), tmp_path / n)
    (tmp_path / "run").mkdir()
    for n in ["cert_K2M1.txt", "prop7_certificates.txt"]:
        shutil.copy(os.path.join(P, n), tmp_path / "run" / n)
    for prog in ["check_seven", "check_ck", "check_qi_char", "check_torsion"]:
        assert run(tmp_path, prog + ".py").stdout == (tmp_path / (prog + ".txt")).read_text()
    out = run(tmp_path, "check_remark1.py").stdout.splitlines()
    assert {x: out.count(x) for x in out} == {"MAIN -> MAIN": 5, "SEVEN -> EXTRA": 8,
                                               "EXTRA -> not among the 13 (1,1)-vectors": 16}
    out = run(tmp_path, "clip_window.py", "1", "1", "2/7", "run/prop7_certificates.txt").stdout
    assert "(K,M)=(1,1) N=65 r0=2/7: 13 index vectors with nonempty polygons" in out
    assert "same set of index vectors as the certificate: True" in out


def test_tight_walks(tmp_path):
    """Corollary on finite witnesses: the eight certificates of value 2/7 give positive integer relations among tight
    rotations (closed walks of length 140 or 42), type q has relations 3, 4, 5 in G_5, and the Moser spindle has
    chi = 4 and maps to K_{7/2}; the program reproduces its stored output."""
    d = stage(tmp_path, ["tight_walks.py", "cert_K2M1_seven.txt"])
    out = run(d, "tight_walks.py").stdout
    assert out == open(os.path.join(P, "tight_walks.txt")).read()
    assert out.count("multiplicities [63, 52, 25], length 140") == 4
    assert out.count("multiplicities [13, 14, 15], length 42") == 4
    assert out.count("multiplicities (3, 4, 5)") == 8
    assert "proper 3-colourings: 0; a homomorphism to K_7/2: (0, 2, 4, 6, 3, 5, 1)" in out


def test_finite_ball(tmp_path):
    """Remark after the corollary on finite witnesses: the ball of radius 2 in the Cayley graph of the 140 vectors of
    Theorem C is bipartite (19 600 edges), and the unit-distance graph it induces (216 more edges) has a 5-cycle, no
    triangle and a homomorphism to K_{5/2}; the program reproduces its stored output."""
    d = stage(tmp_path, ["finite_ball.py"])
    out = run(d, "finite_ball.py").stdout
    assert out == open(os.path.join(P, "finite_ball.txt")).read()
    assert "9941 points; unit-distance pairs 19816; Cayley edges 19600; others 216" in out
    assert "Cayley ball bipartite: True" in out and "no triangle" in out
    assert "circular chromatic number 5/2" in out


def test_finite_witness_referee(tmp_path):
    """The referee of the corollary on finite witnesses: its checks of the certificates, of type q, of the tight
    relations of the 7-adic colourings, of the Moser spindle and of the tight-cycle and winding lemmas reproduce."""
    src = os.path.join(P, "indep_FW")
    d = tmp_path / "indep_FW"
    shutil.copytree(src, d)
    for n in ["cert_K2M1.txt", "cert_K2M1_seven.txt"]:
        shutil.copy(os.path.join(P, n), tmp_path / n)
    for prog in ["check_seven_certs", "check_type_q", "seven_adic_tight", "moser", "test_lemmas_AB"]:
        assert run(d, prog + ".py").stdout == (d / (prog + ".out")).read_text()


def test_finite_connection_sets_referee(tmp_path):
    """The referee of Lemma 22 (finite connection sets): the radii of Lemma 9 (sup r_theta = sqrt10/14 < 0.23 on
    (2/7, 1/3]) and the 7-adic characters on the 140 vectors of Theorem C, whose tight sets carry positive relations."""
    d = tmp_path / "indep_FW2"
    shutil.copytree(os.path.join(P, "indep_FW2"), d)
    out = run(d, "check_radii.py").stdout
    assert "max |y|^2 over X = 2" in out and "sup s*sqrt10/3 = 0.2258769757" in out
    out = run(d, "thmC_tight7.py").stdout
    assert out.count("all values in {2..5}/7: True;  tight-set sizes: [35]  positive relation found & verified exactly for all: True") == 2

