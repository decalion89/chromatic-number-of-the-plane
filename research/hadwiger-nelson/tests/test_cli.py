"""The command line of hn/cli.py: `verify` accepts the stored certificates it handles and rejects a
tampered one, never calls an unchecked proof verified, and `demo` rebuilds the Moser spindle
certificate outside the repository."""
import json
import os
import shutil
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DRAT = shutil.which("drat-trim")


def _hn(*args, cwd=ROOT, **env):
    return subprocess.run([sys.executable, "-m", "hn.cli", *args], cwd=cwd, capture_output=True, text=True,
                          timeout=600, env=dict(os.environ, PYTHONPATH=ROOT, **env))


def _check_uncolourable(r):
    """With drat-trim on PATH the proof is checked. Without it the solver's answer is reported as
    such, with exit status 2, and not as verified."""
    if DRAT:
        assert r.returncode == 0 and r.stdout.startswith("VERIFIED: no proper"), r.stdout + r.stderr
    else:
        assert r.returncode == 2 and r.stdout.startswith("UNSAT for k="), r.stdout + r.stderr
        assert "VERIFIED" not in r.stdout, r.stdout


def test_verify_accepts_the_moser_spindle_certificates():
    r = _hn("verify", os.path.join("certificates", "moser_spindle_4coloring.json"))
    assert r.returncode == 0 and r.stdout.startswith("VERIFIED: ["), r.stdout + r.stderr
    assert r.stdout.count("VERIFIED") == 1, r.stdout
    _check_uncolourable(_hn("verify", os.path.join("certificates", "moser_spindle_no3coloring.json")))


def test_verify_rejects_a_tampered_colouring(tmp_path):
    d = json.load(open(os.path.join(ROOT, "certificates", "moser_spindle_4coloring.json")))
    d["coloring"][1] = d["coloring"][0]      # vertices 0 and 1 are at distance 1 in the spindle
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps(d))
    r = _hn("verify", str(bad))
    assert r.returncode == 1 and r.stdout.startswith("REJECTED"), r.stdout + r.stderr


def test_demo_writes_a_certificate_that_verifies(tmp_path):
    r = _hn("demo", "--out", str(tmp_path))
    assert r.returncode == 0, r.stdout + r.stderr
    _check_uncolourable(_hn("verify", str(tmp_path / "moser_spindle_no3coloring.json")))


def test_demo_writes_to_hn_out_by_default(tmp_path):
    r = _hn("demo", HN_OUT=str(tmp_path))
    assert r.returncode == 0, r.stdout + r.stderr
    assert (tmp_path / "moser_spindle_no3coloring.json").exists()
