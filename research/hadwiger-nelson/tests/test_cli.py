"""The command line of hn/cli.py: `verify` accepts the stored certificates it handles and rejects a
tampered one, and `demo` rebuilds the Moser spindle certificate."""
import json
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _hn(*args, cwd=ROOT):
    return subprocess.run([sys.executable, "-m", "hn.cli", *args], cwd=cwd, capture_output=True, text=True,
                          timeout=600, env=dict(os.environ, PYTHONPATH=ROOT))


def test_verify_accepts_the_moser_spindle_certificates():
    for name in ("moser_spindle_no3coloring.json", "moser_spindle_4coloring.json"):
        r = _hn("verify", os.path.join("certificates", name))
        assert r.returncode == 0 and r.stdout.startswith("VERIFIED"), r.stdout + r.stderr


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
    cert = tmp_path / "moser_spindle_no3coloring.json"
    r = _hn("verify", str(cert))
    assert r.returncode == 0 and r.stdout.startswith("VERIFIED"), r.stdout + r.stderr
