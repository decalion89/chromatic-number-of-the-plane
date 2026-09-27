"""The stored three-point certificates (data/threepoint/*.npz) prove the lower bounds of
notes/local_colourings.md, section 14.

Each certificate is a dual solution of the programme of scripts/threepoint.py; scripts/threepoint_verify.py
rebuilds the blocks in interval arithmetic, proves the dual matrices positive definite by an exact rational
LDL^T and bounds alpha.  The plane is vertex-transitive, so chi >= q^2 / alpha >= q^2 / floor(bound).
Marked slow: up to four minutes each.
"""
import os, sys, json
import numpy as np
import pytest

HERE = os.path.dirname(__file__)
sys.path.append(os.path.join(HERE, '..', 'scripts'))   # appended, so that no script can shadow an installed package

CERTS = os.path.join(HERE, '..', 'data', 'threepoint')

# certificate: (q, kind, the lower bound on chi it proves)
CLAIMS = {
    'inert13.npz': (13, 'inert', 5),
    'inert29.npz': (29, 'inert', 6),
    'std37.npz': (37, 'std', 6),
    'inert37.npz': (37, 'inert', 6),
    'std41.npz': (41, 'std', 6),
    'inert41.npz': (41, 'inert', 6),
    'std43.npz': (43, 'std', 6),
    'std47.npz': (47, 'std', 6),
}


def test_every_stored_certificate_has_a_claim():
    assert sorted(f for f in os.listdir(CERTS) if f.endswith('.npz')) == sorted(CLAIMS)


@pytest.mark.slow
@pytest.mark.parametrize("name", sorted(CLAIMS))
def test_certificate(name):
    from threepoint_verify import main
    q, kind, chi = CLAIMS[name]
    meta = json.loads(str(np.load(os.path.join(CERTS, name))['meta']))
    assert (meta['q'], meta['kind']) == (q, kind)
    n = q * q
    bound, Lmax = main(os.path.join(CERTS, name), {})
    amax = int(bound.b)                   # alpha <= floor of the rigorous upper bound
    assert (chi - 1) * amax < n           # so chi(plane) >= n / alpha > chi - 1
