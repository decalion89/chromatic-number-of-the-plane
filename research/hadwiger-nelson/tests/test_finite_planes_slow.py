"""chi(F_q^2) >= 5 for q = 29, 31, 41, 43 (notes/local_colourings.md, section 12): no proper 4-colouring.

These planes have no unit triangle (q = 5, 7 mod 12), so only an edge is pinned. Each solve takes 15 to 60
seconds, and about 8 minutes for q = 41, which keeps the file out of CI.
"""
import sys, os
import pytest
sys.path.insert(0, os.path.dirname(__file__))

from test_biquadratic_bounds import colourable
from test_finite_planes import plane


@pytest.mark.slow
@pytest.mark.parametrize("q", [29, 31, 41, 43])
def test_the_plane_needs_five_colours(q):
    U, E = plane(q)
    edge = next(e for e in E if e[0] == 0)
    assert not colourable(q * q, E, 4, pin=edge)
