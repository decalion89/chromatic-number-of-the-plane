# Exact clique shapes (vertex coordinates as expressions in sqrt, sin, cos, pi and rationals).
import numpy as np
SHAPES = {
    'tri111':  [('0', '0'), ('1', '0'), ('1/2', 'sqrt(3)/2')],
    'tri222':  [('0', '0'), ('2', '0'), ('1', 'sqrt(3)')],
    'col112':  [('0', '0'), ('1', '0'), ('2', '0')],
    'tri122':  [('0', '0'), ('1', '0'), ('1/2', 'sqrt(15)/2')],
    'tri33':   [('0', '0'), ('sqrt(3)', '0'), ('sqrt(3)/2', '3/2')],
    'tri11s3': [('0', '0'), ('1', '0'), ('3/2', 'sqrt(3)/2')],
    'tri1s3s3':[('0', '0'), ('1', '0'), ('1/2', 'sqrt(11)/2')],
    'rhomb1':  [('0', '0'), ('1', '0'), ('1/2', 'sqrt(3)/2'), ('1/2', '-sqrt(3)/2')],
    # golden ratio phi = (1+sqrt(5))/2
    'trippp':  [('0', '0'), ('(1+sqrt(5))/2', '0'), ('(1+sqrt(5))/4', 'sqrt(3)*(1+sqrt(5))/4')],
    'tri11p':  [('0', '0'), ('(1+sqrt(5))/2', '0'), ('(1+sqrt(5))/4', 'sqrt(1-((1+sqrt(5))/4)**2)')],
    'tri1pp':  [('0', '0'), ('1', '0'), ('1/2', 'sqrt(((1+sqrt(5))/2)**2-1/4)')],
    'pent1':   [('cos(2*pi*%d/5)/(2*sin(pi/5))' % k, 'sin(2*pi*%d/5)/(2*sin(pi/5))' % k) for k in range(5)],
}
def ev(e):
    return float(eval(e, {'sqrt': np.sqrt, 'sin': np.sin, 'cos': np.cos, 'pi': np.pi}))

def add_d_shapes(dexpr, tag):
    """Triangles with sides in {1, d}: (d,d,d), (1,1,d) if d < 2, (1,d,d) if 2d > 1."""
    d = '(' + dexpr + ')'
    SHAPES['tri_ddd_' + tag] = [('0', '0'), (d, '0'), (d + '/2', 'sqrt(3)*' + d + '/2')]
    if ev(dexpr) < 2:
        SHAPES['tri_11d_' + tag] = [('0', '0'), (d, '0'), (d + '/2', 'sqrt(1-' + d + '**2/4)')]
    SHAPES['tri_1dd_' + tag] = [('0', '0'), ('1', '0'), ('1/2', 'sqrt(' + d + '**2-1/4)')]
    return [n for n in SHAPES if n.endswith('_' + tag)]

import os
if os.environ.get('DSHAPE'):
    add_d_shapes(os.environ['DSHAPE'], os.environ.get('DTAG', 'x'))
