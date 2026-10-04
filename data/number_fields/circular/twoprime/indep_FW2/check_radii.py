"""Check Lemma 9 (lem:radii) bounds used for r_theta < 0.23 on (2/7, 1/3]:
 (i) max |x| over P_1 (s = 1) equals sqrt(10)/3;  (ii) the hexagon X has vertices of modulus <= sqrt(2);
 and max over theta in (2/7,1/3] of r_theta = max(s*sqrt(10)/3, sqrt(2)*delta) < 0.23."""
from fractions import Fraction as Fr
from itertools import combinations
import math
rho = (Fr(3,5), Fr(4,5))
def cmul(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def conj(a): return (a[0], -a[1])
def rpow(j):
    g=(Fr(1),Fr(0)); base = rho if j>=0 else conj(rho)
    for _ in range(abs(j)): g=cmul(g,base)
    return g
# (i) P_1 with s=1: constraints  +-(xbar rho^j)_e <= 1, j in {-1,0,1}, e in {1,2}; x = (x1,x2), xbar=(x1,-x2)
# (xbar*g)_1 = x1 g1 + x2 g2 ; (xbar*g)_2 = x1 g2 - x2 g1
cons=[]
for j in (-1,0,1):
    g=rpow(j)
    for row in ((g[0],g[1]),(g[1],-g[0])):
        for sg in (1,-1):
            cons.append(((sg*row[0], sg*row[1]), Fr(1)))
def vertices(cons):
    V=[]
    for (a,b1),(c,b2) in combinations(cons,2):
        det=a[0]*c[1]-a[1]*c[0]
        if det==0: continue
        x=((b1*c[1]-a[1]*b2)/det, (a[0]*b2-b1*c[0])/det)
        if all(r[0]*x[0]+r[1]*x[1] <= rhs for r,rhs in cons):
            V.append(x)
    return set(V)
VP=vertices(cons)
print('P_1 vertices:', sorted((str(v[0]),str(v[1])) for v in VP))
print('max |x|^2 over P_1 =', max(v[0]**2+v[1]**2 for v in VP), ' (10/9 =', Fr(10,9),')')
# (ii) hexagon X: eta_{j,e} (ybar rho^j)_e <= 1 for |j|<=1, eta_-1=1+i, eta_0=1-i, eta_1=-1-i
eta={-1:(1,1),0:(1,-1),1:(-1,-1)}
cons2=[]
for j in (-1,0,1):
    g=rpow(j)
    rows=((g[0],g[1]),(g[1],-g[0]))
    for e in (0,1):
        cons2.append(((eta[j][e]*rows[e][0], eta[j][e]*rows[e][1]), Fr(1)))
VX=vertices(cons2)
print('X vertices:', sorted((str(v[0]),str(v[1])) for v in VX))
print('max |y|^2 over X =', max(v[0]**2+v[1]**2 for v in VX))
# r_theta sup on (2/7,1/3]: s<3/14, delta<1/21
s_sup=Fr(3,14); d_sup=Fr(1,21)
print('sup s*sqrt10/3 =', float(s_sup)*math.sqrt(10)/3, ' sup sqrt2*delta =', math.sqrt(2)*float(d_sup))
