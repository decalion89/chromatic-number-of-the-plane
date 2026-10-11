"""TikZ picture of a configuration: the points of X23 in black, the new points in red, unit-distance edges in grey
(edges at a new point in red). usage: python3 fig_config.py ../../data/forbidden_distances/m1/configurations/x30.json fig-config.tex"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'forbidden_distances', 'm1'))
import exactconf as ec
pts, _ = ec.load(sys.argv[1])
n = len(pts)
xy = [(float(ec.kval(re)), float(ec.kval(im))) for re, im in pts]
d2, edges = ec.graph(pts)
L = [r'\begin{tikzpicture}[scale=4.6,every node/.style={font=\scriptsize}]']
for i, j in edges:
    new = i >= 23 or j >= 23
    st = 'red!70!black,line width=0.5pt' if new else 'black!35,line width=0.4pt'
    L.append(rf'\draw[{st}] ({xy[i][0]:.5f},{xy[i][1]:.5f}) -- ({xy[j][0]:.5f},{xy[j][1]:.5f});')
cx = sum(x for x, _ in xy) / n; cy = sum(y for _, y in xy) / n
for k, (x, y) in enumerate(xy):
    col = 'red!80!black' if k >= 23 else 'black'
    # label away from the centre of the configuration
    dx, dy = x - cx, y - cy; r = (dx * dx + dy * dy) ** 0.5 or 1
    L.append(rf'\fill[{col}] ({x:.5f},{y:.5f}) circle (0.022);')
    L.append(rf'\node[{col},inner sep=1pt] at ({x + 0.075 * dx / r:.5f},{y + 0.075 * dy / r:.5f}) {{{k + 1}}};')
L.append(r'\end{tikzpicture}')
open(sys.argv[2], 'w').write('\n'.join(L) + '\n')
print(n, 'points', len(edges), 'edges')
