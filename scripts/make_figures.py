"""Draw the README figures from exact data: the 2-adic 4-colourings in action, and the two lower bounds.

Writes SVG files to docs/figures/:
  plane_q311.svg       unit-distance graph on a piece of Q(sqrt3, sqrt11)^2, 4-coloured by hn.adelic.q311_colour
  lower_bounds.svg     the Moser spindle (Q(sqrt3, sqrt11)) and the 10-vertex rhombus chain (Q(sqrt2, sqrt3)),
                       each properly 4-coloured by its 2-adic colouring
  quadratic_q11.svg    the 94-vertex graph over Q(sqrt11) of notes/quadratic_planes.md, with no 3-colouring,
                       coloured by its stored 4-colouring (data/quadratic_planes/q11.json)

Every edge is found by exact arithmetic (hn.graph.build_graph), and the script asserts that no edge
is monochromatic before drawing. Colours are the four residue pairs (rho(alpha), rho(beta)) in F_2^2,
each with its own marker shape so identity never rests on colour alone.

usage: python3 scripts/make_figures.py
"""
import json
import os
import sys
from fractions import Fraction as Fr

HN_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HN_DIR)

from hn.adelic import q23_colour, q311_colour  # noqa: E402
from hn.field import Field  # noqa: E402
from hn.geometry import Point, Rotation  # noqa: E402
from hn.graph import build_graph  # noqa: E402

OUT = os.path.join(HN_DIR, "docs", "figures")

# Categorical slots 1, 2, 3 and 7 of the reference palette; this set passes the all-pairs
# colour-vision checks on the light surface, and every class also has its own marker shape.
COLOURS = ["#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"]
SHAPES = ["circle", "square", "triangle", "diamond"]
LABELS = ["(0, 0)", "(1, 0)", "(0, 1)", "(1, 1)"]
SURFACE, INK, INK2, EDGE = "#fcfcfb", "#0b0b0b", "#52514e", "#c9c8c1"
FONT = "system-ui, -apple-system, 'Segoe UI', sans-serif"


def marker(shape, x, y, r, fill):
    ring = f'stroke="{SURFACE}" stroke-width="1.5"'
    if shape == "circle":
        return f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{r:.2f}" fill="{fill}" {ring}/>'
    if shape == "square":
        s = r * 0.9
        return f'<rect x="{x - s:.2f}" y="{y - s:.2f}" width="{2 * s:.2f}" height="{2 * s:.2f}" rx="1" fill="{fill}" {ring}/>'
    if shape == "triangle":
        h = r * 1.25
        pts = f"{x:.2f},{y - h:.2f} {x + h * 0.95:.2f},{y + h * 0.7:.2f} {x - h * 0.95:.2f},{y + h * 0.7:.2f}"
        return f'<polygon points="{pts}" fill="{fill}" {ring}/>'
    d = r * 1.25
    pts = f"{x:.2f},{y - d:.2f} {x + d:.2f},{y:.2f} {x:.2f},{y + d:.2f} {x - d:.2f},{y:.2f}"
    return f'<polygon points="{pts}" fill="{fill}" {ring}/>'


def panel(g, colours, x0, y0, w, h, r, title, subtitle, highlight=(), labels=None):
    """One panel: edges, then markers, fitted into the box (x0, y0, w, h)."""
    xs = [p.fx for p in g.vertices]
    ys = [p.fy for p in g.vertices]
    pad = 18
    sx = (w - 2 * pad) / max(max(xs) - min(xs), 1e-9)
    sy = (h - 2 * pad - 34) / max(max(ys) - min(ys), 1e-9)
    s = min(sx, sy)
    cx = x0 + w / 2 - s * (max(xs) + min(xs)) / 2
    cy = y0 + 34 + (h - 34) / 2 + s * (max(ys) + min(ys)) / 2
    X = [cx + s * x for x in xs]
    Y = [cy - s * y for y in ys]
    out = [f'<text x="{x0 + w / 2:.1f}" y="{y0 + 14:.1f}" text-anchor="middle" font-size="13" '
           f'font-weight="600" fill="{INK}">{title}</text>',
           f'<text x="{x0 + w / 2:.1f}" y="{y0 + 30:.1f}" text-anchor="middle" font-size="11" '
           f'fill="{INK2}">{subtitle}</text>']
    for i in range(len(X)):
        for j in g.adj[i]:
            if i < j:
                hl = (i, j) in highlight or (j, i) in highlight
                style = f'stroke="{INK2}" stroke-width="2" stroke-dasharray="5 3"' if hl else \
                    f'stroke="{EDGE}" stroke-width="1"'
                out.append(f'<line x1="{X[i]:.2f}" y1="{Y[i]:.2f}" x2="{X[j]:.2f}" y2="{Y[j]:.2f}" {style}/>')
    for i in range(len(X)):
        c = colours[i]
        out.append(marker(SHAPES[c], X[i], Y[i], r, COLOURS[c]))
    for i, (text, dx, dy) in (labels or {}).items():
        out.append(f'<text x="{X[i] + dx:.1f}" y="{Y[i] + dy:.1f}" font-size="12" font-style="italic" '
                   f'fill="{INK}">{text}</text>')
    return out


def legend(x0, y0):
    out = [f'<text x="{x0}" y="{y0}" font-size="11" fill="{INK2}">colour = residue pair (ρ(α), ρ(β)):</text>']
    x = x0 + 196
    for c in range(4):
        out.append(marker(SHAPES[c], x, y0 - 4, 5, COLOURS[c]))
        out.append(f'<text x="{x + 9}" y="{y0}" font-size="11" fill="{INK}">{LABELS[c]}</text>')
        x += 62
    return out


def svg(width, height, body, label):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" width="{width}" '
            f'height="{height}" role="img" aria-label="{label}" font-family="{FONT}">\n'
            f'<rect width="{width}" height="{height}" rx="10" fill="{SURFACE}"/>\n' + "\n".join(body) + "\n</svg>\n")


def proper(g, colours):
    return all(colours[i] != colours[j] for i in range(g.n) for j in g.adj[i])


def plane_figure():
    F = Field((3, 11))
    half, s3, s11 = F.rational(Fr(1, 2)), F.sqrt(3), F.sqrt(11)
    rot60 = Rotation(half, s3 * half)
    moser = Rotation(F.rational(Fr(5, 6)), s11 * F.rational(Fr(1, 6)))
    units, u = [], Point(F.one(), F.zero())
    for _ in range(6):
        units += [u, moser(u), moser.inverse()(u)]
        u = rot60(u)
    origin = Point(F.zero(), F.zero())
    pts = {origin} | set(units) | {a + b for a in units for b in units}
    g = build_graph(pts)
    colours = [q311_colour(p) for p in g.vertices]
    assert proper(g, colours), "the 2-adic colouring must be proper"
    body = panel(g, colours, 0, 8, 720, 560, 4.2,
                 "A piece of the plane over ℚ(√3, √11), coloured with four colours",
                 f"{g.n} points, {g.m} unit-distance edges; no edge joins two points of the same colour")
    body += legend(150, 588)
    return svg(720, 604, body, "Unit-distance graph on a piece of the plane over Q(sqrt3, sqrt11), 4-coloured"), g


def lower_bounds_figure():
    # Moser spindle in Q(sqrt3, sqrt11)
    F = Field((3, 11))
    half, s3, s11 = F.rational(Fr(1, 2)), F.sqrt(3), F.sqrt(11)
    rhombus = [Point(F.zero(), F.zero()), Point(F.one(), F.zero()), Point(half, s3 * half),
               Point(F.rational(Fr(3, 2)), s3 * half)]
    moser = Rotation(F.rational(Fr(5, 6)), s11 * F.rational(Fr(1, 6)))
    g1 = build_graph(rhombus + [moser(p) for p in rhombus])
    c1 = [q311_colour(p) for p in g1.vertices]
    assert g1.n == 7 and g1.m == 11 and proper(g1, c1)
    tips = [g1.vertices.index(Point(F.rational(Fr(3, 2)), s3 * half)),
            g1.vertices.index(moser(Point(F.rational(Fr(3, 2)), s3 * half)))]
    # The rhombus chain in Q(sqrt2, sqrt3)
    d = json.load(open(os.path.join(HN_DIR, "data", "chain23.json")))
    G = Field((2, 3))
    pts = [Point(G.element([Fr(a, b) for a, b in x]), G.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
    g2 = build_graph(pts)
    c2 = [q23_colour(p, g2.vertices[0]) for p in g2.vertices]
    assert g2.n == 10 and g2.m == 16 and proper(g2, c2)
    o = g2.vertices.index(Point(G.zero(), G.zero()))
    s2, s3b = G.sqrt(2), G.sqrt(3)
    u3 = Point(G.rational(Fr(-2, 3)) - s2 * G.rational(Fr(1, 6)), G.rational(Fr(-2, 3)) + s2 * G.rational(Fr(1, 6)))
    end = [g2.vertices.index(Point(s3b * (G.one() + u3.x), s3b * (G.one() + u3.y)))]
    assert end[0] in g2.adj[o]
    body = panel(g1, c1, 0, 8, 360, 330, 6, "Four colours are needed: ℚ(√3, √11)",
                 "the Moser spindle, 7 vertices", highlight=[tuple(tips)],
                 labels={tips[0]: ("A", 9, 4), tips[1]: ("B", 9, 4)})
    body += panel(g2, c2, 370, 8, 360, 330, 6, "Four colours are needed: ℚ(√2, √3)",
                  "a chain of three unit rhombi, 10 vertices", highlight=[(o, end[0])],
                  labels={o: ("P₀", -22, 5), end[0]: ("P₃", -22, 0)})
    body.append(f'<text x="365" y="352" text-anchor="middle" font-size="11" fill="{INK2}">Dashed: an edge whose two '
                f'ends every 3-colouring forces to share a colour, so three colours cannot suffice.</text>')
    body += legend(150, 376)
    return svg(730, 392, body, "The Moser spindle and the 10-vertex rhombus chain, each 4-coloured"), (g1, g2)


def quadratic_figure():
    """The vertex-critical graph over Q(sqrt11): points from the exact data, edges checked to length 1 (in floating
    point here; the checker does it exactly), the stored 4-colouring asserted proper."""
    import math
    g = json.load(open(os.path.join(HN_DIR, "data", "quadratic_planes", "q11.json")))
    r = math.sqrt(g["d"])
    P = [((a + b * r) / g["D"], (c + e * r) / g["D"]) for a, b, c, e in g["points"]]
    E = [tuple(e) for e in g["edges"]]
    col = [int(c) for c in g["four_colouring"]]
    assert all(abs(math.dist(P[a], P[b]) - 1) < 1e-9 for a, b in E)
    assert all(col[a] != col[b] for a, b in E), "the stored 4-colouring must be proper"
    W, H, pad, top = 720, 640, 22, 46
    xs, ys = [p[0] for p in P], [p[1] for p in P]
    s = min((W - 2 * pad) / (max(xs) - min(xs)), (H - top - 2 * pad - 30) / (max(ys) - min(ys)))
    cx = W / 2 - s * (max(xs) + min(xs)) / 2
    cy = top + pad + (H - top - 2 * pad - 30) / 2 + s * (max(ys) + min(ys)) / 2
    X = [cx + s * x for x in xs]
    Y = [cy - s * y for y in ys]
    body = [f'<text x="{W / 2:.1f}" y="22" text-anchor="middle" font-size="13" font-weight="600" fill="{INK}">'
            f'A plane over a real quadratic field that needs four colours: ℚ(√11)</text>',
            f'<text x="{W / 2:.1f}" y="38" text-anchor="middle" font-size="11" fill="{INK2}">{len(P)} points, '
            f'{len(E)} segments of length exactly 1; no 3-colouring exists. One 4-colouring is shown.</text>',
            f'<circle cx="{cx:.2f}" cy="{cy:.2f}" r="{s:.2f}" fill="none" stroke="{INK2}" stroke-width="1" '
            f'stroke-dasharray="4 4" opacity="0.6"/>']
    body += [f'<line x1="{X[a]:.2f}" y1="{Y[a]:.2f}" x2="{X[b]:.2f}" y2="{Y[b]:.2f}" stroke="{EDGE}" '
             f'stroke-width="1"/>' for a, b in E]
    body += [marker(SHAPES[col[v]], X[v], Y[v], 4.2, COLOURS[col[v]]) for v in range(len(P))]
    body.append(f'<text x="{W / 2:.1f}" y="{H - 32}" text-anchor="middle" font-size="11" fill="{INK2}">Dashed: the '
                f'unit circle around the point at the origin, which has 32 neighbours on it.</text>')
    body.append(f'<text x="{W / 2:.1f}" y="{H - 16}" text-anchor="middle" font-size="11" fill="{INK2}">Shapes and '
                f'colours: the four colour classes of a proper 4-colouring.</text>')
    return svg(W, H, body, "The 94-vertex unit-distance graph over Q(sqrt11), which needs four colours"), (len(P), len(E))


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    s, g = plane_figure()
    open(os.path.join(OUT, "plane_q311.svg"), "w").write(s)
    t, _ = lower_bounds_figure()
    open(os.path.join(OUT, "lower_bounds.svg"), "w").write(t)
    q, (n, m) = quadratic_figure()
    open(os.path.join(OUT, "quadratic_q11.svg"), "w").write(q)
    print(f"plane_q311.svg: {g.n} points, {g.m} edges; lower_bounds.svg; quadratic_q11.svg: {n} points, {m} "
          f"edges; written to {OUT}")
