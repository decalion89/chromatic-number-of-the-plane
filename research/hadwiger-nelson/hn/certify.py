"""Certificates: what a third party has to check, and nothing else.

A claim about chi(R^2) is only worth what an independent checker can confirm.
This module produces and validates the two kinds of evidence.

  Upper bound on a graph  -- a colouring.  Checked by re-deriving every
    unit-distance pair *exactly* from the published coordinates, in O(n^2),
    without the spatial hash that built the graph, and confirming no such pair
    is monochromatic.

  Lower bound on a graph  -- a DRAT proof of unsatisfiability, checked by
    drat-trim.  For this to mean anything the checker must also confirm that
    the CNF really encodes the graph, so `verify_uncolorable` re-derives the
    edges from coordinates and rebuilds the formula before checking the proof.
    A proof of the wrong formula proves nothing.

`verify_certificate` runs the whole chain from a JSON file, and is deliberately
written to depend on as little of this package as possible.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import tempfile
from fractions import Fraction
from typing import Dict, List, Optional, Sequence, Tuple

from .field import Field, FieldElement
from .geometry import Point
from .graph import UnitDistanceGraph

__all__ = [
    "save_certificate",
    "load_certificate",
    "exact_edges",
    "verify_coloring",
    "write_dimacs",
    "verify_uncolorable",
    "verify_certificate",
]


# -- serialisation --------------------------------------------------------

def save_certificate(
    graph: UnitDistanceGraph,
    path: str,
    k: int,
    claim: str,
    coloring: Optional[Sequence[int]] = None,
    notes: Optional[dict] = None,
) -> dict:
    """Write the graph, the claim, and any colouring as exact rationals."""
    field = graph.vertices[0].field
    doc = {
        "claim": claim,
        "k": k,
        "field_generators": list(field.gens),
        "basis": ["1"] + [f"sqrt({field._prod[m]})" for m in range(1, field.dim)],
        "n": graph.n,
        "m": graph.m,
        "vertices": [
            {
                "x": [str(c) for c in v.x.c],
                "y": [str(c) for c in v.y.c],
            }
            for v in graph.vertices
        ],
        "coloring": list(coloring) if coloring is not None else None,
        "notes": notes or {},
    }
    with open(path, "w") as f:
        json.dump(doc, f, indent=1)
    return doc


def load_certificate(path: str) -> Tuple[List[Point], dict]:
    with open(path) as f:
        doc = json.load(f)
    field = Field(doc["field_generators"])
    pts = [
        Point(
            field.element([Fraction(s) for s in v["x"]]),
            field.element([Fraction(s) for s in v["y"]]),
        )
        for v in doc["vertices"]
    ]
    return pts, doc


# -- independent re-derivation -------------------------------------------

def exact_edges(points: Sequence[Point], progress: bool = False) -> List[Tuple[int, int]]:
    """Every unit-distance pair, by exact O(n^2) comparison.

    Deliberately brute force: it shares no code path with the spatial hash used
    to build graphs, so agreement between the two is real evidence.
    """
    n = len(points)
    out = []
    for i in range(n):
        pi = points[i]
        for j in range(i + 1, n):
            if pi.is_unit_apart(points[j]):
                out.append((i, j))
        if progress and i % 200 == 0:
            print(f"    exact edges: {i}/{n}", flush=True)
    return out


def verify_coloring(points: Sequence[Point], coloring: Sequence[int], k: int) -> Tuple[bool, str]:
    """Check a proper k-colouring against exactly re-derived edges."""
    if len(coloring) != len(points):
        return False, f"colouring has {len(coloring)} entries for {len(points)} vertices"
    bad = [c for c in coloring if not (0 <= c < k)]
    if bad:
        return False, f"colour {bad[0]} outside 0..{k-1}"
    edges = exact_edges(points)
    for i, j in edges:
        if coloring[i] == coloring[j]:
            return False, f"vertices {i} and {j} are at distance 1 and share colour {coloring[i]}"
    return True, f"proper {k}-colouring; {len(edges)} unit-distance pairs checked exactly"


# -- SAT side --------------------------------------------------------------

def write_dimacs(n: int, edges: Sequence[Tuple[int, int]], k: int, path: str) -> int:
    """The direct k-colouring encoding, written from scratch.

    x[v][c] = 1 + v*k + c.  One 'at least one colour' clause per vertex and one
    'not both' clause per edge and colour.  No symmetry breaking and no vertex
    selectors: a certificate formula must be the plain statement of the claim,
    with nothing extra that a checker would have to trust.
    """
    def x(v, c):
        return 1 + v * k + c

    clauses = []
    for v in range(n):
        clauses.append([x(v, c) for c in range(k)])
    for u, v in edges:
        for c in range(k):
            clauses.append([-x(u, c), -x(v, c)])
    with open(path, "w") as f:
        f.write(f"p cnf {n * k} {len(clauses)}\n")
        for cl in clauses:
            f.write(" ".join(map(str, cl)) + " 0\n")
    return len(clauses)


def verify_uncolorable(
    points: Sequence[Point],
    k: int,
    workdir: Optional[str] = None,
    drat_trim: Optional[str] = None,
    timeout: int = 7200,
) -> Tuple[bool, str]:
    """Prove, checkably, that these points admit no proper k-colouring.

    Edges are re-derived exactly from the coordinates, the CNF is rebuilt from
    those edges, a solver emits a DRAT proof, and drat-trim verifies it against
    that same CNF.  If drat-trim is unavailable the SAT result is still
    reported, clearly marked as unverified.
    """
    from pysat.formula import CNF
    from pysat.solvers import Solver

    tmp = workdir or tempfile.mkdtemp(prefix="hn-cert-")
    os.makedirs(tmp, exist_ok=True)
    cnf_path = os.path.join(tmp, "graph.cnf")
    proof_path = os.path.join(tmp, "graph.drat")

    edges = exact_edges(points)
    n_clauses = write_dimacs(len(points), edges, k, cnf_path)

    cnf = CNF(from_file=cnf_path)
    # Glucose emits the terminating empty clause; CaDiCaL leaves it implicit
    # and drat-trim then refuses the proof, so certificates use Glucose.
    s = Solver(name="g4", bootstrap_with=cnf, with_proof=True)
    try:
        sat = s.solve()
        if sat:
            return False, f"the graph IS {k}-colourable -- no lower bound here"
        proof = s.get_proof()
    finally:
        s.delete()

    with open(proof_path, "w") as f:
        for line in proof:
            f.write(line + "\n")
        # CaDiCaL stops at the last learnt clause and leaves the empty clause
        # implicit; drat-trim needs it stated.  Appending it asserts nothing --
        # drat-trim still has to derive it, and rejects the proof if it cannot.
        if not proof or proof[-1].strip() != "0":
            f.write("0\n")

    checker = drat_trim or shutil.which("drat-trim")
    if not checker:
        return True, (
            f"UNSAT for k={k} over {len(points)} vertices and {len(edges)} exact edges "
            f"({n_clauses} clauses); drat-trim NOT FOUND, so the proof is unverified"
        )
    try:
        r = subprocess.run(
            [checker, cnf_path, proof_path], capture_output=True, text=True, timeout=timeout
        )
    except subprocess.TimeoutExpired:
        return False, "drat-trim timed out"
    out = r.stdout + r.stderr
    if "s VERIFIED" in out:
        return True, (
            f"VERIFIED: no proper {k}-colouring of {len(points)} points with "
            f"{len(edges)} exactly-confirmed unit distances (drat-trim checked "
            f"the DRAT proof of the {n_clauses}-clause formula)"
        )
    return False, f"drat-trim did not verify the proof:\n{out[-2000:]}"


def verify_certificate(path: str, drat_trim: Optional[str] = None) -> Tuple[bool, str]:
    """Check a certificate file end to end."""
    points, doc = load_certificate(path)
    k = doc["k"]
    seen = set()
    for i, p in enumerate(points):
        if p in seen:
            return False, f"vertex {i} is a duplicate"
        seen.add(p)
    if doc.get("coloring") is not None:
        ok, msg = verify_coloring(points, doc["coloring"], k)
        return ok, f"[{doc['claim']}] {msg}"
    return verify_uncolorable(points, k, drat_trim=drat_trim)
