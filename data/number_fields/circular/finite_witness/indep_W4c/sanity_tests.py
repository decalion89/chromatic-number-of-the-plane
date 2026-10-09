#!/usr/bin/env python3
"""Mutation tests: the geometry checker must REJECT corrupted witnesses.

Usage: python3 sanity_tests.py WITNESS.json.gz
Each test mutates a deep copy of the witness and asserts that check_all
reports at least one failing check, and (where meaningful) that the
expected specific check is among the failures.
"""
import copy
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_geom import load, check_all  # noqa: E402

random.seed(2024)


def run(data, **kw):
    res, stats = check_all(data, do_unit=False, log=lambda *a: None, **kw)
    return {k for k, v in res.items() if not v}, stats


def main():
    path = sys.argv[1]
    base = load(path)
    fails0, _ = run(base)
    print("baseline failures:", sorted(fails0))
    assert not fails0, "baseline must pass"
    n = len(base["points"])
    E = set((min(e), max(e)) for e in base["edges"])
    P = [tuple(p) for p in base["points"]]
    gpm = set()
    for g in base["generators"]:
        gpm.add(tuple(g))
        gpm.add(tuple(-x for x in g))
    adj = {i: set() for i in range(n)}
    for (i, j) in E:
        adj[i].add(j)
        adj[j].add(i)
    fv = base["fixed_vertex"]
    results = []

    def expect(name, data, must_fail):
        fails, _ = run(data)
        ok = bool(fails) and all(m in fails for m in must_fail)
        results.append((name, ok, sorted(fails)))
        print("%-70s %s   failing checks: %s" % (name, "REJECTED (good)" if ok else "NOT REJECTED / WRONG REASON (bad)", sorted(fails)))

    # 1. moved points
    for trial in range(3):
        v = random.choice([i for i in range(n) if i != fv])
        coord = random.randrange(8)
        d = copy.deepcopy(base)
        d["points"][v][coord] += 1
        expect("moved point %d (coordinate %d += 1)" % (v, coord), d,
               ["claim3_listed_edges_are_generator_differences"])
    # moved point by a generator to an unoccupied place
    for trial in range(200):
        v = random.choice([i for i in range(n) if i != fv])
        g = random.choice(sorted(gpm))
        q = tuple(a + b for a, b in zip(P[v], g))
        if q not in set(P):
            d = copy.deepcopy(base)
            d["points"][v] = list(q)
            expect("moved point %d by a generator to an empty site" % v, d,
                   ["claim3_listed_edges_are_generator_differences"])
            break
    d = copy.deepcopy(base)
    d["points"][fv][0] += 1
    expect("moved the fixed vertex (origin)", d, ["claim3_fixed_vertex_is_origin"])
    # moved point onto another point
    d = copy.deepcopy(base)
    v, w = 0, 1
    d["points"][v] = list(d["points"][w])
    expect("point 0 moved onto point 1 (duplicate)", d, ["points_distinct"])

    # 2. an edge whose difference is not a generator
    for trial in range(1000):
        i, j = random.randrange(n), random.randrange(n)
        if i != j and (min(i, j), max(i, j)) not in E:
            dd = tuple(b - a for a, b in zip(P[i], P[j]))
            assert dd not in gpm
            d = copy.deepcopy(base)
            d["edges"].append([i, j])
            expect("added non-generator edge (%d,%d)" % (i, j), d,
                   ["claim3_listed_edges_are_generator_differences"])
            break
    # an added edge at unit distance whose difference is not a generator
    from biquad import Field, int_sq_norm_8
    A, B = (3, 11) if "11" in base["field"] else (2, 3)
    F = Field(A, B)
    D = base["denominator"]
    found = None
    for i in range(n):
        for j in range(i + 1, n):
            dd = tuple(b - a for a, b in zip(P[i], P[j]))
            if dd not in gpm and int_sq_norm_8(F, dd) == (D * D, 0, 0, 0):
                found = (i, j)
                break
        if found:
            break
    if found:
        d = copy.deepcopy(base)
        d["edges"].append(list(found))
        expect("added UNIT-DISTANCE non-generator edge %s" % (found,), d,
               ["claim3_listed_edges_are_generator_differences"])
    # self loop
    d = copy.deepcopy(base)
    d["edges"].append([5, 5])
    expect("added self-loop (5,5)", d, ["edges_well_formed"])

    # 3. missing Cayley edge
    for trial in range(3):
        k = random.randrange(len(base["edges"]))
        d = copy.deepcopy(base)
        e = d["edges"].pop(k)
        expect("removed Cayley edge %s" % (e,), d, ["claim3_every_cayley_pair_is_listed"])

    # 4. improper colouring
    for trial in range(3):
        (i, j) = random.choice(sorted(E))
        d = copy.deepcopy(base)
        d["colouring"][i] = d["colouring"][j]
        expect("improper colouring on edge (%d,%d)" % (i, j), d, ["claim4_proper_4_colouring"])
    d = copy.deepcopy(base)
    d["colouring"][3] = 4
    expect("colour value 4 out of range", d, ["claim4_proper_4_colouring"])
    d = copy.deepcopy(base)
    d["colouring"].pop()
    expect("colouring list one too short", d, ["claim4_proper_4_colouring"])

    # 5. cycle through a non-edge
    for trial in range(3):
        ci = random.randrange(len(base["cycles"]))
        cyc = base["cycles"][ci]
        m = len(cyc)
        t = random.randrange(m)
        # replace vertex t by some vertex that is not adjacent to its successor
        nxt = cyc[(t + 1) % m]
        cand = [u for u in range(n) if u not in cyc and u not in adj[nxt]]
        u = random.choice(cand)
        d = copy.deepcopy(base)
        d["cycles"][ci] = cyc[:t] + [u] + cyc[t + 1:]
        expect("cycle %d with vertex %d replaced by %d (non-edge step)" % (ci, cyc[t], u), d,
               ["claim5_cycles_valid"])
    # a 4-"cycle" made of a path i-j-k-l whose closing pair is a non-edge
    d = copy.deepcopy(base)
    i = 0
    j = sorted(adj[i])[0]
    k = [x for x in sorted(adj[j]) if x != i][0]
    l = [x for x in sorted(adj[k]) if x not in (i, j) and x not in adj[i]][0]
    d["cycles"].append([i, j, k, l])
    expect("appended walk %s whose closing step is a non-edge" % ([i, j, k, l],), d, ["claim5_cycles_valid"])
    # repeated vertex / bad length / empty cycle
    c0 = base["cycles"][0]
    d = copy.deepcopy(base)
    d["cycles"].append(c0 + c0)
    expect("cycle traversed twice (repeated vertices)", d, ["claim5_cycles_valid"])
    d = copy.deepcopy(base)
    d["cycles"].append([])
    expect("empty cycle", d, ["claim5_cycles_valid"])
    # a closed walk of length 6 along edges (if a hexagon exists we would use it; use a triangle-free check: build length-6 walk u v u' ... not simple)
    # instead: take a 4-cycle and insert a back-and-forth detour -> repeated vertex & length 6
    c4 = next(c for c in base["cycles"] if len(c) == 4)
    d = copy.deepcopy(base)
    d["cycles"].append([c4[0], c4[1], c4[0], c4[1], c4[2], c4[3]])
    expect("length-6 closed walk with repeats", d, ["claim5_cycles_valid"])

    # 6. generator corruptions
    d = copy.deepcopy(base)
    d["generators"][0] = [x * 2 for x in d["generators"][0]]
    expect("generator 0 doubled (length 2)", d, ["claim3_generators_unit_length", "claim12_generators_match_construction"])
    if found:
        dd = tuple(b - a for a, b in zip(P[found[0]], P[found[1]]))
        d = copy.deepcopy(base)
        d["generators"][-1] = list(dd)
        expect("last generator replaced by another unit vector", d, ["claim12_generators_match_construction"])
    d = copy.deepcopy(base)
    d["generators"].pop()
    expect("one generator dropped", d, ["claim12_generators_match_construction"])

    # 7. unreachable point
    d = copy.deepcopy(base)
    far = [10 ** 6 + 1] + [0] * 7
    d["points"].append(far)
    d["colouring"].append(0)
    expect("added isolated far point", d, ["claim3_all_points_reachable"])
    # 8. fixed vertex not the origin (point the index elsewhere)
    d = copy.deepcopy(base)
    d["fixed_vertex"] = (fv + 1) % n
    expect("fixed_vertex index moved to another point", d, ["claim3_fixed_vertex_is_origin"])

    nbad = sum(1 for r in results if not r[1])
    print("SANITY SUMMARY: %d tests, %d correctly rejected, %d not rejected" % (len(results), len(results) - nbad, nbad))
    sys.exit(0 if nbad == 0 else 1)


if __name__ == "__main__":
    main()
