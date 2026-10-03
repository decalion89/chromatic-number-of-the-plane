"""fastcert.py: exact certificate builder, faster variant of certify_w2.py.
Same certificate format (checked by check_w.py).  Differences: per node one LP feasibility test plus LP ranges only
for a sample of the unfixed relations (the K shortest by |r|_1 and K random), and children enumerated in order of
distance from the LP-midpoint; every leaf still carries an exact Farkas vector."""
import sys, json, time, math, random
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from certify_w2 import relation_basis, box_range, farkas, lp_feasible, lp_range, LO, HI
K = int(sys.argv[3]) if len(sys.argv) > 3 else 8
def main():
    inp = json.load(open(sys.argv[1])); V = [tuple(u) for u in inp["units"]]; n = len(V)
    R = relation_basis(V); R.sort(key=lambda r: sum(abs(c) for c in r))
    rng = random.Random(1)
    st = {"nodes": 0, "leaves": 0, "last": time.time()}; t0 = time.time()
    def node(fixed):
        st["nodes"] += 1
        if time.time() - st["last"] > 30:
            st["last"] = time.time()
            print(f"  nodes {st['nodes']} leaves {st['leaves']} depth {len(fixed)} [{time.time()-t0:.0f}s]", flush=True)
        rows = [R[j] for j, _ in fixed]; rhs = [z for _, z in fixed]
        if rows and not lp_feasible(rows, rhs, n):
            y = farkas(rows, rhs)
            if y is None: raise RuntimeError("no exact Farkas vector")
            st["leaves"] += 1
            return {"leaf": [str(v) for v in y]}
        done = {j for j, _ in fixed}
        free = [j for j in range(len(R)) if j not in done]
        if not free: raise RuntimeError("feasible leaf")
        cand = free[:K] + rng.sample(free[K:], min(K, max(0, len(free) - K)))
        best = None
        for j in cand:
            rg = lp_range(rows, rhs, n, R[j])
            cnt = 0 if rg is None else max(0, math.floor(rg[1] + 1e-9) - math.ceil(rg[0] - 1e-9) + 1)
            if best is None or cnt < best[0]: best = (cnt, j)
            if cnt <= 1: break
        j = best[1]; lo, hi = box_range(R[j])
        kids = {}
        for z in range(lo, hi + 1):
            kids[str(z)] = node(fixed + [(j, z)])
        return {"branch": j, "lo": lo, "hi": hi, "kids": kids}
    tree = node([])
    json.dump({"d": inp["d"], "D": inp["D"], "units": [list(u) for u in V], "relations": R, "tree": tree}, open(sys.argv[2], "w"))
    print(f"certificate written: {st['nodes']} nodes, {st['leaves']} leaves, {len(R)} relations  [{time.time()-t0:.0f}s]")
main()
