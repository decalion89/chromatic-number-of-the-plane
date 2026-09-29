"""tabu.py -- python wrapper for ./tabucol"""
import os, subprocess, tempfile
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def tabucol(n, E, init=None, fixed=(), k=4, maxiter=2_000_000, seed=1, timeout=None):
    """E: (m,2) int array; init: length-n colours (-1 = none); fixed: [(v,c)]. Returns colour array or None."""
    if init is None:
        init = -np.ones(n, dtype=np.int64)
    fd, path = tempfile.mkstemp(dir="/dev/shm", prefix="tabu_", suffix=".txt")
    with os.fdopen(fd, "w") as f:
        f.write(f"{n} {len(E)} {k} {maxiter} {seed} {len(fixed)}\n")
        for v, c in fixed:
            f.write(f"{v} {c}\n")
        f.write("\n".join(str(int(c)) for c in init) + "\n")
        np.savetxt(f, np.asarray(E, dtype=np.int64), fmt="%d")
    try:
        out = subprocess.run(["nice", "-n", "19", os.path.join(HERE, "tabucol"), path], capture_output=True,
                             text=True, timeout=timeout).stdout.split()
    except subprocess.TimeoutExpired:
        out = ["FAIL", "-1"]
    finally:
        os.unlink(path)
    if out and out[0] == "OK":
        col = np.array(list(map(int, out[1:])), dtype=np.int64)
        assert len(col) == n
        E = np.asarray(E)
        assert np.all(col[E[:, 0]] != col[E[:, 1]])
        return col
    return None


if __name__ == "__main__":
    import sys, time, json
    from flat import directions, build_edges
    D, CD, U = directions()
    P = np.load(sys.argv[1])
    st = json.load(open(sys.argv[2])) if len(sys.argv) > 2 else None
    E, J = build_edges(P, U)
    fixed = [(v, i) for i, v in enumerate(st["triangle"])] if st else []
    for seed in range(1, 4):
        t = time.time()
        col = tabucol(len(P), E, fixed=fixed, maxiter=int(sys.argv[3]) if len(sys.argv) > 3 else 5_000_000, seed=seed)
        print(f"n={len(P)} m={len(E)} seed {seed}: {'found' if col is not None else 'FAIL'} in {time.time() - t:.1f}s", flush=True)
