"""Command line entry points.

    python -m hn.cli verify   CERT [--drat-trim PATH]
    python -m hn.cli spindle  --k K [--m-max M] [--steps S] [--radius R]
    python -m hn.cli demo
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys

from .certify import save_certificate, verify_certificate
from .coloring import chromatic_number, find_uncolorable_core, is_k_colorable
from .field import QSQRT3_11 as FIELD
from .generate import unit_vectors
from .geometry import SPINDLE, eisenstein, origin
from .graph import build_graph


def cmd_verify(args) -> int:
    checker = args.drat_trim or shutil.which("drat-trim")
    ok, msg = verify_certificate(args.certificate, drat_trim=checker)
    print(("VERIFIED: " if ok else "REJECTED: ") + msg)
    return 0 if ok else 1


def cmd_demo(args) -> int:
    """Rebuild the classic results from scratch, as a self-check."""
    rhombus = [origin(), eisenstein(1, 0), eisenstein(0, 1), eisenstein(1, 1)]
    g = build_graph(rhombus + [SPINDLE(p) for p in rhombus])
    print(f"Moser spindle: {g}")
    print(f"  chi = {chromatic_number(g)[0]}  (expected 4)")
    core = find_uncolorable_core(g, 3, verbose=False)
    print(f"  minimal non-3-colourable subgraph: {core}")
    out = args.out or "certificates"
    os.makedirs(out, exist_ok=True)
    p = os.path.join(out, "moser_spindle_no3coloring.json")
    save_certificate(g, p, 3, "chi(R^2) >= 4: this unit-distance graph has no proper 3-colouring")
    print(f"  wrote {p}")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="hn", description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    v = sub.add_parser("verify", help="check a certificate file end to end")
    v.add_argument("certificate")
    v.add_argument("--drat-trim", default=None, help="path to the drat-trim binary")
    v.set_defaults(func=cmd_verify)

    dg = sub.add_parser("degrey", help="rebuild de Grey's 1581-vertex graph and certify chi >= 5")
    dg.add_argument("--out", default=None)
    dg.add_argument("--timeout", type=float, default=3600)
    dg.set_defaults(func=cmd_degrey)

    d = sub.add_parser("demo", help="rebuild and certify the classic chi >= 4 result")
    d.add_argument("--out", default=None)
    d.set_defaults(func=cmd_demo)

    args = ap.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())


def cmd_degrey(args) -> int:
    """Rebuild de Grey's 1581-vertex graph and certify chi(R^2) >= 5."""
    import os

    from .certify import save_certificate
    from .coloring import is_k_colorable
    from .degrey import build_G

    g = build_G()
    print(f"de Grey 1581: {g}  (paper: n=1581, m=7877)")
    sat, _ = is_k_colorable(g, 4, timeout=args.timeout)
    print(f"  4-colourable: {sat}")
    if sat is not False:
        print("  reconstruction did NOT reproduce the result")
        return 1
    out = args.out or "certificates"
    os.makedirs(out, exist_ok=True)
    p = os.path.join(out, "degrey_1581_no4coloring.json")
    save_certificate(
        g, p, 4,
        "chi(R^2) >= 5: this 1581-vertex unit-distance graph has no proper 4-colouring",
        notes={"source": "de Grey, arXiv:1804.02385", "rebuilt_from": "published recipe"},
    )
    print(f"  wrote {p}")
    return 0
