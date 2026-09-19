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

    d = sub.add_parser("demo", help="rebuild and certify the classic chi >= 4 result")
    d.add_argument("--out", default=None)
    d.set_defaults(func=cmd_demo)

    args = ap.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
