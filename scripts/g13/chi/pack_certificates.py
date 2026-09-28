"""Pack the logs of the proof that G_13 has no proper 5-colouring (notes/g13_chi.md) into one archive.

usage: python3 scripts/g13/chi/pack_certificates.py ASSEMBLY [--out-dir DIR]

ASSEMBLY is a merged copy of this directory: every share's cases/*.certlog, the cube files of the re-splits
(cases/*_leaf*.icnf, cases/*_split*.icnf) and logs/E37_A.log, merged as README.md ("Checking the result") says.
The script writes two files into DIR (by default certificates/ of the repository):
  - g13_chi_certlogs.tar.gz, the archive of exactly those files, under the same names;
  - g13_chi_SHA256SUMS.txt, the SHA-256 of the archive and then of every file in it, in the format of sha256sum
    (`sha256sum -c --ignore-missing` checks the archive in DIR, and the files where the archive is unpacked).
The archive is deterministic: the names sorted, regular files only, mode 644, owner and group 0 with no names, the
modification time 2026-09-27 00:00 UTC, GNU tar format, and gzip at level 9 with no file name and no time stamp. So
the same files give the same archive, with the same zlib.
Every line of every log must be a line of certify.py, and name a leaf of that log's cube file. The four base cube
files (cases/E37_B.icnf, F36.icnf, F35.icnf, F34.icnf) and cases/SHA256SUMS are not packed: they are in the
repository, and the assembly must have the same ones. Formulas (*.cnf) are not packed: `share_driver.py regen`
writes them again.
Check the result with: python3 scripts/verify_g13_chi.py DIR/g13_chi_certlogs.tar.gz --stats
"""
import argparse
import gzip
import hashlib
import io
import os
import re
import sys
import tarfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
ARCHIVE, SUMS = "g13_chi_certlogs.tar.gz", "g13_chi_SHA256SUMS.txt"
MTIME = 1790467200                                      # 2026-09-27 00:00:00 UTC
BASE = ["cases/E37_B.icnf", "cases/F36.icnf", "cases/F35.icnf", "cases/F34.icnf", "cases/SHA256SUMS"]
CASE = r"(?:E37_B|F36|F35|F34)"
PACKED = [re.compile(rf"cases/{CASE}(?:_(?:leaf|split)\d+)*\.certlog"),      # the logs of certify.py
          re.compile(rf"cases/{CASE}(?:_(?:leaf|split)\d+)+\.icnf"),          # the cube files of the re-splits
          re.compile(r"logs/E37_A\.log")]
LINE = re.compile(r"(\S+): sha256 [0-9a-f]{64}; kissat (?:UNSAT|SAT|UNKNOWN) in [\d.]+ s \(proof [\d.]+ MB\); "
                  r"drat-trim (?:VERIFIED|NOT VERIFIED|not run)(?:, \d+ of \d+ lemmas in core)? in [\d.]+ s$")


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def members(assembly):
    """the files of the assembly that go into the archive (names relative to it, sorted), and the others"""
    found, other = [], []
    for sub in ("cases", "logs"):
        folder = os.path.join(assembly, sub)
        for name in sorted(os.listdir(folder)) if os.path.isdir(folder) else []:
            rel = f"{sub}/{name}"
            if not os.path.isfile(os.path.join(assembly, rel)):
                other.append(rel)
            elif any(p.fullmatch(rel) for p in PACKED):
                found.append(rel)
            elif rel not in BASE:
                other.append(rel)
    return sorted(found), other


def check_log(rel, data):
    """every line is a line of certify.py, for a leaf of this log's cube file (or, in logs/E37_A.log, for one of
    the four formulas E37_A_c)"""
    stem = os.path.basename(rel)[:-len(".certlog")] if rel.endswith(".certlog") else None
    for k, line in enumerate(data.decode().splitlines(), 1):
        m = LINE.match(line)
        if not m:
            sys.exit(f"{rel}, line {k}: not a line of certify.py: {line!r}")
        name = m.group(1)
        ok = (re.fullmatch(rf"{re.escape(stem)}_leaf\d+", name) if stem
              else name in ("E37_A6.cnf", "E37_A7.cnf", "E37_A9.cnf", "E37_A11.cnf"))
        if not ok:
            sys.exit(f"{rel}, line {k}: a line for {name}, which is not a leaf of this log")


def main():
    ap = argparse.ArgumentParser(description="Pack the logs of plan D into certificates/g13_chi_certlogs.tar.gz.")
    ap.add_argument("assembly", help="a merged copy of scripts/g13/chi with every share's logs")
    ap.add_argument("--out-dir", default=os.path.join(ROOT, "certificates"))
    a = ap.parse_args()
    for rel in BASE:                    # the assembly splits as the repository does
        mine, theirs = os.path.join(HERE, rel), os.path.join(a.assembly, rel)
        if not os.path.isfile(theirs) or open(theirs, "rb").read() != open(mine, "rb").read():
            sys.exit(f"{theirs} is missing or differs from {mine}")
    names, other = members(a.assembly)
    if "logs/E37_A.log" not in names or not any(n.endswith(".certlog") for n in names):
        sys.exit(f"{a.assembly}: no logs/E37_A.log or no cases/*.certlog: not a merged assembly")
    ignored = [n for n in other if not re.fullmatch(r"cases/.*\.cnf(\.tmp)?", n)]
    if ignored:
        print(f"not packed ({len(ignored)}): {', '.join(ignored[:10])}" + (" ..." if len(ignored) > 10 else ""))
    datas = {}
    for rel in names:
        with open(os.path.join(a.assembly, rel), "rb") as fh:
            datas[rel] = fh.read()
        if not rel.endswith(".icnf"):
            check_log(rel, datas[rel])
    os.makedirs(a.out_dir, exist_ok=True)
    out = os.path.join(a.out_dir, ARCHIVE)
    with open(out + ".tmp", "wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=9, mtime=0) as gz:
            with tarfile.open(fileobj=gz, mode="w", format=tarfile.GNU_FORMAT) as tar:
                for rel in names:
                    info = tarfile.TarInfo(rel)
                    info.size, info.mtime, info.mode, info.type = len(datas[rel]), MTIME, 0o644, tarfile.REGTYPE
                    info.uid = info.gid = 0
                    info.uname = info.gname = ""
                    tar.addfile(info, io.BytesIO(datas[rel]))
    os.replace(out + ".tmp", out)
    digest = sha256(open(out, "rb").read())
    with open(os.path.join(a.out_dir, SUMS), "w") as fh:
        fh.write(f"{digest}  {ARCHIVE}\n")
        for rel in names:
            fh.write(f"{sha256(datas[rel])}  {rel}\n")
    logs = sum(1 for n in names if not n.endswith(".icnf"))
    lines = sum(datas[n].count(b"\n") for n in names if not n.endswith(".icnf"))
    print(f"{out}: {len(names)} files ({logs} logs with {lines} lines, {len(names) - logs} cube files of re-splits), "
          f"{sum(map(len, datas.values())) / 1e6:.1f} MB unpacked, {os.path.getsize(out) / 1e6:.1f} MB; "
          f"sha256 {digest}")
    print(f"{os.path.join(a.out_dir, SUMS)}: the SHA-256 of the archive and of the {len(names)} files")


if __name__ == "__main__":
    main()
