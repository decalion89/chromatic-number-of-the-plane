#!/bin/sh
# Rebuild OpenAI's transfer for a set of forbidden distances (OAI/Ext) on top of their library, and print the axioms.
#
# usage: sh build_distances.sh MATH_LEAN_DIR LEAN_BIN LEAN_PATH WORK_DIR [JOBS]
#
# MATH_LEAN_DIR  the lean/ directory of a checkout of https://github.com/openai/math at commit
#                fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb (it holds OAI/...); it is not modified
# LEAN_BIN       the lean binary of toolchain v4.34.1
# LEAN_PATH      the Mathlib build, as for ../build_closure.py
# WORK_DIR       a new directory: the sources go to WORK_DIR/src, the .olean files to WORK_DIR/build
# JOBS           parallel compilations (default 2; one compilation can take 6 GB of memory)
#
# OpenAI's 33 modules in the import closure are compiled from their unmodified sources, then the four files of
# OAI/Ext/ here. The last lines printed are the output of `#print axioms` in OAI/Ext/Check.lean.
set -eu
here=$(cd "$(dirname "$0")" && pwd)
math=$1; lean=$2; lp=$3; work=$4; jobs=${5:-2}
if [ -e "$work" ]; then echo "$work exists; give a new directory" >&2; exit 2; fi
if [ -e "$math/OAI/Ext" ]; then echo "$math/OAI/Ext exists; give a clean checkout of openai/math" >&2; exit 2; fi
mkdir -p "$work/src"
cp -R "$math/OAI" "$work/src/"
mkdir -p "$work/src/OAI/Ext"
cp "$here"/OAI/Ext/*.lean "$work/src/OAI/Ext/"
python3 "$here/../build_closure.py" "$work/src" "$lean" "$lp" "$work/build" "$jobs" OAI.Ext.Check
cat "$work/build/OAI/Ext/Check.out"
