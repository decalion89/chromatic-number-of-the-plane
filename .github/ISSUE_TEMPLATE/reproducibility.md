---
name: Result does not reproduce
about: A test, certificate or command that fails or gives a different answer
title: "Reproduction: "
labels: reproducibility
---

**Command.** Exactly what you ran, from which directory.

**Expected.** What the documentation says should happen.

**Observed.** The output or error (the last lines are enough).

**Environment.** Operating system, `python3 --version`, and `python3 -m pip freeze | grep -i -E "numpy|scipy|sympy|mpmath|python-sat|cvxopt|clarabel"`. The versions used for the results are in `research/hadwiger-nelson/requirements-lock.txt`.
