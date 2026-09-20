"""The threshold, and the orbit-wise reduction.

N = {x : x.sigma(x) = 1} factors along the sigma-orbits of the primes above 5,
because sigma permutes the factors of A = O/5 within orbits and the condition
is componentwise.  So a functional supported on ONE orbit and zero elsewhere
is nonzero on all of N as soon as it is nonzero on that orbit's N_i:

    IF ONE ORBIT ADMITS AN EVERYWHERE-NONZERO FUNCTIONAL, THE WHOLE FIELD DOES.

Blocking therefore needs EVERY orbit to block on its own, and the orbit types
are indexed by one number: the residue degree f of the prime, with sigma
either fixing it (f even, sigma = Frob^{f/2}) or swapping a pair.  Degrees 2,
4 and 6 said: fixed f = 6 and swapped f = 3 block, nothing smaller does.  Both
of those are residue degree 3 over the REAL subfield.  This checks degree 8 --
fixed f = 8 and swapped f = 4, i.e. residue degree 4 over F -- and the mixed
type f = 6 + f = 2, which the reduction above predicts must fail.
"""
import sys
sys.path.insert(0, "/tmp/claude-0/-home-user-darwin-50/"
                   "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad")
from fieldtypes import decide

print("predicted to block (residue degree 4 over the real subfield):",
      flush=True)
decide([("fix", 8)])
decide([("swap", 4)])
print("predicted to FAIL by the orbit-wise reduction, despite containing a "
      "blocking orbit:", flush=True)
decide([("fix", 6), ("fix", 2)])
decide([("swap", 3), ("swap", 1)])
