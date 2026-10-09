"""A machine-checkable attack on the Hadwiger-Nelson problem.

chi(R^2), the least number of colours needed to paint the plane so that no two
points at distance exactly 1 share a colour, is 6 or 7: OpenAI proved chi >= 6 in
September 2026, without exhibiting a graph.  The lower bound moved from 4 to 5 in
2018 (de Grey), by exhibiting a finite unit-distance graph with no proper
4-colouring.  An explicit witness for 6 is a finite unit-distance graph with no
proper 5-colouring -- a finite object, and therefore searchable.

This package builds such graphs with exact arithmetic and decides their
colourability with a SAT solver, emitting certificates a third party can check
without trusting any of this code.
"""

__version__ = "1.0.0"
