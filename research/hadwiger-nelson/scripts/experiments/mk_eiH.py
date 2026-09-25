import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, json
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
exec(open(HN_DIR + "/scripts/screen_lambda.py").read().split("for name, pts in")[0])
ia, ib = VH.index(P(*extra[0])), VH.index(P(*extra[1]))
json.dump({"field_generators": [3, 11], "note": "Exoo-Ismailescu H (214 points), as a point set; A, B their forced pair at distance 5",
           "A": ia, "B": ib,
           "points": [[[[t.numerator, t.denominator] for t in q.x.c], [[t.numerator, t.denominator] for t in q.y.c]] for q in VH]},
          open(HN_DIR + "/data/ei_H214.json", "w"))
print("A, B =", ia, ib, "n =", len(VH))
