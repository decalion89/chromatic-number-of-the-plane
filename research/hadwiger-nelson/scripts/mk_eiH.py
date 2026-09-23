import sys, json
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
exec(open("/home/user/darwin-50/research/hadwiger-nelson/scripts/screen_lambda.py").read().split("for name, pts in")[0])
ia, ib = VH.index(P(*extra[0])), VH.index(P(*extra[1]))
json.dump({"field_generators": [3, 11], "note": "Exoo-Ismailescu H (214 points), as a point set; A, B their forced pair at distance 5",
           "A": ia, "B": ib,
           "points": [[[[t.numerator, t.denominator] for t in q.x.c], [[t.numerator, t.denominator] for t in q.y.c]] for q in VH]},
          open("/home/user/darwin-50/research/hadwiger-nelson/data/ei_H214.json", "w"))
print("A, B =", ia, ib, "n =", len(VH))
