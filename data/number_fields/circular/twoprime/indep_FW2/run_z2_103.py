import test_finite_pieces as T
T.EXAMPLES = [("Z^2, S={(-1,2),(-1,1),(0,1),(3,0)} (10/3)", 2, 1, [(-1, 2, 0), (-1, 1, 0), (0, 1, 0), (3, 0, 0)],
               [[(0, n), (0, n)] for n in range(1, 14)] )]
res = T.run(cross_check_max=40)
print(res)
