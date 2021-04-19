from vbsupertree import *

X = Clade("ABCDEF")
X1 = Clade("ABCDE")
X2 = Clade("ABCDF")

flat_sbn = SBN.random(X, 10)
trees_raw = [flat_sbn.random_tree() for _ in range(20)]
support = PCSPSupport.from_trees(trees_raw)

support1 = support.restrict(X1)
support2 = support.restrict(X2)

visited = set()
mutual_support = support1.mutualize(support2, verbose=True, visited=visited)

for a, b, c in sorted(visited, key=lambda x: (len(x[1]), x[1].clade()), reverse=True):
    print(f"{a}, {b}, {c}")

triplet_list = []
for a,b,c in visited:
    cl = a.clade
    W1 = Clade(X1 & cl)
    W2 = Clade(X2 & cl)
    b2 = SubsplitClade(b, W1)
    c2 = SubsplitClade(c, W2)
    triplet_list.append((a, b2, c2))

for a, b, c in sorted(triplet_list, key=lambda x: (len(x[1].clade), x[1].clade, x[1].other_clade, len(x[2].clade), x[2].clade, x[2].other_clade), reverse=True):
    print(f"{a}, {b}, {c}")

