from vbsupertree import *

# X = Clade("ABCDEF")
# X1 = Clade("ABCDE")
# X2 = Clade("ABCDF")
# X1 = Clade("ABCD")
# X2 = Clade("ABEF")
# X1 = Clade("ABC")
# X2 = Clade("DEF")

# X = Clade("ABCDEFGH")
# X1 = Clade("EFGH")
# X2 = Clade("ABCD")

X = Clade("ABCDEFGHIJKL")
X1 = Clade("ABCDEKL")
X2 = Clade("FGHIJKL")

# flat_sbn = SBN.random(X, 10)
# trees_raw = [flat_sbn.random_tree() for _ in range(16)]
trees_raw = [MyTree.random(X) for _ in range(256)]
support = PCSPSupport.from_trees(trees_raw)

support1 = support.restrict(X1)
support2 = support.restrict(X2)

visited = set()
mutual_support = support1.mutualize(support2, verbose=True, visited=visited)

triplet_list = sorted(visited, key=lambda x: (len(x[1].subsplit.clade()), x[1].subsplit.clade()), reverse=True)
for a, b, c in triplet_list:
    print(f"{a}, {b}, {c}")

# triplet_list = []
# for a,b,c in visited:
#     cl = a.clade
#     W1 = Clade(X1 & cl)
#     W2 = Clade(X2 & cl)
#     b2 = SubsplitClade(b, W1)
#     c2 = SubsplitClade(c, W2)
#     triplet_list.append((a, b2, c2))

for a, b, c in sorted(triplet_list, key=lambda x: (len(x[1].clade), x[1].clade, x[1].other_clade, len(x[2].clade), x[2].clade, x[2].other_clade), reverse=True):
    print(f"{a}, {b}, {c}")

trip_dict = dict()
for a,b,c in triplet_list:
    bc = (b, c)
    if bc not in trip_dict:
        trip_dict[bc] = []
    trip_dict[bc].append(a)


for key, val in trip_dict.items():
    b, c = key
    print(f"* ({b},{c})")
    print(len(val))
    for item in val:
        print(f"  {item}")

for key, val in trip_dict.items():
    b, c = key
    if len(val) < 2:
        continue
    print(f"* ({b},{c}) has {len(val)}")
    for item in val:
        print(f"  {item}")

hits = {1: [], 2: [], "simultaneous": []}
for (sscl1, sscl2), sscl_list in trip_dict.items():
    print(f"Examining pair ({sscl1}, {sscl2}):")
    for big_sscl in sscl_list:
        print(f"  Examining sscl {big_sscl}")
        sscl_rest1 = big_sscl.restrict(X1)
        print(f"    X1 restriction: {sscl_rest1}")
        if sscl_rest1 == sscl1:
            print(f"      Same as X1 reference {sscl1}")
        elif sscl_rest1.clade == sscl1.clade and len(sscl_rest1.other_clade) == 0:
            print(f"      Trivial descendant of X1 reference {sscl1}")
        else:
            hits[1].append((sscl1, sscl2, sscl_list))
            print("      HIT!!!!!")
        sscl_rest2 = big_sscl.restrict(X2)
        print(f"    X2 restriction: {sscl_rest2}")
        if sscl_rest2 == sscl2:
            print(f"      Same as X2 reference {sscl2}")
        elif sscl_rest2.clade == sscl2.clade and len(sscl_rest2.other_clade) == 0:
            print(f"      Trivial descendant of X2 reference {sscl2}")
        else:
            hits[2].append((sscl1, sscl2, sscl_list))
            print("      HIT!!!!!")
        if sscl_rest1.clade == sscl1.clade and len(sscl_rest1.other_clade) == 0 and \
                sscl_rest2.clade == sscl2.clade and len(sscl_rest2.other_clade) == 0:
            print("      Simultaneous triviality!!!")
            hits["simultaneous"].append((sscl1, sscl2, sscl_list))

hits

len_dict = dict()
for key, val in trip_dict.items():
    b, c = key
    l = len(val)
    if l not in len_dict:
        len_dict[l] = []
    len_dict[l].append(f"({b},{c})")

for key, val in len_dict.items():
    print(f"* {key}")
    print(f"{val}")

len_dict

out = []
A = X1 - X2
B = X2 - X1
C = X1 & X2
for (sscl1, sscl2), sscl_list in trip_dict.items():
    n = len(sscl_list)
    cl1 = sscl1.clade
    ocl1 = sscl1.other_clade
    cl2 = sscl2.clade
    ocl2 = sscl2.other_clade
    a = len(cl1 & A) > 0
    b = len(cl2 & B) > 0
    c1 = len(cl1 & C) > 0
    c2 = len(cl2 & C) > 0
    astar = len(ocl1 & A) > 0
    bstar = len(ocl2 & B) > 0
    cstar1 = len(ocl1 & C) > 0
    cstar2 = len(ocl2 & C) > 0
    out.append((str(sscl1), str(sscl2), n, a, b, c1, c2, astar, bstar, cstar1, cstar2))

import csv

header = ["sscl1", "sscl2", "n", "a", "b", "c1", "c2", "astar", "bstar", "cstar1", "cstar2"]

with open("data/pcsp_mutualization_set_arith2.csv", 'w') as csv_file:
    csv_writer = csv.writer(csv_file)
    csv_writer.writerow(header)
    csv_writer.writerows(out)



for sscl, chilren in support.data.items():
    print(sscl)



# Time complexity experiments

# Scenario 1

from vbsupertree import *
import timeit

X = Clade("ABCDEF")
X1 = Clade("ABCDE")
X2 = Clade("ABCDF")

flat_sbn = SBN.random(X, 10)
trees_raw = [flat_sbn.random_tree() for _ in range(1)]
support = PCSPSupport.from_trees(trees_raw)
support1 = support.restrict(X1)
support2 = support.restrict(X2)

len(support1)
len(support2)

# %timeit support1.mutualize(support2)
timeit.timeit("support1.mutualize(support2)", number=100, globals=globals())/100

# 1 trees 4 x 4 = 16: 690 µs

# Automation

from vbsupertree import *
import timeit
import itertools


def leave_k_out(seq, k):
    if k < 0 or k > len(seq) // 2:
        raise ValueError("Argument k not appropriate length.")
    subset1 = seq[k:]
    subset2 = seq[:k] + seq[2*k:]
    return subset1, subset2


def complexity_experiment(taxa, k_out, n_trees, reps=5):
    print(f"Starting {taxa}, leave {k_out} out")
    X1, X2 = leave_k_out(taxa, k_out)
    trees_raw = [MyTree.random(taxa) for _ in range(n_trees)]
    print(f" Training trees ({n_trees}) sampled.")
    support = PCSPSupport.from_trees(trees_raw)
    support1 = support.restrict(X1)
    support2 = support.restrict(X2)
    size1 = len(support1)
    size2 = len(support2)
    print(f" {size1} * {size2} = {size1 * size2} in time: ", end="")
    time = timeit.timeit("support1.mutualize(support2)", number=reps, globals={"support1": support1, "support2": support2}) / reps
    print(time)
    mutual_support = support1.mutualize(support2)
    mut_size = len(mutual_support)
    return size1, size2, mut_size, time


def outer_loop(taxon_sets, ns, reps=5):
    result = []
    for taxon_set, n_trees in itertools.product(taxon_sets, ns):
        for k_out in range(1, (len(taxon_set) // 2) + 1):
            size1, size2, mut_size, time = complexity_experiment(taxon_set, k_out, n_trees, reps)
            result.append((taxon_set, k_out, n_trees, size1, size2, size1*size2, mut_size, time))
    return result


# complexity_experiment("ABCDEFGH", 3, 10, 25)

taxon_sets = ["ABCDEFGHIJ", "ABCDEFGHIJK", "ABCDEFGHIJKL", "ABCDEFGHIJKLM"]
ns = [1, 4, 9, 16, 25, 36, 49, 64]

results = outer_loop(taxon_sets, ns)


import csv

header = ["taxa", "k_out", "n_trees", "size1", "size2", "prod", "mut_size", "time"]

with open("data/pcsp_mutualization_runtimes3.csv", 'w') as csv_file:
    csv_writer = csv.writer(csv_file)
    csv_writer.writerow(header)
    csv_writer.writerows(results)

# Automation 2

from vbsupertree import *
import itertools


def leave_k_out(seq, k):
    if k < 0 or k > len(seq) // 2:
        raise ValueError("Argument k not appropriate length.")
    subset1 = seq[k:]
    subset2 = seq[:k] + seq[2*k:]
    return subset1, subset2


def count_complexity_experiment(taxa, k_out, n_trees):
    print(f"Starting {taxa}, leave {k_out} out")
    X1, X2 = leave_k_out(taxa, k_out)
    trees_raw = [MyTree.random(taxa) for _ in range(n_trees)]
    print(f" Training trees ({n_trees}) sampled.")
    support = PCSPSupport.from_trees(trees_raw)
    support1 = support.restrict(X1)
    support2 = support.restrict(X2)
    size1 = len(support1)
    size2 = len(support2)
    print(f" {size1} * {size2} = {size1 * size2}")
    stats = dict()
    mutual_support = support1.mutualize(support2, stats=stats)
    mut_size = len(mutual_support)
    visited_triplets=stats.get('visited_triplets', 0)
    skipped_due_to_visited = stats.get('skipped_due_to_visited', 0)
    potential_children_generated = stats.get('potential_children_generated', 0)
    pcsp_generated = stats.get('pcsp_generated', 0)
    childs_child_clades = stats.get('childs_child_clades', 0)
    print(f" visited_triplets={visited_triplets}")
    print(f" skipped_due_to_visited={skipped_due_to_visited}")
    print(f" potential_children_generated={potential_children_generated}")
    print(f" pcsp_generated={pcsp_generated}")
    print(f" childs_child_clades={childs_child_clades}")
    return size1, size2, mut_size, visited_triplets, skipped_due_to_visited, potential_children_generated, pcsp_generated, childs_child_clades


def count_outer_loop(taxon_sets, ns):
    result = []
    for taxon_set, n_trees in itertools.product(taxon_sets, ns):
        for k_out in range(1, (len(taxon_set) // 2) + 1):
            size1, size2, mut_size, visited_triplets, skipped_due_to_visited, potential_children_generated, pcsp_generated, childs_child_clades = count_complexity_experiment(taxon_set, k_out, n_trees)
            result.append((taxon_set, k_out, n_trees, size1, size2, size1*size2, mut_size, visited_triplets, skipped_due_to_visited, potential_children_generated, pcsp_generated, childs_child_clades))
    return result


taxon_sets = ["ABCDEFGHIJ", "ABCDEFGHIJK", "ABCDEFGHIJKL", "ABCDEFGHIJKLM"]
ns = [1, 4, 9, 16, 25, 36, 49, 64]

results = count_outer_loop(taxon_sets, ns)


import csv

header = ["taxa", "k_out", "n_trees", "size1", "size2", "prod", "mut_size", "visited_triplets", "skipped_due_to_visited", "potential_children_generated", "pcsp_generated", "childs_child_clades"]

with open("data/pcsp_mutualization_counts2.csv", 'w') as csv_file:
    csv_writer = csv.writer(csv_file)
    csv_writer.writerow(header)
    csv_writer.writerows(results)

