from vbsupertree import *


newick_str = "((C:0.67,(A:0.85,B:0.66):1.00):0.76,(D:0.70,(F:0.93,E:0.63):0.87):0.94);"
tree = MyTree(newick_str)
print(tree)

big_sup = PCSPSupport.from_tree(tree)
print(big_sup)

abd_sup = big_sup.restrict(Clade("ABD"))
print(abd_sup)

cdef_sup = big_sup.restrict(Clade("CDEF"))
print(cdef_sup)

cdef_sup.mutualize(abd_sup)

