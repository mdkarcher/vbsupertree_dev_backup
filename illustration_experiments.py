from vbsupertree import *

X = Clade("ABCD")
flat_sbn = SBN.random(X, 10)
true_tree = flat_sbn.random_tree()
print(true_tree)
true_tree.tree.write(format=9)
# '(A,(B,(C,D)));'
true_tree = MyTree('(A,(B,(C,D)));')
# true_tree = MyTree('(((A,B),C),(D,E));')

X1 = Clade("ABD")
X2 = Clade("ACD")
tree1 = true_tree.restrict(X1)
tree2 = true_tree.restrict(X2)
print(tree1)
print(tree2)
true_tree.to_pcsp_support().to_string_set()
tree1.to_pcsp_support().to_string_set()
tree2.to_pcsp_support().to_string_set()

psup1 = PCSPSupport.from_tree(tree1)
psup2 = PCSPSupport.from_tree(tree2)
psupm = psup1.mutualize(psup2, verbose=True)

trees = list(psupm.all_trees(verbose=False))
len(trees)

for tree in trees:
    print(tree)


X1 = Clade("ABC")
X2 = Clade("BCD")

flat_sbn1 = SBN.random(X1, 10)
flat_sbn2 = SBN.random(X2, 10)

tree1 = flat_sbn1.random_tree()
print(tree1)
tree2 = flat_sbn2.random_tree()
print(tree2)

ref1 = tree1.to_pcsp_support()
ref2 = tree2.to_pcsp_support()

mut = ref1.mutualize(ref2)
len(list(mut.all_trees()))

