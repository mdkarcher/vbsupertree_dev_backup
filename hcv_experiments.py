import pickle

from vbsupertree import *
from parse_beast import *


# exploring HCV data

all_trees = parse_beast_nexus("data/HCV/hcv-1826.trees")
print(all_trees[-1])
tree_dist = TreeDistribution.from_list(all_trees[5001:])
len(tree_dist)
sbn = SBN.from_tree_distribution(tree_dist)
len(sbn)

# supertree HCV experiment

# all_trees = parse_beast_nexus("data/HCV/hcv.trees")
all_trees_1826 = parse_beast_nexus("data/HCV/hcv-1826.trees")
all_trees_1879 = parse_beast_nexus("data/HCV/hcv-1879.trees")
# tree_dist = TreeDistribution.from_list(all_trees[5001:])
tree_dist_1826 = TreeDistribution.from_list(all_trees_1826[5001:])
tree_dist_1879 = TreeDistribution.from_list(all_trees_1879[5001:])
# sbn = SBN.from_tree_distribution(tree_dist)
sbn_1826 = SBN.from_tree_distribution(tree_dist_1826)
sbn_1879 = SBN.from_tree_distribution(tree_dist_1879)
# support = sbn.support()
support_1826 = sbn_1826.support()
support_1879 = sbn_1879.support()
mutual_support = support_1826.mutualize(support_1879)

mutual_support.is_complete()

# ResolutionTrim HCV

tree_dict = parse_topology_count("data/HCV/ResolutionTrim/hcv_topology_count.txt", dust=1)
tree_dict_1844 = parse_topology_count("data/HCV/ResolutionTrim/hcv-1844_topology_count.txt", dust=1)
tree_dict_1879 = parse_topology_count("data/HCV/ResolutionTrim/hcv-1879_topology_count.txt", dust=1)
tree_dist = TreeDistribution(tree_dict)
tree_dist.normalize()
tree_dist_1844 = TreeDistribution(tree_dict_1844)
tree_dist_1844.normalize()
tree_dist_1879 = TreeDistribution(tree_dict_1879)
tree_dist_1879.normalize()
sbn = SBN.from_tree_distribution(tree_dist)
sbn_1844 = SBN.from_tree_distribution(tree_dist_1844)
sbn_1879 = SBN.from_tree_distribution(tree_dist_1879)
support = sbn.support()
support_1844 = sbn_1844.support()
support_1879 = sbn_1879.support()
pcsp_probabilities = sbn.pcsp_probabilities()
pcsp_probabilities_1844 = sbn_1844.pcsp_probabilities()
pcsp_probabilities_1879 = sbn_1879.pcsp_probabilities()
mutual_support = support_1844.mutualize(support_1879)
# with open("data/HCV/ResolutionTrim/saved_states/mutual_support.pkl", "wb") as ff:
#     pickle.dump(mutual_support, ff)
# with open("data/HCV/ResolutionTrim/saved_states/mutual_support.pkl", "rb") as ff:
#     foo = pickle.load(ff)
mutual_support.is_complete()
mutual_support = mutual_support.prune(verbose=True)

uncovered_support_1844 = support_1844.to_set() - mutual_support.restrict(sbn_1844.root_clade()).to_set()
uncovered_support_1879 = support_1879.to_set() - mutual_support.restrict(sbn_1879.root_clade()).to_set()
len(uncovered_support_1844)
len(uncovered_support_1879)

total_uncovered_probabilities_1844 = sum(pcsp_probabilities_1844[pcsp] for pcsp in uncovered_support_1844)
total_uncovered_probabilities_1879 = sum(pcsp_probabilities_1879[pcsp] for pcsp in uncovered_support_1879)
print(total_uncovered_probabilities_1844)
print(total_uncovered_probabilities_1879)

trimmed_sbn_1844 = sbn_1844.copy()
trimmed_sbn_1844.remove_many(uncovered_support_1844)
trimmed_sbn_1844 = trimmed_sbn_1844.prune()
trimmed_sbn_1844.normalize()

trimmed_sbn_1879 = sbn_1879.copy()
trimmed_sbn_1879.remove_many(uncovered_support_1879)
trimmed_sbn_1879 = trimmed_sbn_1879.prune()
trimmed_sbn_1879.normalize()

starting_sbn = SBN.random_from_support(support=mutual_support, concentration=10)
# len(starting_sbn)
# len(mutual_support)
# sum(pcsp_probabilities.get(pcsp, 0.0) for pcsp in mutual_support.to_set() - starting_sbn.support().to_set())

true_sbn_trim = sbn.copy()
tst_support = true_sbn_trim.support()
tst_uncovered = tst_support.to_set() - starting_sbn.support().to_set()
sum(pcsp_probabilities[pcsp] for pcsp in tst_uncovered)
# 0.1371
true_sbn_trim.remove_many(tst_uncovered)
true_sbn_trim.is_complete()
true_sbn_trim = true_sbn_trim.prune()
true_sbn_trim.normalize()
true_sbn_trim.support().to_set() - starting_sbn.support().to_set()

true_sbn_trim.kl_divergence(starting_sbn)

supertree_sbn, kl_list, true_kl_list = starting_sbn.gradient_descent(
    references=[trimmed_sbn_1844, trimmed_sbn_1879],
    starting_gamma=2.0, max_iteration=50, true_reference=true_sbn_trim
)

# plotting

from matplotlib import pyplot as plt
import seaborn as sns


plt.rcParams.update({'font.size': 100})

sns.set_context("poster")
fig, (ax_kl, ax_true_kl) = plt.subplots(1, 2, figsize=(10, 5), sharey='none', sharex='none', constrained_layout=True)
ax_kl.set(title="Supertree Loss", xlabel="Iteration", ylabel="Nats", yscale="log")
ax_kl.plot(kl_list)
ax_true_kl.set(title="KL vs. Truth", xlabel="Iteration", ylabel="", yscale="log")
ax_true_kl.plot(true_kl_list)
plt.savefig(f"figures/hcv/vbsupertree_hcv_1844_1879_horiz.pdf", format="pdf")

sns.set_context("poster")
fig, (ax_kl, ax_true_kl) = plt.subplots(2, 1, figsize=(5, 8), sharey='none', sharex='none', constrained_layout=True)
ax_kl.set(ylabel="Supertree Loss", yscale="log")
ax_kl.plot(kl_list)
ax_true_kl.set(xlabel="Iteration", ylabel="KL vs. Truth", yscale="log")
ax_true_kl.plot(true_kl_list)
plt.savefig(f"figures/hcv/vbsupertree_hcv_1844_1879_vert.pdf", format="pdf")

