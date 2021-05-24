setwd("~/Documents/Research/Matsen/Experimental/vbsupertree")

library(tidyverse)
library(readr)
library(broom)
library(rpart)
library(rpart.plot)

# runtimes <- read_csv("data/pcsp_mutualization_runtimes.csv", 
                     # col_types = cols(k_out = col_integer(),
                                      # n_trees = col_integer(), size1 = col_integer(), 
                                      # size2 = col_integer(), prod = col_integer()))
runtimes <- read_csv("data/pcsp_mutualization_runtimes3.csv", 
                     col_types = cols(k_out = col_integer(), 
                                      n_trees = col_integer(), size1 = col_integer(), 
                                      size2 = col_integer(), prod = col_integer(), 
                                      mut_size = col_integer()))
View(runtimes)

subset = runtimes %>% filter(taxa == "ABCDEFGHIJK", k_out == 4)
plot(log(time) ~ log(prod), data=subset)
mod = lm(log(time) ~ log(prod), data=subset)
mod

results = runtimes %>% filter(n_trees > 4) %>% group_by(taxa, k_out) %>% do(tidy(lm(log(time) ~ log(prod), data=.)))
View(results)

results2 = runtimes %>% group_by(taxa, k_out) %>% do(tidy(lm(log(time) ~ log(mut_size), data=.)))
View(results2)


counts <- read_csv("data/pcsp_mutualization_counts2.csv", 
                   col_types = cols(k_out = col_integer(), 
                                    n_trees = col_integer(), size1 = col_integer(), 
                                    size2 = col_integer(), prod = col_integer(), 
                                    mut_size = col_integer(), visited_triplets = col_integer(), 
                                    skipped_due_to_visited = col_integer(), 
                                    potential_children_generated = col_integer(),
                                    pcsp_generated = col_integer(), childs_child_clades = col_integer()))
View(counts)

subset = counts %>% filter(taxa == "ABCDEFGHIJKL", k_out == 4)
plot(log(pcsp_generated) ~ log(prod), data=subset)

# results = counts %>% filter(n_trees > 0) %>% group_by(taxa, k_out) %>% do(tidy(lm(log(pcsp_generated) ~ log(prod), data=.))) %>% filter(term == "log(prod)")
results = counts %>% filter(n_trees > 0) %>% group_by(taxa, k_out) %>% do(tidy(lm(log(potential_children_generated) ~ log(prod), data=.))) %>% filter(term == "log(prod)")
View(results)

sum(counts["childs_child_clades"]/counts["pcsp_generated"] > 2)


set_arith <- read_csv("data/pcsp_mutualization_set_arith2.csv", 
                      col_types = cols(n = col_character()))
View(set_arith)

threes = set_arith %>% filter(n == "3")
twos = set_arith %>% filter(n == "2")
ones = set_arith %>% filter(n == "1")

set_arith %>% group_by(n) %>% summarize_at(vars(a:cstar2), mean)

set_arith %>% filter(n == "3", c1 == FALSE)
set_arith %>% filter(!a & !b)

rtree = rpart(n ~ a + b + c1 + astar + bstar + cstar1 + cstar2, data=set_arith)
rpart.plot(rtree)
