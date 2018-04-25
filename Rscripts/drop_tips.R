# drop tip(s)

# dependencies
library(ape)
library(phytools)
library(dplyr)

# command line args
args = commandArgs(trailingOnly=TRUE)

tree_file <- args[1]
outgroup_file <- args[2]
out_file <- args[3]

# read tree
tree <- midpoint.root(read.tree(tree_file))

# read file containing tips to remove
outgroups <- read.delim(outgroup_file, sep = "\t", header=FALSE, stringsAsFactors = FALSE)

# remove outgroup(s)
outgroup.removed <- drop.tip(tree, tree$tip.label[tree$tip.label %in% outgroups$V1] )

# write tree
write.tree(outgroup.removed , file = out_file)
