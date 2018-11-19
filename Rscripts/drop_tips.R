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
tree <- read.tree(tree_file)

# read file containing tips to remove
outgroups <- read.delim(outgroup_file, sep = "\t", header=FALSE, stringsAsFactors = T)

# feedback
no_samples_b4 <- length(tree$tip.label)
print(paste(no_samples_b4,"tips in input tree"))

# find list of samples to remove
remove <- NULL
if (length(args)==3) {
  remove <- tree$tip.label[tree$tip.label %in% outgroups$V1]
}else{
  remove <- tree$tip.label[!(tree$tip.label %in% outgroups$V1)]
}

print(paste(length(remove),"samples to remove"))

# remove outgroup(s)
outgroup.removed <- tree
if ( no_samples_b4 == length(remove) ) {
  print(" ERROR - cannot remove all samples from the tree")
} else if ( 0 == length(remove) ) {
  print(" WARNING - no samples to remove")
} else {
  outgroup.removed <- drop.tip(tree, tree$tip.label[tree$tip.label %in% remove] )
}

# write tree
write.tree(outgroup.removed , file = out_file)
