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

#tree_file <- "/home/sbayliss/Desktop/projects/Bayliss_Neisseria/phylogeny/CC23/CFML/embranch_masked.rooted_on_outgroup.contam_removed.tre"
#outgroup_file <- "/home/sbayliss/Desktop/projects/Bayliss_Neisseria/phylogeny/CC23/CFML/CC23_samples_rename.txt"

# read tree
tree <- read.tree(tree_file)

# read file containing tips to remove
outgroups <- read.delim(outgroup_file, sep = "\t", header=FALSE, stringsAsFactors = T)
colnames(outgroups) <- c("id", "replace")
outgroups$id <- as.factor(outgroups$id)

# join tip labels
meta <- data.frame(id = tree$tip.label) %>%
  left_join(outgroups)

# replace tip labels
tree$tip.label <- meta$replace

# write tree
write.tree(tree , file = out_file)
