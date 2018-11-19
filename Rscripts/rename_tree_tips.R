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
colnames(outgroups) <- c("id", "replace")
outgroups$id <- as.factor(outgroups$id)

# join tip labels
meta <- data.frame(id = tree$tip.label) %>%
  left_join(outgroups)

# replace tip labels
tree$tip.label <- meta$replace

# write tree
write.tree(tree , file = out_file)
