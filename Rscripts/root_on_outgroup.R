# root on outgroup(s)

# dependencies
library(ape)
library(phytools)
library(dplyr)

# subfunctions
oldest.mrca<-function(tree,tips){
  H<-nodeHeights(tree)
  X<-mrca(tree)
  n<-length(tips)
  nodes<-height<-vector(); k<-1
  for(i in 1:(n-1)) for(j in (i+1):n){
    nodes[k]<-X[tips[i],tips[j]]
    height[k]<-H[match(nodes[k],tree$edge[,1]),1]
    k<-k+1
  }
  z<-match(min(height),height)
  return(nodes[z])
}

# command line args

args = commandArgs(trailingOnly=TRUE)

tree_file <- args[1]
outgroup_file <- args[2]
out_file <- args[3]

# read tree
tree <- midpoint.root(read.tree(tree_file))

# read outgroup file
outgroups <- read.delim(outgroup_file, sep = "\t", header=FALSE, stringsAsFactors = FALSE)

# root on outgroup(s)
rooted.tree <- tree
if(length(outgroups$V1)>1){
  
  # multiple - find outgroup deepest node to root on
  outgroup.node <- oldest.mrca( tree, tree$tip.label[tree$tip.label %in% outgroups$V1] )
  rooted.tree <- root(tree, node=outgroup.node, resolve.root = TRUE)

}else{
  
  # single - root on outgroup
  rooted.tree <- root(tree, as.character(outgroups$V1[1]), resolve.root = TRUE)

}

# optional: remove outgroups
if(length(args)>3){
  rooted.tree <- drop.tip(rooted.tree, tree$tip.label[tree$tip.label %in% outgroups$V1], trim.internal = TRUE, rooted=1)
}

# write tree
ape::write.tree(rooted.tree , file = out_file)
