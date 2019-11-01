
# plot gwas hits for multiple variants in one figure
library(grid)
library(scales)
library(ggplot2)
library(dplyr)
library(phangorn)
library(ggtree)
library(ggpubr)
library(gridExtra)

# functions
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

pallete <- c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c")

# output 
output = "/home/sbayliss/Desktop/projects/Pascoe_IBS/treeWAS/results_all_vars/"

# files
hit_file <- "/home/sbayliss/Desktop/projects/Pascoe_IBS/treeWAS/results_all_vars/combined.all_loci.annotated.freq.tab"
pangenome_file<- "/home/sbayliss/Desktop/projects/Pascoe_IBS/PIRATE/PIRATE.gene_families.ordered.tsv"
tree_file <- "/home/sbayliss/Desktop/projects/Pascoe_IBS/tree/RAxML_bestTree.IBS_maskedML.midpoint.nwk"
var_file <- "/home/sbayliss/Desktop/projects/Pascoe_IBS/treeWAS/inputs_all/combined_variants.0.01.tab"
metadata_file <- "/home/sbayliss/Desktop/projects/Pascoe_IBS/metadata/all/phenotypes.txt"

# read data 
pangenome <- read.delim(pangenome_file, sep = "\t", header = T, comment.char = '#')
hits_raw <- read.delim(hit_file, sep = "\t", header = T, comment.char = '#')
colnames(hits_raw)[1] <- "variant"
variants <- read.delim(var_file, sep = "\t", header = T, comment.char = '#')
tree_raw <- read.tree(tree_file)
meta <- read.delim(metadata_file, sep = "\t", header = T, comment.char = '#')

# get pangenome order
order <- as.factor(pangenome$gene_family)

# combine hits and variants 
hits <- hits_raw %>% left_join(variants)

# check intersection between tree and variants 
samples <- gsub("^X", "", colnames(hits)[24:length(colnames(hits))])
tips <- tree_raw$tip.label
remove <- tips[!(tips %in% intersect(tips, samples))]

# prepare tree
tree_raw$tip.label <- gsub(pattern = ":", replacement = "", tree_raw$tip.label)
tree <- tree_raw
if (length(remove)>0){
  tree <- drop.tip(tree_raw, remove)
}

# make manhattan plot 

# select type and index
type <- "subsequent" # "terminal", "simultaneous", "subsequent"
counti <- 1:length(colnames(hits))
ind <-  counti[colnames(hits) == type ]

# make test varianble
hits$test <- hits[,ind]

# get hit per gene - sort by min p
uniq <- hits %>% group_by(gene_family) %>% slice(which.min(test)) %>% arrange(test)

# filter hits for SNPS and core
hits_sub <- hits %>% filter(variant_type %in% c("core_allele", "snp"))

# filter order on genes present in hits 
order_filt <- as.factor(pangenome$gene_family[pangenome$gene_family %in% hits_sub$gene_family])

# plot as manhattan plot
mh <- ggplot(hits_sub, aes(x=gene_family,y=test) ) + 
  geom_point(colour = "black") + # aes(colour=variant_type)) + 
  scale_y_continuous("GWAS p-value", trans=reverselog_trans(10)) +
  scale_x_discrete (limits=order_filt) +
  scale_color_manual(values = pallete) +
  theme(legend.position = "none", 
        axis.text.x = element_blank()) +
  xlab("Position in Genome") 
mh

## add frequency plot

# print 
ggsave (plot=mh, filename = sprintf("%s/%s.manhattan.tiff", output, type), device = "tiff", height = 7, width = 12) 


# plot tree and top variants

# top X hits
top_n <- 10
top <- uniq[1:top_n,]

# plot tree with label 
t_A <- ggtree(tree) %<+% meta + 
  geom_tiplab(align = T) +
  geom_tippoint(aes(colour = as.factor(IBS))) +
  scale_color_manual(values = c("0"="blue", "1"="red"))
t_A

# select top hist and make hm table 
l <- length(top[1,])
hmp <- t(top[,24:l])
colnames(hmp) <- top$variant
rownames(hmp) <- gsub("^X", "", rownames(hmp))
hmp[hmp==1] <- "Present"
hmp[hmp==0] <- "Absent"
hmp <- as.data.frame(hmp)
hm <- NULL
hm <- gheatmap(t_A, hmp, colnames_position = NULL, colnames = T, colnames_angle = 90, 
               colnames_offset_y = -5, color = "black", font.size = 2, offset = 0.0002) +
  scale_fill_manual(values = c("Present"="firebrick", "Absent"="lightgrey") ) + 
  theme(plot.margin=unit(c(0.1,0.1,1,0.8),"cm")) + theme(legend.position = 0)
hm

# print 
ggsave (hm, filename = sprintf("%s/%s.heatmap.tiff", output, type), device = "tiff", height = 14, width = 12) 

