
# split input args 
args = commandArgs(trailingOnly=TRUE)

# files
gwas_file <- args[1]
#gwas_file <- "/home/sbayliss/Desktop/projects/multispecies/pyseer/lmm/c_coli/porcine/1.results.tsv"

# data
data <- read.delim(gwas_file, header = T, sep = "\t")

# function https://www.gettinggeneticsdone.com/2010/07/qq-plots-of-p-values-in-r-using-base.html
ggd.qqplot = function(pvector, main=NULL, ...) {
  o = -log10(sort(pvector,decreasing=F))
  e = -log10( 1:length(o)/length(o) )
  plot(e,o,pch=19,cex=1, main=main, ...,
       xlab=expression(Expected~~-log[10](italic(p))),
       ylab=expression(Observed~~-log[10](italic(p))),
       xlim=c(0,max(e)), ylim=c(0,max(o)))
  lines(e,e,col="red")
}

# plot and save
png(sprintf("%s.png", gwas_file), width = 960, height = 960)
ggd.qqplot(data$lrt.pvalue)
dev.off()

