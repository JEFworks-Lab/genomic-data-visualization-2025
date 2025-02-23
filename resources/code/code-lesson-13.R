data <- read.csv('~/Desktop/genomic-data-visualization-2025/data/eevee.csv.gz', row.names=1)
head(data)
pos <- data[,2:3] ## what happens if you try to deconvolve spatial positions?
colnames(pos) <- c('x', 'y')
gexp <- data[, 4:ncol(data)]
rownames(pos) <- rownames(gexp) <- data$barcode
head(pos)
gexp[1:5,1:5]

## you may need to install fgsea from Bioconductor's BiocManager and remotes
# install.packages("package") ## installing from CRAN
# BiocManager::install("fgsea") ## install from Bioconductor 
# require(remotes) 
# remotes::install_github('JEFworks-Lab/STdeconvolve') ## installing from Github
library(STdeconvolve)
## remove pixels with too few genes
counts <- cleanCounts(t(gexp), min.lib.size = 100, verbose=TRUE, plot=TRUE)
## feature select for genes
corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05, nTopOD=1000) ## changing parameters?
## choose optimal number of cell-types
ldas <- fitLDA(t(as.matrix(corpus)), Ks = seq(5, 8)) ## pick a different K from me
## get best model results
optLDA <- optimalModel(models = ldas, opt = "8")
## extract deconvolved cell-type proportions (theta) and transcriptional profiles (beta)
results <- getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1000)
deconProp <- results$theta
deconGexp <- results$beta
## visualize deconvolved cell-type proportions
vizAllTopics(deconProp, pos,
             r=8, lwd=0)	 

## run through this yourself, with a different K
## see if you can figure out where are B cells? (MS4A1+)




