file <- '~/Desktop/genomic-data-visualization-2025/data/pikachu.csv.gz'
data <- read.csv(file)
data[1:5,1:10]

pos <- data[, 5:6]
rownames(pos) <- data$cell_id
gexp <- data[, 7:ncol(data)]
rownames(gexp) <- data$barcode
head(gexp)
head(pos)

## if you normalize, if you log transform
## you can use different transformations of the gene expression
## for different parts of the analysis
loggexp <- log10(gexp+1)

com <- kmeans(loggexp, centers=5)
clusters <- com$cluster
clusters <- as.factor(clusters) ## tell R it's a categorical variable
names(clusters) <- rownames(gexp)
head(clusters)

pcs <- prcomp(loggexp)
df <- data.frame(pcs$x, clusters)
ggplot(df, aes(x=PC1, y=PC2, col=clusters)) + geom_point()

df <- data.frame(pcs$x, clusters, gene = gexp[, 'ERBB2'])
ggplot(df, aes(x=PC1, y=PC2, col=gene)) + geom_point()

## differential expression
?t.test
?wilcox.test

interest <- 5
cellsOfInterest <- names(clusters)[clusters == interest]
otherCells <- names(clusters)[clusters != interest]
i <- 'ERBB2'
genetest <- gexp[,i]
names(genetest) <- rownames(gexp)
genetest[cellsOfInterest]
genetest[otherCells]
t.test(genetest[cellsOfInterest], genetest[otherCells], alternative = 'two.sided')
t.test(genetest[cellsOfInterest], genetest[otherCells], alternative = 'greater')
t.test(genetest[cellsOfInterest], genetest[otherCells], alternative = 'less')

## same syntax
wilcox.test(genetest[cellsOfInterest], genetest[otherCells], alternative = 'two.sided')

results <- sapply(1:ncol(gexp), function(i) {
  genetest <- gexp[,i]
  names(genetest) <- rownames(gexp)
  out <- t.test(genetest[cellsOfInterest], genetest[otherCells], alternative = 'greater')
  out$p.value
})
names(results) <- colnames(gexp)

head(sort(results, decreasing = FALSE))
df <- data.frame(pcs$x, clusters, gene = loggexp[, 'BASP1'])
ggplot(df, aes(x=PC1, y=PC2, col=gene)) + geom_point()

## try t test on log transformed, and compare t test on linear
## try t test vs wilcox test, compare the p-values -> -log10(p values)
-log10(results)
## visualize on a tSNE instead of PCs

which(results < 0.05) ## technically should do multiple testing correction

## which of these clusters is a T cell cluster?
colnames(gexp)[grepl('^CD', colnames(gexp))]
i <- "CD3E"
genetest <- gexp[,i]
names(genetest) <- rownames(gexp)
results <- sapply(1:5, function(interest) {
  cellsOfInterest <- names(clusters)[clusters == interest]
  otherCells <- names(clusters)[clusters != interest]
  out <- t.test(genetest[cellsOfInterest], genetest[otherCells], alternative = 'greater')
  out$p.value
})
names(results) <- 1:5
results

## does this cluster 1 also upregulate CD4?
i <- "CD4"
genetest <- gexp[,i]
names(genetest) <- rownames(gexp)
results <- sapply(1:5, function(interest) {
  cellsOfInterest <- names(clusters)[clusters == interest]
  otherCells <- names(clusters)[clusters != interest]
  out <- t.test(genetest[cellsOfInterest], genetest[otherCells], alternative = 'greater')
  out$p.value
})
names(results) <- 1:5
results

emb <- Rtsne::Rtsne(pcs$x[,1:10])$Y
head(emb)
df <- data.frame(emb, clusters, gene = loggexp[, 'CD3E'])
ggplot(df, aes(x=X1, y=X2, col=gene)) + geom_point()
ggplot(df, aes(x=X1, y=X2, col=clusters)) + geom_point()

df <- data.frame(emb, clusters, gene = loggexp[, 'CD4'])
ggplot(df, aes(x=X1, y=X2, col=gene)) + geom_point()


