data <- read.csv('~/Desktop/genomic-data-visualization-2025/data/mystery_data.csv.gz', row.names=1)
data[1:5,1:5]

pos <- data[,1:2]
exp <- data[, 3:ncol(data)]
head(pos)
head(exp)

dim(exp)

norm <- log10(exp/rowSums(exp) * 1e6 + 1)
norm[1:5,1:5]

com <- kmeans(norm, centers=5)
clusters <- as.factor(com$cluster)
head(clusters)

library(ggplot2)
df <- data.frame(pos, clusters)
g1 <- ggplot(df, aes(x=x, y=y, col=clusters)) + geom_point()

ct1 <- names(clusters)[which(clusters == 1)]
ctother <- names(clusters)[which(clusters != 1)]
results <- sapply(colnames(norm), function(i) {
  wilcox.test(norm[ct1, i], norm[ctother, i])$p.value ## two sided test
})
names(results) <- colnames(norm)
sort(results[results < 0.05/ncol(norm)])

## Fam216a testis and cerebellum 
## Prkca involved in neuronal signal transduction

df <- data.frame(pos, gene=norm[,'Fam216a'])
g2 <- ggplot(df, aes(x=x, y=y, col=gene)) + geom_point()
library(patchwork)
g1 + g2

sort(colSums(norm))
## brain cytoplasm marker
df <- data.frame(pos, gene=norm[,'Bc1'])
ggplot(df, aes(x=x, y=y, col=gene)) + geom_point() + scale_color_viridis_c()

emb <- Rtsne::Rtsne(norm)
df <- data.frame(emb$Y, clusters)
ggplot(df, aes(x = X1, y=X2, col=clusters)) + geom_point()

com2 <- kmeans(emb$Y, centers=5)
clusters2 <- as.factor(com2$cluster)
names(clusters2) <- rownames(norm)
df <- data.frame(emb$Y, clusters2)
ggplot(df, aes(x = X1, y=X2, col=clusters2)) + geom_point()
df <- data.frame(pos, clusters2)
ggplot(df, aes(x = x, y=y, col=clusters2)) + geom_point()

elbow <- sapply(2:20, function(k) {
  out <- kmeans(norm, centers=k)
  out$tot.withinss
})
plot(2:20, elbow)

pcs <- prcomp(norm)
plot(pcs$sdev)
elbow <- sapply(2:20, function(k) {
  out <- kmeans(pcs$x[,1:2], centers=k)
  out$tot.withinss
})
plot(2:20, elbow)

com <- as.factor(kmeans(pcs$x[,1:2], centers=8)$cluster)
df <- data.frame(pcs$x[,1:2], com)
ggplot(df, aes(x = PC1, y=PC2, col=com)) + geom_point()
df <- data.frame(pos, com)
ggplot(df, aes(x = x, y=y, col=com)) + geom_point()

ct1 <- names(com)[which(com == 5)]
ctother <- names(com)[which(com != 5)]
results <- sapply(colnames(norm), function(i) {
  wilcox.test(norm[ct1, i], norm[ctother, i], alternative = 'greater')$p.value 
})
names(results) <- colnames(norm)
sort(results[results < 0.05/ncol(colnames)])

## Ttr - transport thyroid hormone
df <- data.frame(pcs$x[,1:2], gene = norm[,'Ttr'])
ggplot(df, aes(x = PC1, y=PC2, col=gene)) + geom_point()

## guess: neurons on dish 
## guess: organ on a chip
## guess: dissociated cells

## answer: simulated from a spot-based spatial transcriptomics dataset
## originally a brain tissue
## shuffled the genes across spots
## magnitude of genes is preserved, but the coherence is disrupted
## this data is rubish
## similar to real data from failed experiments, failed processing 
## conclusion: always be skeptical 
## conclusion: always visualize your data (hint)
exp[1:5,1:5]



