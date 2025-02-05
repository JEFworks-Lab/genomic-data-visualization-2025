# file <- '~/Desktop/genomic-data-visualization-2025/data/eevee.csv.gz'
# data <- read.csv(file)
# data[1:5,1:10]
# 
# pos <- data[, 3:4]
# rownames(pos) <- data$barcode
# head(pos)
# gexp <- data[, 5:ncol(data)]
# rownames(gexp) <- data$barcode
# gexp[1:5,1:5]
# dim(gexp)

file <- '~/Desktop/genomic-data-visualization-2025/data/pikachu.csv.gz'
data <- read.csv(file)
data[1:5,1:10]

pos <- data[, 5:6]
rownames(pos) <- data$cell_id
gexp <- data[, 7:ncol(data)]
rownames(gexp) <- data$barcode

## limiting to top 1000 most highly expressed genes
## to explore: what happens if you use top 2000, or all genes?
#topgenes <- names(sort(colSums(gexp), decreasing=TRUE)[1:2000])
#gexpsub <- gexp[,topgenes]
#gexpsub[1:5,1:5]
#dim(gexpsub)
gexpsub <- gexp

## to explore: what happens if we normalize?
norm <- gexpsub/log10(data$cell_area+1) * 3
norm[1:5,1:5]
## what happens if we do a PCA, take the top N PCs, and then do a tSNE?
pcs <- prcomp(norm)
## what if you transpose this matrix t(gexpsub)

library(Rtsne)
#emb <- Rtsne(pcs$x[,1:5])
emb <- Rtsne(norm)
names(emb)
head(emb$Y)

library(ggplot2)
df <- data.frame(tSNE1 = emb$Y[,1], tSNE2 = emb$Y[,2],
                 gene = gexpsub[, 'ERBB2'], totgexp = rowSums(gexpsub), area = data$cell_area)
g3 <- ggplot(df, aes(x = tSNE1, y = tSNE2, col=log10(area+1))) + geom_point(size=0.01)
g3

library(patchwork)
g2 + g3

## how do different genes relate to each other?
## how do cells relate to each other in this lower dimensional representation
## how does that relate to space?

## idea: use color saturation to encode x and y spatial positions
df <- data.frame(tSNE1 = emb$Y[,1], tSNE2 = emb$Y[,2], x = data$aligned_x, y = data$aligned_y)
ggplot(df, aes(x = tSNE1, y = tSNE2, col=y)) + geom_point(size=0.01) + 
  scale_color_gradient(low = 'lightgrey', high='red')

ggplot(df, aes(x = x, y = y, col=tSNE2)) + geom_point(size=0.01) + 
  scale_color_gradient(low = 'lightgrey', high='red')

## idea: compare linear and non-linear dimensionality reduction
df <- data.frame(tSNE1 = emb$Y[,1], tSNE2 = emb$Y[,2],
                 PC1 = pcs$x[,1], PC2 = pcs$x[,2],
                 x = data$aligned_x, y = data$aligned_y)
head(df)
p1 <- ggplot(df, aes(x=x, y=y, col=tSNE1)) + geom_point(size=0.01)
p2 <- ggplot(df, aes(x=x, y=y, col=PC1)) + geom_point(size=0.01)
p3 <- ggplot(df, aes(x=tSNE1, y=PC1)) + geom_point(size=0.01)
p1+p2 + p3

## look at loading values on PC1
head(pcs$rotation[,1])
head(sort(pcs$rotation[,1], decreasing=TRUE))
head(sort(pcs$rotation[,1], decreasing=FALSE))
head(sort(pcs$rotation[,2], decreasing=TRUE))
head(sort(pcs$rotation[,2], decreasing=FALSE))

df <- data.frame(gene1 = norm[, 'POSTN'], gene2= norm[, 'KRT7'], gene3 = norm[,'ERBB2'],
                 tSNE1 = emb$Y[,1], tSNE2 = emb$Y[,2],
                 PC1 = pcs$x[,1], PC2 = pcs$x[,2])
ggplot(df, aes(x=gene1, gene3)) + geom_point()
ggplot(df, aes(x = tSNE1, y = tSNE2, col=log10(gene3+1))) + geom_point(size = 0.01)
ggplot(df, aes(x = PC1, y = PC2, col=log10(gene3+1))) + geom_point(size = 0.01)


