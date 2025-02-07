file <- '~/Desktop/genomic-data-visualization-2025/data/pikachu.csv.gz'
data <- read.csv(file)
data[1:5,1:10]

pos <- data[, 5:6]
rownames(pos) <- data$cell_id
gexp <- data[, 7:ncol(data)]
rownames(gexp) <- data$barcode

?kmeans
## try many ks
ks = c(5, 10, 15, 20)
totw <- sapply(ks, function(k) {
  print(k)
  com <- kmeans(gexp, centers=k)
  return(com$tot.withinss)
})
## how can I use this information?
plot(ks, totw) ## change this to ggplot

com <- kmeans(log10(gexp+1), centers=5)
clusters <- com$cluster
clusters <- as.factor(clusters) ## tell R it's a categorical variable
names(clusters) <- rownames(gexp)
head(clusters)

pcs <- prcomp(gexp)
df <- data.frame(pcs$x, clusters)
ggplot(df, aes(x=PC3, y=PC4, col=clusters)) + geom_point()

## what would have happened if I did kmeans clustering in PC-space?
com2 <- kmeans(pcs$x[,1:4], centers=5)
clusters2 <- as.factor(com2$cluster)
df <- data.frame(pcs$x, clusters=clusters2)
ggplot(df, aes(x=PC1, y=PC2, col=clusters)) + geom_point()

## what happens if I did kmeans on tSNE?
## what if I visualize the clusters in physical space?
## what happens if I normalize the data first?

emb <- Rtsne::Rtsne(gexp)
head(emb$Y)
df <- data.frame(emb$Y)
ggplot(df, aes(x=X1, y=X2)) + geom_point()
com3 <- kmeans(emb$Y, centers = 6)
clusters3 <- as.factor(com3$cluster)
df <- data.frame(emb$Y, clusters=clusters3)
ggplot(df, aes(x=X1, y=X2, col=clusters)) + geom_point(size=0.01)

df <- data.frame(emb$Y, clusters=clusters, totgexp = rowSums(gexp)) ## from gene expression
ggplot(df, aes(x=X1, y=X2, col=log10(totgexp+1))) + geom_point(size=0.01)

df <- data.frame(emb$Y, clusters=clusters)
ggplot(df, aes(x=X1, y=X2, col=clusters)) + geom_point(size=0.01)

## normalizing data? log transform data?


file <- '~/Desktop/genomic-data-visualization-2025/data/eevee.csv.gz'
data <- read.csv(file)
data[1:5,1:10]

pos <- data[, 3:4]
rownames(pos) <- data$barcode
head(pos)
gexp <- data[, 5:ncol(data)]
rownames(gexp) <- data$barcode
gexp[1:5,1:5]
dim(gexp)

## limiting to top 1000 most highly expressed genes
## to explore: what happens if you use top 2000, or all genes?
topgenes <- names(sort(colSums(gexp), decreasing=TRUE)[1:2000])
gexpsub <- gexp[,topgenes]
gexpsub[1:5,1:5]
dim(gexpsub)

## normalize by total expression
norm <- gexpsub/rowSums(gexpsub) * 10000
norm[1:5,1:5]
norm <- log10(norm + 1) ## does this matter?

## kmeans
ks = 1:25
totw <- sapply(ks, function(k) {
  print(k)
  com <- kmeans(norm, centers=k)
  return(com$tot.withinss)
})
## how can I use this information?
plot(ks, totw) ## change this to ggplot

com <- kmeans(norm, centers=15)
clusters <- as.factor(com$cluster)
head(clusters)

## dimensional reduction
#pcs <- prcomp(norm)
#df <- data.frame(pcs$x, clusters, pos)
emb <- Rtsne(norm)
df <- data.frame(emb$Y, clusters, pos)
ggplot(df, aes(x=X1, y=X2, col=clusters)) + geom_point()
ggplot(df, aes(x=aligned_x, y=aligned_y, col=clusters)) + geom_point(size=2)

## next week: differential expression analysis
## if we want to prove clusters are transcriptionally distinct populations of cells
## they must have some differentially expressed genes



