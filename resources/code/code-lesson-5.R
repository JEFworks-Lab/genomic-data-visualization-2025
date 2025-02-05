file <- '~/Desktop/genomic-data-visualization-2025/data/pikachu.csv.gz'
data <- read.csv(file)
data[1:5,1:10]

pos <- data[, 5:6]
rownames(pos) <- data$cell_id
gexp <- data[, 7:ncol(data)]
rownames(gexp) <- data$barcode

## for eevee data, to make things faster, limit to top 1000 most highly expressed genes
topgenes <- names(sort(colSums(gexp), decreasing=TRUE)[1:1000])
gexpsub <- gexp[,topgenes]
dim(gexp)
dim(gexpsub)

?prcomp
pcs <- prcomp(gexpsub)
names(pcs)
pcs$sdev
pcs$rotation[1:5,1:5] ## beta values representing PCs as linear combinations fo genes
pcs$x[1:5,1:5]

## run PCA on your data
## goal 1: visualize a scree plot
df <- data.frame(sdev = pcs$sdev, index=1:length(pcs$sdev))
ggplot(df[1:10,], aes(x = index, y = sdev)) + geom_line()
## goal 2: explore the relationship between gene loadings, cells expression of genes, and cells PC eigenvector value
## goal 3: relationship between total gene expression and PCs
df <- data.frame(pcs$x, totalgene = rowSums(gexp), gene = gexp[, 'POSTN'])
ggplot(df, aes(x = PC1, y = PC2, col=totalgene)) + geom_point() +
  scale_color_gradient(low = 'lightgrey', high='red')
ggplot(df, aes(x = PC1, y = PC2, col=log10(gene+1))) + geom_point() +
  scale_color_gradient(low = 'lightgrey', high='red')
ggplot(df, aes(x = PC2, y = PC3, col=log10(gene+1))) + geom_point() +
  scale_color_gradient(low = 'lightgrey', high='red')
## goal 4: explore prcomp parameters with scale=TRUE

## normalize data
## 'CPM normalization'
## normalize to 10,000 (counts per 10,000)
## normalize to median 
hist(rowSums(gexp))
median(rowSums(gexp))
## need to divide each row by its row sum
norm <- gexp/rowSums(gexp) * 10000
rowSums(norm) ## sanity check

#pcs <- prcomp(norm[,topgenes])
pcs <- prcomp(norm)
df <- data.frame(sdev = pcs$sdev, index=1:length(pcs$sdev))
ggplot(df[1:100,], aes(x = index, y = sdev)) + geom_line()
df <- data.frame(pcs$x, totalgene = rowSums(gexp), gene = gexp[, 'POSTN'])
ggplot(df, aes(x = PC1, y = PC2, col=totalgene)) + geom_point() +
  scale_color_gradient(low = 'lightgrey', high='red')
ggplot(df, aes(x = PC1, y = PC2, col=log10(gene+1))) + geom_point() +
  scale_color_gradient(low = 'lightgrey', high='red')
## how does this all relate to physical space?


