data <- read.csv('~/Desktop/genomic-data-visualization-2025/data/eevee.csv.gz', row.names=1)
pos <- data[,2:3] ## what happens if you try to deconvolve spatial positions?
colnames(pos) <- c('x', 'y')
gexp <- data[, 4:ncol(data)]
rownames(pos) <- rownames(gexp) <- data$barcode
head(pos)
gexp[1:5,1:5]

library(MERINGUE)
W <- MERINGUE::getSpatialNeighbors(pos, filterDist = 55)
MERINGUE::plotNetwork(pos, W)
class(W)
dim(W)
W[1:5,1:5]
## for grad students: try euclidean or KNN adjacency definition

g <- 'DSP'
g <- 'IGKC'
## leave it to you to normalize; think about how results may differ
gene <- gexp[, g]; names(gene) <- rownames(gexp)
set.seed(0) ## for reproducibility
gene <- sample(gene)
names(gene) <- rownames(gexp)
df <- data.frame(pos, g)
ggplot(df, aes(x=x,y=y,col=gene)) + geom_point(size=3) + scale_color_viridis_c()

## positive I, a negative I, or some other option?
?moranTest
moranTest(gene, W)
## shuffled is close to expected value

## try it out for yourself
## writing a loop to test all the genes
## see if you can find the gene with the most positive I
## gene with the most negative I
## recall ?sapply

results <- do.call(rbind, 
                   lapply(colnames(gexp), function(g) { 
                     gene <- gexp[, g]
                     names(gene) <- rownames(gexp)
                     moranTest(gene, W)
                     }
                     ))
rownames(results) <- colnames(gexp)
head(results)

g <- 'MMP11'
results[g,]
head(results[order(results[,1], decreasing=TRUE),])
gene <- gexp[, g]; names(gene) <- rownames(gexp)
df <- data.frame(pos, g)
ggplot(df, aes(x=x,y=y,col=gene)) + geom_point(size=3) + scale_color_viridis_c()

g <- 'ACTB'
results[g,]
head(results[order(results[,1], decreasing=TRUE),])
gene <- gexp[, g]; names(gene) <- rownames(gexp)
df <- data.frame(pos, g)
ggplot(df, aes(x=x,y=y,col=gene)) + geom_point(size=3) + scale_color_viridis_c()

## food for thought, 
## what if you did clustering or PCA or TSNE on just spatially variable genes?

g <- 'BMP15'
results[g,]
head(results[order(results[,1], decreasing=FALSE),])
gene <- gexp[, g]; names(gene) <- rownames(gexp)
df <- data.frame(pos, g)
ggplot(df, aes(x=x,y=y,col=gene)) + geom_point(size=3) + scale_color_viridis_c()

moranTest(gene, W, alternative='two.sided')

## given a uniformly expressed house keeping gene
## in a tissue where cell density is variable 
## we are using a spot-based detection approach
## what happens if we run Moran's I on the unnormalized data?
## what happens if we normalize for cell density? (equivalent to dividing by total genes detected)

## robustness in density using voronoi tesselation W, robust to density of pos
ppos <- pos
ppos[,1] <- ppos[,1]^10
plot(ppos)
df <- data.frame(ppos, g)
ggplot(df, aes(x=x,y=y,col=gene)) + geom_point(size=3) + scale_color_viridis_c()
## try it out for yourself
