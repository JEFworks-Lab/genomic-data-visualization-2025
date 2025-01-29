file <- '~/Desktop/genomic-data-visualization-2025/data/pikachu.csv.gz'
data <- read.csv(file)
data[1:5,1:10]

pos <- data[, 5:6]
rownames(pos) <- data$cell_id
gexp <- data[, 7:ncol(data)]
rownames(gexp) <- data$cell_id

## access the ith row
i = 1
gexp[i,]
gexp[c('6032', '6033'),] ## synatactic sugar
## access the jth column
j = 1
gexp[,j]
## regex is worth learning
gexp[, colnames(gexp)[grepl('^SLC', colnames(gexp))]]

## writing loops
x = Sys.time()
out <- c()
for(i in 1:10) {
  out <- c(out, i)
}
Sys.time() - x
out <- sapply(1:10, function(i) { 
  return(i)
})

mean(gexp[, 1])
var(gexp[, 1])

## loop through our gene expression matrix
## compute the variance across the first 100 genes
## output the result as a vector 
## repeat for the mean
## visualize the resulting relationship between mean and variance
## stretch goal for pros: use Sys.time() to compare runtimes

allvars1 <- sapply(gexp[,1:100], function(i) {
  var(i)
})
allvars2 <- sapply(1:100, function(i) {
  var(gexp[, i])
})
names(allvars2) <- colnames(gexp)[1:100]
max(allvars1)
allvars1[which(allvars1 == max(allvars1))]
sort(allvars1, decreasing=TRUE)[1]

allmeans <- sapply(1:ncol(gexp), function(i) {
  mean(gexp[, i])
})
names(allmeans) <- colnames(gexp)

## weird thing in R called do.call, row bind
df <- do.call(rbind, lapply(1:100, function(i) {
  return(data.frame(var = var(gexp[, i]), 
             mean = mean(gexp[, i]), 
             sum = sum(gexp[, i])))
}))
rownames(df) <- colnames(gexp)[1:100]
df

## want to explore relationship between mean and variance
library(ggplot2)
ggplot(df) + geom_point(aes(x = log10(var), y = log10(mean)))
ggplot(df) + geom_point(aes(x = var, y = mean)) + scale_x_log10() + scale_y_log10()

## geometric primitive of line
## y positional information to encode mean
## color hue encode the var
df$gene_name = rownames(df)
df = df[order(df$var, decreasing=TRUE),]
ggplot(df, aes(x = 1:100, y = mean, fill = var)) + 
  geom_bar(stat = "identity") 

#What is the max total genes are detected per cell/spot?
#What is the max unique genes are detected per cell/spot?
#What is the max total expression for a gene across all cells/spots?
#What is the max average expression for a gene across all cells/spots?








