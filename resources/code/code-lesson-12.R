## install.packages("gganimate")
## if you are having trouble rendering, try
## install.packages('gifski')
library(gganimate)

file <- '~/Desktop/genomic-data-visualization-2025/data/pikachu.csv.gz'
data <- read.csv(file)
data[1:5,1:10]

pos <- data[, 5:6]
rownames(pos) <- data$cell_id
gexp <- data[, 7:ncol(data)]
rownames(gexp) <- data$barcode
head(gexp)
head(pos)

loggexp <- log10(gexp+1)

com <- kmeans(loggexp, centers=7)
clusters <- com$cluster
clusters <- as.factor(clusters) 
names(clusters) <- rownames(gexp)
head(clusters)

pcs <- prcomp(loggexp)
emb <- Rtsne::Rtsne(loggexp)$Y
head(emb)

## make individual plots
library(ggplot2)
df1 <- data.frame(emb, clusters)
colnames(df1) <- c('x', 'y', 'clusters')
ggplot(df1, aes(x=x, y=y, col=clusters)) + geom_point(size=0.01)

df2 <- data.frame(pos, clusters)
colnames(df2) <- c('x', 'y', 'clusters')
ggplot(df2, aes(x = x, y = y, col=clusters)) + geom_point(size=0.01, alpha=0.5)

## in order to animate, need to make new data frame with all the information
df <- rbind(
  cbind(df1, order = 1),
  cbind(df2, order = 2)
)
head(df)
p <- ggplot(df, aes(x=x, y=y, col=clusters)) + geom_point(size = 0.01, alpha=0.5)

library(gganimate)
anim <- p + transition_states(order) + view_follow() + ease_aes('linear')
animate(anim, height=300, width=300) 
## https://gganimate.com/
## make your own, see if you can help improve this
## figure out how many ductal epithelial clusters there are?

df1 <- data.frame(emb, gene = loggexp[, 'CD4'])
df2 <- data.frame(emb, gene = loggexp[, 'POSTN'])
colnames(df2) <- colnames(df1) <- c('x', 'y', 'gene')
df <- rbind(
  cbind(df1, order = 1, name = 'CD4'),
  cbind(df2, order = 2, name = 'POSTN')
)
head(df)
p <- ggplot(df, aes(x=x, y=y, col=gene)) + geom_point(size = 0.01, alpha=0.5) +
  scale_color_viridis_c()

anim <- p + transition_states(order) + view_follow() + ease_aes('linear') +
  labs(title = '{closest_state}') 
animate(anim, height=300, width=300) 

## not effective use of animation
df1 <- data.frame(pos, col = clusters == 1)
df2 <- data.frame(pos, col = clusters == 2)
colnames(df2) <- colnames(df1) <- c('x', 'y', 'cluster')
df <- rbind(
  cbind(df1, order = 1),
  cbind(df2, order = 2)
)
head(df)
p <- ggplot(df, aes(x=x, y=y, col=cluster)) + geom_point(size = 0.01, alpha=0.5)
anim <- p + transition_states(order) + view_follow() + ease_aes('linear') +
  labs(title = '{closest_state}') 
animate(anim, height=300, width=300) 
