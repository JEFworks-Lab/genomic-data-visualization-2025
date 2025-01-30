# install.packages('ggplot2')
# library(ggplot2)
# 
# data(iris)
# class(iris)
# iris 
# head(iris)
# 
# iris[1:5,]
# iris[,1:5]
# iris[1:5,1:5]
# 
# iris$Species
# table(iris$Species)
# 
# ggplot(data = iris, aes(x = Sepal.Length, y = Petal.Width, shape = Species)) +
#   geom_point(size = 3) +
#   labs(title = "Sepal Dimensions by Species",
#        x = "Sepal Length",
#        y = "Sepal Width") +
#   theme_minimal()

file <- '~/code/genomic-data-visualization-2025/data/pikachu.csv.gz'
data <- read.csv(file)
data[1:5, 1:20]

dim(data)
ncol(data)

#relationship of two genes
# head(sort(colSums(data), decreasing=TRUE), n=20)
library(ggplot2)
ggplot(data) + geom_point(aes(x=POSTN, y = LUM, col=ERBB2)) + scale_color_gradient(low='lightgrey', high = 'red')

ggplot(data) + geom_point(aes(x=aligned_x, y = aligned_y, col=ERBB2)) + scale_color_gradient(low='black', high = 'red') + theme_bw()

#what is postn? erbb2