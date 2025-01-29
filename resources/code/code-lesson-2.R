## pikachu - imaging
## eevee - sequencing
file <- '~/Desktop/genomic-data-visualization-2025/data/pikachu.csv.gz'
data <- read.csv(file)
data[1:5,1:10]

## how many genes are in your dataset?
## how many cells are in your dataset?
dim(data)
?ncol
?nrow
?rownames
?colnames

head(sort(colSums(data), decreasing=TRUE), n=20)
library(ggplot2)
## geometrive primitive = point
## xy = position
## saturation = expression magnitude ERBB2

ggplot(data) + geom_point(aes(x = aligned_x, y = aligned_y,
                              col = ERBB2)) +
  scale_color_gradient(low = 'lightgrey', high = 'red') +
  theme_bw()

## if you want to make more salient the fact that 
## the eevee dataset has positions that are crooked
## how would you do it?

ggplot(data) + geom_line(aes(x = POSTN, 
                              y = ERBB2,
                             col = cell_area)) +
  scale_color_gradient(low = 'lightgrey', high = 'red') +
  theme_bw()

ggplot(data) + geom_point(aes(x = aligned_x, 
                              y = aligned_y,
                              col = ERBB2,
                              shape = POSTN > 10),
                          size = 0.1) +
  scale_color_gradient(low = 'lightgrey', high = 'red') +
  theme_bw()
