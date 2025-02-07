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

# compute overall dataset statistics
num_rows <- nrow(data)
num_cols <- ncol(data)
column_names <- colnames(data)
data_types <- sapply(data, class)

print(paste("Number of Rows:", num_rows))
print(paste("Number of Columns:", num_cols))
print("Column Names:")
print(column_names)
print("Data Types:")
print(data_types)

#the rows are individual cells it seems and the columns are a mix of mainly gene expression (integers) and statistics like cell area, x, y coordinates too


# dim(data)
# ncol(data)

# head(sort(colSums(data), decreasing=TRUE), n=20)
library(ggplot2)
# ggplot(data) + geom_point(aes(x=POSTN, y = LUM, col=ERBB2)) + scale_color_gradient(low='lightgrey', high = 'red')

# ggplot(data) + geom_point(aes(x=aligned_x, y = aligned_y, col=cell_area)) + scale_color_gradient(low='lightgrey', high = 'red') + theme_bw()

ggplot(data) +
  geom_point(aes(x = aligned_x, y = aligned_y, col = cell_area)) +
  scale_color_gradient(low = 'lightgrey', high = 'red') +
  stat_density_2d(aes(x = aligned_x, y = aligned_y, fill = ..density..), geom = "raster", contour = FALSE, alpha = 0.5) +
  theme_bw() +
  labs(title = "Cell Size and Density Clustering",
       x = "Aligned X",
       y = "Aligned Y",
       color = "Cell Area",
       fill = "Density")


cor_x <- cor(data$cell_area, data$aligned_x, use = "complete.obs")
cor_y <- cor(data$cell_area, data$aligned_y, use = "complete.obs")

print(paste("Correlation between cell area and X coordinate:", round(cor_x, 3)))
print(paste("Correlation between cell area and Y coordinate:", round(cor_y, 3)))

#a next step I would like to take is to relate the cell size to the density more clearly through partitioning out sections of the x,y grid, 
#assigning densities to those regions and then assigning the region's density to the cells inside the region. Therefore, each cell will have two attributes: 
# a density and a cell size. I can then do a simple scatter from there to visually inspect if there is a correlation. 

