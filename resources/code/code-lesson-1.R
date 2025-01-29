## install ggplot if you have not already done so (https://ggplot2.tidyverse.org/)
#install.packages('ggplot2')
library(ggplot2)

data(iris) ## load built in dataset
class(iris)
iris ## would not recommend
head(iris) ## look at first few lines
iris[1:5,] ## first 5 rows
iris[, 1:5] ## first 5 columns
iris[1:5,1:5] ## top corner

iris$Species
table(iris$Species)

ggplot(iris) + 
  geom_point(aes(x = Sepal.Length, 
                 y = Petal.Length, 
                 shape = Species)) +
  theme_classic()

?geom_line

## goal for hands on component
## apply data visualization design principles learned today
## use continuity? 
## explore different visual channels?
## see if you can design a whole new data visualization? use lines? use areas? 
