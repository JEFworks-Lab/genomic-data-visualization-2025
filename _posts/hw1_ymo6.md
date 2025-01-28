---
layout: post
title:  "HW1: Gene expression pattern for GNB1 and HES4"
author: Alexandra Mo
jhed: ymo6
categories: [ HW1 ]
image: homework/hw1/hw1_ymo6.png
featured: false
---

### 1. What data types are you visualizing?

### 2. What data encodings (geometric primitives and visual channels) are you using to visualize these data types?

### 3. What about the data are you trying to make salient through this data visualization? 

### 4. What Gestalt principles or knowledge about perceptiveness of visual encodings are you using to accomplish this?

### 5. Code (paste your code in between the ``` symbols)

```r
library(ggplot2)

file <- "~/Downloads/eevee.csv.gz"
data <- read.csv(file)
data[1:5,1:10]


dim(data)
ncol(data)
nrow(data)
rownames(data)
colnames(data)


ggplot(data) +
geom_point(aes(x = aligned_x, y = aligned_y, col = HES4)) + 
scale_color_gradient(low = 'lightgrey', high = 'red')+
theme_bw()


ggplot(data) + geom_point(aes(x = aligned_x, y = aligned_y,col = GNB1)) + 
scale_color_gradient(low = 'lightgrey',high = 'green')+
theme_bw()

data$color <- ifelse(
  data$HES4 != 0 & data$GNB1 == 0, "red", ifelse(
    data$HES4 == 0 & data$GNB1 != 0, "green", ifelse(
      data$HES4 != 0 & data$GNB1 != 0, "brown",
      "lightgrey"
    )
  )
)

ggplot(data) + geom_point(aes(x = aligned_x, y = aligned_y, col = color)) + 
  scale_color_identity(
    guide = "legend",
    labels = c("HES4", "GNB1", "Both Genes", "Neither Gene"),
    breaks = c("red","green","brown","lightgrey")
  )  + 
  labs(title = "Gene expression pattern for GNB1 and HES4") + theme_bw()
```
