library(ggplot2)
require(gridExtra)
library(dplyr)
library(AnomDetct)
library(parallel)
library(POT)
library(rmutil)

set.seed(112)
gen_gam <- rggamma(n = 6276, s = 1.5, m = 40, f = 1)
hist(gen_gam)

v1 <- runif(100, 65, 66)
v2 <- runif(50, 130, 131)
v3 <- runif(25, 195, 196)


dat_clust <- c(gen_gam, v1, v2, v3)
df <- as.data.frame(dat_clust)