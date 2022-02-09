### AE Melton, 2020
# PCA stuff

#
library(ggfortify)
data(iris)
head(iris)
#

# First step of a PCA: standardize the data
iris.scaled <- scale(x = iris[,1:4], center = T, scale = T)
head(iris.scaled)
#

# Second step: covariance matrix
cov.mat <- cov(x = iris.scaled)
cov.mat
#

# Third step: eigenvalues and eigenvectors
eigen.mat <- eigen(x = cov.mat)
eigen.mat
#