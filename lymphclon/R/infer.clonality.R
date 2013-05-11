#!/usr/bin/Rscript
library(expm)
library(MASS)
library(Matrix)
library(VGAM)

infer.clonality <- function(
  read.count.matrix, 
  regularization.clones = matrix(0, 0, ncol(read.count.matrix))) {

replicates <- read.count.matrix
n <- nrow(replicates)
num.replicates <- ncol(replicates)

for (i in 1:num.replicates)
{
    replicates[, i] <- read.count.matrix[, i] / sum(read.count.matrix[, i])
}

clonality.matrix <- t(replicates) %*% replicates
num.estimators <- num.replicates * (num.replicates-1) / 2

## simple model where each read is independent
reads.per.replicate <- as.numeric(apply(read.count.matrix, 2, sum))
simple.precision.weights <- matrix(reads.per.replicate, nrow=num.replicates) %*% matrix(reads.per.replicate, ncol=num.replicates)
simple.precision.weights <- lower.tri(simple.precision.weights) * simple.precision.weights
simple.precision.clonality <- sum(simple.precision.weights * clonality.matrix) / sum(simple.precision.weights)
#print(simple.precision.weights)
#print(clonality.matrix)
#print(simple.precision.weights * clonality.matrix)

estimators.rownums <- diag(c(1:num.replicates))  %*% matrix(1, num.replicates, num.replicates)
estimators.colnums <- t(estimators.rownums)
lowtri.indx <- as.vector(lower.tri(clonality.matrix))
estimators.vec <- as.vector(clonality.matrix)[lowtri.indx]
rownums.vec <- as.vector(estimators.rownums)[lowtri.indx]
colnums.vec <- as.vector(estimators.colnums)[lowtri.indx]

estimator.table.raw <- data.frame(theta = estimators.vec, rn = rownums.vec, cn = colnums.vec)

# the abbreviation etr2 refers to the self cross product of estimator.table.raw
etr2 <- merge(estimator.table.raw, estimator.table.raw, by = c()) 
# next step is to generate a unique entry for each unordered pairing of estimators that shared exactly one source replicate

indx <- etr2$cn.x == etr2$rn.y
estimator.table.pairs1 <- data.frame(cn1 = etr2$cn.x[indx], rn1 = etr2$rn.x[indx], cn2 = etr2$cn.y[indx], rn2 = etr2$rn.y[indx], shared = etr2$cn.x[indx], rn = etr2$rn.x[indx], cn = etr2$cn.y[indx])
indx <- (etr2$rn.x == etr2$rn.y) & (etr2$cn.x < etr2$cn.y)
estimator.table.pairs2 <- data.frame(cn1 = etr2$cn.x[indx], rn1 = etr2$rn.x[indx], cn2 = etr2$cn.y[indx], rn2 = etr2$rn.y[indx], shared = etr2$rn.x[indx], rn = etr2$cn.y[indx], cn = etr2$cn.x[indx])
indx <- (etr2$cn.x == etr2$cn.y) & (etr2$rn.x < etr2$rn.y)
estimator.table.pairs3 <- data.frame(cn1 = etr2$cn.x[indx], rn1 = etr2$rn.x[indx], cn2 = etr2$cn.y[indx], rn2 = etr2$rn.y[indx], shared = etr2$cn.x[indx], rn = etr2$rn.y[indx], cn = etr2$rn.x[indx])

estimator.table.pairs <- rbind(estimator.table.pairs1, estimator.table.pairs2, estimator.table.pairs3)
estimator.table.pairs$covariance <- diag(clonality.matrix)[estimator.table.pairs$shared]

row.col.to.indx <- function(r, c) {return (((r-1)*(r-2))/2 + c)} 
# this is for the lower half of cov mat, ie where r<c  r-1 choose 2 + c

# the reverse will be pre computed
indx.to.row <- rep(0, num.replicates*num.replicates)
indx.to.col <- rep(0, num.replicates*num.replicates)
curr.indx <- 0

# indx.to.row maps each integer indexing the n-choose-2 estimators to an integer row
for (r in 2:num.replicates) {
  for (c in 1:(r - 1)) {
    curr.indx <- row.col.to.indx(r, c) # could have been "++1" but just to be safe...
    indx.to.row[curr.indx] <- r
    indx.to.col[curr.indx] <- c 
  }
}

estimator.table.pairs$cov.r <- row.col.to.indx(estimator.table.pairs$rn1, estimator.table.pairs$cn1)
estimator.table.pairs$cov.c <- row.col.to.indx(estimator.table.pairs$rn2, estimator.table.pairs$cn2)
num.pairs <- num.replicates* (num.replicates-1) /2

cov.matrix.halfdone <- as.matrix(sparseMatrix(i = estimator.table.pairs$cov.r, j = estimator.table.pairs$cov.c, x = estimator.table.pairs$covariance))
# not necessarily upper or lower triangular -- just one of each (i, j) or (j, i) entry for i \ne j has been filled
cov.matrix <- as.matrix(cov.matrix.halfdone + t(cov.matrix.halfdone))
estimator.vec.forcov <- rep(0, num.pairs)
# fill out the variances
for (r in 2:num.replicates) {
  for (c in 1:(r-1)) {
    # could have been "++1" instead, but this code clarifies the intended use case for row.col.to.indx
    curr.indx <- row.col.to.indx(r, c) 
    cov.matrix[curr.indx, curr.indx] <- clonality.matrix[r, r] + clonality.matrix[c, c]
    estimator.vec.forcov[curr.indx] <- clonality.matrix[r, c]
  }
}

#use the covariance matrix to compute the mvg MLE estimate, given the n-choose-2 estimators
root.prec.matrix <- sqrtm(ginv(cov.matrix))
numerator <- rep(1, num.pairs) %*% root.prec.matrix %*% estimator.vec.forcov
denominator <- t(rep(1, num.pairs)) %*% root.prec.matrix %*% rep(1, num.pairs)
rao.blackwell.mvg.clonality <- abs(numerator / denominator)
return (c(simple.precision.clonality = simple.precision.clonality, rao.blackwell.mvg.clonality = rao.blackwell.mvg.clonality))
}
