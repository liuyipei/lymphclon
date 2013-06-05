#!/usr/bin/Rscript
library(expm)
library(MASS)
library(Matrix)
library(VGAM)
library(corpcor)

infer.clonality <- function(
  read.count.matrix, 
  regularization.clones = matrix(0, 0, ncol(read.count.matrix)),
  variance.method = 'loo.2',
  estimate.abundances = F,
  loo.squared.err.est = c(),
  use.squared.err.est = c()) {

replicates <- rbind(read.count.matrix, regularization.clones)
n <- nrow(replicates)
num.replicates <- ncol(replicates)
num.pairs <- num.replicates * (num.replicates - 1) / 2

for (i in 1:num.replicates)
{
    replicates[, i] <- read.count.matrix[, i] / sum(read.count.matrix[, i])
}

clonality.matrix <- t(replicates) %*% replicates

if (length(use.squared.err.est > 0)) {
  use.squared.err.est <- as.matrix(use.squared.err.est)
}

# usr.1: argument use.squared.err.est specifies 1st order variances
# usr.2: argument use.squared.err.est specifies 2nd order variances
# mle.1: use maximum likelihood estimate of 1st order variances
# mle.2: use maximum likelihood estimate of 2nd order variances
# loo.1: use leave-one-out estimate of 1st order variances
# loo.2: use leave-one-out estimate of 2nd order variances
# corpcor.1: use shrinkage-based estimate of 1st order variances (see package corpcor)
# corpcor.2: use shrinkage-based estimate of 2nd order variances (see package corpcor)

if (length(use.squared.err.est) == num.replicates) {
  variance.method <- 'usr.1'
} else if (length(use.squared.err.est) == num.replicates * num.replicates) {
  variance.method <- 'usr.2'
} else if (length(loo.squared.err.est) == 1) { 
  variance.method <- ifelse(loo.squared.err.est, 'loo.2', 'mle.1')
}

Lambda.matrix <- matrix(data = NA, nrow = num.replicates, ncol = num.replicates)
ptinv.Lambda.matrix <- matrix(data = NA, nrow = num.replicates, ncol = num.replicates)
if (variance.method == 'usr.1') {
  epsilon.vec <- use.squared.err.est
  inv.eps.vec <- 1 / epsilon.vec
  diag(Lambda.matrix) <- inv.eps.vec
  ptinv.Lambda.matrix <- 1 / Lambda.matrix
} else if (variance.method == 'usr.2') {
  # make this appear like Lambda
  epsilon.vec <- diag(use.squared.err.est)
  inv.eps.vec <- 1 / epsilon.vec
  Lambda.matrix <- -use.squared.err.est
  diag(Lambda.matrix) <- inv.eps.vec
  ptinv.Lambda.matrix <- 1 / Lambda.matrix
} else if (variance.method %in% c('mle.1', 'corpcor.1')) {
  if (variance.method == 'mle.1') {
    inv.eps.vec <- diag(ginv(clonality.matrix))
  } else { 
    # variance.method == 'corpcor.1'
    capture.output(inv.eps.vec <- diag(invcov.shrink(replicates) / n))
  }

  epsilon.vec <- 1 / inv.eps.vec # the estimated conditional errors associated with each replicate
  diag(Lambda.matrix) <- inv.eps.vec
  diag(ptinv.Lambda.matrix) <- epsilon.vec
} else if (variance.method %in% c('mle.2', 'corpcor.2')) { 
  if (variance.method == 'mle.2') {
    Lambda.matrix <- ginv(clonality.matrix)
  } else {# variance.method == 'corpcor.2'
    capture.output(Lambda.matrix <- invcov.shrink(replicates) / n)
  }
 
  ptinv.Lambda.matrix <- 1 / Lambda.matrix
  inv.eps.vec <- diag(Lambda.matrix)
  epsilon.vec <- 1 / inv.eps.vec # the estimated conditional errors associated with each replicate
} else if (variance.method %in% c('loo.1', 'loo.2')) { 
  loo.inv.eps.vec <- rep(0, num.replicates)
  loo.eps.vec <- rep(0, num.replicates)

  loo.Lambda.matrix <- matrix(data = 0, nrow = num.replicates, ncol = num.replicates)
  loo.ptinv.Lambda.matrix <- matrix(data = 0, nrow = num.replicates, ncol = num.replicates)  
  for (i in (1:num.replicates)) {
    curr.sub.Lambda.matrix <- ginv(clonality.matrix[-i, -i])
    loo.Lambda.matrix[-i, -i] <- loo.Lambda.matrix[-i, -i] + curr.sub.Lambda.matrix
    loo.ptinv.Lambda.matrix[-i, -i] <- loo.ptinv.Lambda.matrix[-i, -i] + (1 / curr.sub.Lambda.matrix)

    curr.inv.eps <- diag(curr.sub.Lambda.matrix)
    loo.inv.eps.vec[-i] <- loo.inv.eps.vec[-i] + curr.inv.eps
    loo.eps.vec[-i] <- loo.eps.vec[-i] + (1 / curr.inv.eps)
  }

  loo.Lambda.matrix <- loo.Lambda.matrix / (num.replicates - 1)
  loo.ptinv.Lambda.matrix <- loo.ptinv.Lambda.matrix / (num.replicates - 1)

  loo.inv.eps.vec <- loo.inv.eps.vec / (num.replicates - 1)
  loo.eps.vec <- loo.eps.vec / (num.replicates - 1)

  epsilon.vec <- loo.eps.vec
  inv.eps.vec <- loo.inv.eps.vec

  # "default" value if variance.method == loo.1  
  diag(Lambda.matrix) <- diag(loo.Lambda.matrix)
  diag(ptinv.Lambda.matrix) <- diag(loo.ptinv.Lambda.matrix)

  if (variance.method == 'loo.2') {
    Lambda.matrix <- loo.Lambda.matrix
    ptinv.Lambda.matrix <- loo.ptinv.Lambda.matrix
  }
} else {
  print(sprintf('unknown variance method: %s\n', variance.method))
  print(sprintf('list of valid methods: usr.1, usr.2, mle.1, mle.2, loo.1, loo.2, corpcor.1, corpcor.2'))
}
conditional.replicate.cov.matrix <- -Lambda.matrix # negative conditional covariances
diag(conditional.replicate.cov.matrix) <- diag(ptinv.Lambda.matrix) # inverse conditional variances

num.estimators <- num.replicates * (num.replicates - 1) / 2

## simple model where each read is independent
reads.per.replicate <- as.numeric(apply(read.count.matrix, 2, sum))
simple.precision.weights <- matrix(reads.per.replicate, nrow = num.replicates) %*% 
  matrix(reads.per.replicate, ncol = num.replicates)
simple.precision.weights <- lower.tri(simple.precision.weights) * simple.precision.weights
simple.precision.clonality <- sum(simple.precision.weights * clonality.matrix) / sum(simple.precision.weights)
#print(simple.precision.weights)
#print(clonality.matrix)
#print(simple.precision.weights * clonality.matrix)

estimators.rownums <- diag(c(1:num.replicates)) %*% matrix(1, num.replicates, num.replicates)
estimators.colnums <- t(estimators.rownums)
lowtri.indx <- as.vector(lower.tri(clonality.matrix))
estimators.vec <- as.vector(clonality.matrix)[lowtri.indx]
rownums.vec <- as.vector(estimators.rownums)[lowtri.indx]
colnums.vec <- as.vector(estimators.colnums)[lowtri.indx]

row.col.to.indx <- function(r, c) {return (((r - 1) * (r - 2)) / 2 + c)} 
# this is for the lower half of cov mat, ie where r<c  r-1 choose 2 + c

# the reverse will be pre computed
indx.to.row <- rep(0, num.replicates * num.replicates)
indx.to.col <- rep(0, num.replicates * num.replicates)
curr.indx <- 0

# indx.to.row maps each integer indexing the n-choose-2 estimators to an integer row
for (r in 2:num.replicates) {
  for (c in 1:(r - 1)) {
    curr.indx <- row.col.to.indx(r, c) # could have been "++1" but just to be safe...
    indx.to.row[curr.indx] <- r
    indx.to.col[curr.indx] <- c 
  }
}

############### NEW LOGIC
# fill out the actual estimators
estimator.vec.forcov <- rep(0, num.pairs)
for (r in 2:num.replicates) {
  for (c in 1:(r-1)) {
    curr.indx <- row.col.to.indx(r, c) 
    estimator.vec.forcov[curr.indx] <- clonality.matrix[r, c]
  }
}

debug.data <- c()
full.cov <- matrix(data = 0, nrow = num.estimators, ncol = num.estimators)
# (n choose 2) by (n choose 2)
for (r1 in 2:num.replicates) {
  for (c1 in 1:(r1 - 1)) {
    for (r2 in 2:num.replicates) {
      for (c2 in 1:(r2 - 1)) {
        curr.indx1 <- row.col.to.indx(r1, c1) # could have been "++1" but just to be safe...
        curr.indx2 <- row.col.to.indx(r2, c2)
        curr.cov.components <- c(
          conditional.replicate.cov.matrix[r1, r2],
          conditional.replicate.cov.matrix[r1, c2],
          conditional.replicate.cov.matrix[c1, r2],
          conditional.replicate.cov.matrix[c1, c2])
        set.count <- sum(!is.na(curr.cov.components))
        curr.cov.components[is.na(curr.cov.components)] <- 0
        debug.data <- rbind(debug.data, c(curr.indx1, curr.indx2, r1, c1, r2, c2, 
          set.count, sum(curr.cov.components)))
        full.cov[curr.indx1, curr.indx2] <- sum(curr.cov.components)
      }
    }
  }
}

if (T) {
#### OLD LOGIC
estimator.table.raw <- data.frame(theta = estimators.vec, rn = rownums.vec, cn = colnums.vec)

# the abbreviation etr2 refers to the self cross product of estimator.table.raw
etr2 <- merge(estimator.table.raw, estimator.table.raw, by = c()) 
# next step is to generate a unique entry for each unordered pairing of estimators that shared exactly one source replicate

indx <- etr2$cn.x == etr2$rn.y
estimator.table.pairs1 <- data.frame(cn1 = etr2$cn.x[indx], rn1 = etr2$rn.x[indx], cn2 = etr2$cn.y[indx], 
  rn2 = etr2$rn.y[indx], shared = etr2$cn.x[indx], rn = etr2$rn.x[indx], cn = etr2$cn.y[indx])
indx <- (etr2$rn.x == etr2$rn.y) & (etr2$cn.x < etr2$cn.y)
estimator.table.pairs2 <- data.frame(cn1 = etr2$cn.x[indx], rn1 = etr2$rn.x[indx], cn2 = etr2$cn.y[indx], 
  rn2 = etr2$rn.y[indx], shared = etr2$rn.x[indx], rn = etr2$cn.y[indx], cn = etr2$cn.x[indx])
indx <- (etr2$cn.x == etr2$cn.y) & (etr2$rn.x < etr2$rn.y)
estimator.table.pairs3 <- data.frame(cn1 = etr2$cn.x[indx], rn1 = etr2$rn.x[indx], cn2 = etr2$cn.y[indx], 
  rn2 = etr2$rn.y[indx], shared = etr2$cn.x[indx], rn = etr2$rn.y[indx], cn = etr2$rn.x[indx])

estimator.table.pairs <- rbind(estimator.table.pairs1, estimator.table.pairs2, estimator.table.pairs3)
estimator.table.pairs$covariance <- epsilon.vec[estimator.table.pairs$shared]

estimator.table.pairs$cov.r <- row.col.to.indx(estimator.table.pairs$rn1, estimator.table.pairs$cn1)
estimator.table.pairs$cov.c <- row.col.to.indx(estimator.table.pairs$rn2, estimator.table.pairs$cn2)

cov.matrix.halfdone <- as.matrix(sparseMatrix(
  i = estimator.table.pairs$cov.r, 
  j = estimator.table.pairs$cov.c, 
  x = estimator.table.pairs$covariance))
# not necessarily upper or lower triangular -- just one of each (i, j) or (j, i) entry for i \ne j has been filled
cov.matrix <- as.matrix(cov.matrix.halfdone + t(cov.matrix.halfdone))
estimator.vec.forcov <- rep(0, num.pairs)
# fill out the variances
for (r in 2:num.replicates) {
  for (c in 1:(r-1)) {
    # could have been "++1" instead, but this code clarifies the intended use case for row.col.to.indx
    curr.indx <- row.col.to.indx(r, c) 
    cov.matrix[curr.indx, curr.indx] <- epsilon.vec[r] + epsilon.vec[c]
    estimator.vec.forcov[curr.indx] <- clonality.matrix[r, c]
  }
}
root.prec.matrix <- sqrtm(ginv(cov.matrix))
numerator <- rep(1, num.pairs) %*% root.prec.matrix %*% estimator.vec.forcov
denominator <- t(rep(1, num.pairs)) %*% root.prec.matrix %*% rep(1, num.pairs)
mle.unconditioned.clonality <- abs(numerator / denominator)
old.cov.matrix<-cov.matrix
##### end old logic
} # if(F)

cov.matrix <- full.cov

# use the covariance matrix to compute the mvg MLE estimate, given the n-choose-2 estimators
root.prec.matrix <- sqrtm(ginv(cov.matrix))
numerator <- rep(1, num.pairs) %*% root.prec.matrix %*% estimator.vec.forcov
denominator <- t(rep(1, num.pairs)) %*% root.prec.matrix %*% rep(1, num.pairs)
rao.blackwell.mvg.clonality <- abs(numerator / denominator)

if (estimate.abundances) {
  return.results <- list(
    simple.precision.clonality = simple.precision.clonality, 
    mle.unconditioned.clonality = mle.unconditioned.clonality,
    rao.blackwell.mvg.clonality = rao.blackwell.mvg.clonality,
    estimated.abundances = (replicates %*% inv.eps.vec) / sum(inv.eps.vec),
    estimated.squared.errs = epsilon.vec,
    estimated.precisions = inv.eps.vec,
    variance.method = variance.method)
} else {
  return.results <- list(
    simple.precision.clonality = simple.precision.clonality, 
    mle.unconditioned.clonality = mle.unconditioned.clonality,
    rao.blackwell.mvg.clonality = rao.blackwell.mvg.clonality,
    conditional.replicate.cov.matrix,
    variance.method = variance.method)
}

return (return.results)
}
