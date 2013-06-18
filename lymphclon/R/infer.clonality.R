#!/usr/bin/Rscript
library(expm)
library(MASS)
library(Matrix)
library(VGAM)
library(corpcor)

infer.clonality <- function(
  read.count.matrix, 
  variance.method = 'fpc.1',
  estimate.abundances = F,
  use.squared.err.est = c(),
  num.iterations = 1,
  regularization.method = ''
  ) {

# normalize replicates, and get rep.grahm.matrix once
replicates <- as.matrix(rbind(read.count.matrix))

n <- nrow(replicates)
num.replicates <- ncol(replicates)
replicates <- t(t(replicates) / colSums(replicates))
rep.grahm.matrix <- t(replicates) %*% replicates

num.pairs <- num.replicates * (num.replicates - 1) / 2

## simple model where each read is independent
reads.per.replicate <- as.numeric(apply(read.count.matrix, 2, sum))
simple.precision.weights <- matrix(reads.per.replicate, nrow = num.replicates) %*% 
  matrix(reads.per.replicate, ncol = num.replicates)
simple.precision.weights <- lower.tri(simple.precision.weights) * simple.precision.weights
simple.precision.clonality <- sum(simple.precision.weights * rep.grahm.matrix) / sum(simple.precision.weights)

curr.clonality.score.estimate <- simple.precision.clonality
rb.iter.estimates <- c()


if (length(use.squared.err.est) == num.replicates) {
  variance.method <- 'usr.1'
} else if (length(use.squared.err.est) == num.replicates * num.replicates) {
  variance.method <- 'usr.2'
}

for (curr.iter.number in 1:num.iterations) {

if (variance.method %in% c('fpc.1', 'fpc.2'))
{
  replicates.cov.off.diagonals <- 
    (curr.clonality.score.estimate 
    - diag(rep(curr.clonality.score.estimate, num.replicates))
    ) / n

  replicates.cov.diagonals <- 
    ifelse(
      diag(rep.grahm.matrix) > curr.clonality.score.estimate,
      diag(rep.grahm.matrix), 
      2 * curr.clonality.score.estimate - diag(rep.grahm.matrix)) / n

  replicates.cov <- diag(replicates.cov.diagonals) + replicates.cov.off.diagonals      

} else if (variance.method %in% c('mle.1', 'mle.2')) {
  replicates.cov <- cov(replicates)
}

if (length(use.squared.err.est > 0)) {
  use.squared.err.est <- as.matrix(use.squared.err.est)
}

# usr.1: argument use.squared.err.est specifies 1st order variances
# usr.2: argument use.squared.err.est specifies 2nd order variances
# fpc.1: fixed point iteration: 1st order variances
# fpc.2: fixed point iteration: 2nd order variances
# mle.1: use maximum likelihood estimate of 1st order variances
# mle.2: use maximum likelihood estimate of 2nd order variances
# corpcor.1: use shrinkage-based estimate of 1st order variances (see package corpcor)
# corpcor.2: use shrinkage-based estimate of 2nd order variances (see package corpcor)

Lambda.matrix <- matrix(data = NA, nrow = num.replicates, ncol = num.replicates)
ptinv.Lambda.matrix <- matrix(data = NA, nrow = num.replicates, ncol = num.replicates)
if (variance.method == 'usr.1') {
  epsilon.vec <- use.squared.err.est
  inv.eps.vec <- 1 / epsilon.vec
  diag(Lambda.matrix) <- inv.eps.vec
  diag(ptinv.Lambda.matrix) <- epsilon.vec
} else if (variance.method == 'usr.2') {
  # make this appear like Lambda
  epsilon.vec <- diag(use.squared.err.est)
  inv.eps.vec <- 1 / epsilon.vec
  Lambda.matrix <- -use.squared.err.est
  diag(Lambda.matrix) <- inv.eps.vec
  ptinv.Lambda.matrix <- 1 / Lambda.matrix
} else if (variance.method %in% c('fpc.1', 'mle.1', 'corpcor.1')) {
  if (variance.method %in% c('fpc.1', 'mle.1')) {
    inv.eps.vec <- diag(ginv(replicates.cov) / n)
  } else { 
    # variance.method == 'corpcor.1'
    capture.output(inv.eps.vec <- diag(invcov.shrink(replicates) / n))
  }

  epsilon.vec <- 1 / inv.eps.vec # the estimated conditional errors associated with each replicate
  diag(Lambda.matrix) <- inv.eps.vec
  diag(ptinv.Lambda.matrix) <- epsilon.vec
} else if (variance.method %in% c('fpc.2', 'mle.2', 'corpcor.2')) { 
  if (variance.method %in% c('fpc.2', 'mle.2')) {
    Lambda.matrix <- ginv(replicates.cov) / n
  } else {# variance.method == 'corpcor.2'
    capture.output(Lambda.matrix <- invcov.shrink(replicates) / n)
  }
 
  ptinv.Lambda.matrix <- 1 / Lambda.matrix
  inv.eps.vec <- diag(Lambda.matrix)
  epsilon.vec <- 1 / inv.eps.vec # the estimated conditional errors associated with each replicate
} else {
  print(sprintf('unknown variance method: %s\n', variance.method))
  print(sprintf('list of valid methods: usr.1, usr.2, fpc.1, fpc.2, corpcor.1, corpcor.2'))
}

contributions.to.replicate.cov.matrix <- -Lambda.matrix # negative conditional covariances
diag(contributions.to.replicate.cov.matrix) <- diag(ptinv.Lambda.matrix) # inverse conditional variances
#print('contributions.to.replicate.cov.matrix filled')

num.estimators <- num.replicates * (num.replicates - 1) / 2

estimators.rownums <- diag(c(1:num.replicates)) %*% matrix(1, num.replicates, num.replicates)
estimators.colnums <- t(estimators.rownums)
lowtri.indx <- as.vector(lower.tri(rep.grahm.matrix))
estimators.vec <- as.vector(rep.grahm.matrix)[lowtri.indx]
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
    estimator.vec.forcov[curr.indx] <- rep.grahm.matrix[r, c]
  }
}

# debug.data <- c()
#print('about to fill out full.cov')
full.cov <- matrix(data = 0, nrow = num.estimators, ncol = num.estimators)
# (n choose 2) by (n choose 2)
for (r1 in 2:num.replicates) {
  for (c1 in 1:(r1 - 1)) {
    for (r2 in 2:num.replicates) {
      for (c2 in 1:(r2 - 1)) {
        curr.indx1 <- row.col.to.indx(r1, c1) # could have been "++1" but just to be safe...
        curr.indx2 <- row.col.to.indx(r2, c2)
        curr.cov.components <- c(
          contributions.to.replicate.cov.matrix[r1, r2],
          contributions.to.replicate.cov.matrix[r1, c2],
          contributions.to.replicate.cov.matrix[c1, r2],
          contributions.to.replicate.cov.matrix[c1, c2])
        set.count <- sum(!is.na(curr.cov.components))
        curr.cov.components[is.na(curr.cov.components)] <- 0
        # debug.data <- rbind(debug.data, c(curr.indx1, curr.indx2, r1, c1, r2, c2, set.count, sum(curr.cov.components)))
        full.cov[curr.indx1, curr.indx2] <- sum(curr.cov.components)
      }
    }
  }
}
#print('full.cov filled')

mle.unconditioned.clonality = NA
if (F) {
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
    estimator.vec.forcov[curr.indx] <- rep.grahm.matrix[r, c]
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

#regularize the n-choose-2 by n-choose-2 covariance matrix
mean.eps2 <- mean(epsilon.vec)
min.eps2 <- min(epsilon.vec)
reg.coef <- 0.5
target.matrix <- cov.matrix
if (regularization.method == 'ue.zr.full') { # diagonal of unequals, project
  target.matrix <- diag(diag(cov.matrix))
  reg.coef <- 1
} else if (regularization.method == 'eq.zr.half') { # diagonal of equals
  target.matrix <- diag(rep(2 * mean.eps2, nrow(cov.matrix)))
} else if (regularization.method == 'ue.zr.half') { # diagonal of unequals
  target.matrix <- diag(diag(cov.matrix))
} else if (regularization.method == 'eq.eq.half') { # diagonals of equals, off diagonal of equals
  target.matrix <- cov.matrix
  target.matrix[cov.matrix > 0] <- mean.eps2
  diag(target.matrix) <- 2 * mean.eps2
} else if (regularization.method == 'ue.eq.half') { # diagnonals of unequals, off diagonals are equals
  target.matrix <- cov.matrix
  target.matrix[cov.matrix > 0] <- mean.eps2
  diag(target.matrix) <- diag(cov.matrix)
} else if (regularization.method == 'ue.mn.half') {
  target.matrix <- cov.matrix
  target.matrix[cov.matrix > 0] <- min.eps2
  diag(target.matrix) <- diag(cov.matrix)
} else if (regularization.method == 'ue.mn.full') { 
  target.matrix <- cov.matrix
  target.matrix[cov.matrix > 0] <- min.eps2
  diag(target.matrix) <- diag(cov.matrix)
  reg.coef <- 1
} else if (regularization.method == 'ue.mn.js1') {
  target.matrix <- cov.matrix
  target.matrix[cov.matrix > 0] <- min.eps2
  diag(target.matrix) <- diag(cov.matrix)
  reg.coef.denom <- sum((epsilon.vec - min.eps2) ^ 2)
  reg.coef.numer <- sum((epsilon.vec - mean(epsilon.vec)) ^ 2)
  reg.coef <- reg.coef.numer / reg.coef.denom
} else if (nchar(regularization.method) > 0) {
  write(sprintf('unrecognized regularization method: %s. Using none.', regularization.method), stderr())
}

cov.matrix <- (target.matrix * reg.coef)
           +  (cov.matrix    * (1 - reg.coef))


# use the covariance matrix to compute the mvg MLE estimate, given the n-choose-2 estimators
root.prec.matrix <- sqrtm(ginv(cov.matrix))
numerator <- rep(1, num.pairs) %*% root.prec.matrix %*% estimator.vec.forcov
denominator <- t(rep(1, num.pairs)) %*% root.prec.matrix %*% rep(1, num.pairs)
rao.blackwell.mvg.clonality <- abs(numerator / denominator)

curr.clonality.score.estimate <- as.numeric(rao.blackwell.mvg.clonality) # update

rb.iter.estimates <- append(rb.iter.estimates, curr.clonality.score.estimate)
} #for (curr.iter.number in 1: num.iterations)

if (estimate.abundances) {
  return.results <- list(
    simple.precision.clonality = simple.precision.clonality, 
    mle.unconditioned.clonality = mle.unconditioned.clonality,
    rao.blackwell.mvg.clonality = rao.blackwell.mvg.clonality,
    estimated.abundances = (replicates %*% inv.eps.vec) / sum(inv.eps.vec),
    estimated.squared.errs = epsilon.vec,
    estimated.precisions = inv.eps.vec,
    variance.method = variance.method,
    rb.iter.estimates = rb.iter.estimates,
    reg.coef = reg.coef,
    regularization.method = regularization.method)
} else {
  return.results <- list(
    simple.precision.clonality = simple.precision.clonality, 
    mle.unconditioned.clonality = mle.unconditioned.clonality,
    rao.blackwell.mvg.clonality = rao.blackwell.mvg.clonality,
    contributions.to.replicate.cov.matrix,
    variance.method = variance.method,
    rb.iter.estimates = rb.iter.estimates,
    reg.coef = reg.coef,
    regularization.method = regularization.method)
}

return (return.results)
}
