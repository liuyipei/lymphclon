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
  num.iterations = 1,
  internal.parameters = list()
  ) {


if (length(internal.parameters$replicates > 0)) {
  replicates <- internal.parameters$replicates
} else {
  # normalize replicates, and get rep.grahm.matrix once
  replicates <- as.matrix(read.count.matrix)
  replicates <- t(t(replicates) / colSums(replicates))
}

if (length(internal.parameters$rep.grahm.matrix) > 0) {
  rep.grahm.matrix <- internal.parameters$rep.grahm.matrix
} else {
  rep.grahm.matrix <- t(replicates) %*% replicates
}

n <- nrow(replicates)
num.replicates <- ncol(replicates)
num.pairs <- num.replicates * (num.replicates - 1) / 2

if (length(internal.parameters$simple.precision.clonality) > 0) {
  simple.precision.clonality <- internal.parameters$simple.precision.clonality
} else {
  ## simple model where each read is independent
  reads.per.replicate <- as.numeric(apply(read.count.matrix, 2, sum))
  simple.precision.weights <- matrix(reads.per.replicate, nrow = num.replicates) %*% 
    matrix(reads.per.replicate, ncol = num.replicates)
  simple.precision.weights <- lower.tri(simple.precision.weights) * simple.precision.weights
  simple.precision.clonality <- sum(simple.precision.weights * rep.grahm.matrix) / sum(simple.precision.weights)
}

if (length(internal.parameters$use.squared.err.est) > 0) {
  use.squared.err.est <- as.matrix(internal.parameters$use.squared.err.est)
}

internal.parameters <- list(
  replicates = replicates,
  rep.grahm.matrix = rep.grahm.matrix,
  simple.precision.clonality = simple.precision.clonality,
  use.squared.err.est = use.squared.err.est
)

curr.clonality.score.estimate <- simple.precision.clonality
fpc.iter.estimates <- c()

for (curr.iter.number in 1:num.iterations) {

if (variance.method %in% c('fpc.1'))
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

} else if (variance.method %in% c('mle.1')) {
  replicates.cov <- cov(replicates)
}

# usr.1: argument use.squared.err.est specifies 1st order variances
# fpc.1: fixed point covariance, based on an unbiased clonality and self inner products
# mle.1: use maximum likelihood estimate
# corpcor.1: corpcor covvariance

Lambda.matrix <- matrix(data = NA, nrow = num.replicates, ncol = num.replicates)
ptinv.Lambda.matrix <- matrix(data = NA, nrow = num.replicates, ncol = num.replicates)
if (variance.method == 'usr.1') {
  epsilon.vec <- use.squared.err.est
  inv.eps.vec <- 1 / epsilon.vec
  diag(Lambda.matrix) <- inv.eps.vec
  diag(ptinv.Lambda.matrix) <- epsilon.vec
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
} else {
  stderr(sprintf('unknown variance method: %s\n', variance.method))
  stderr(sprintf('list of valid methods: usr.1, fpc.1, corpcor.1'))
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

pre.reg.matrix <- full.cov

# use the covariance matrix to compute the mvg MLE estimate, given the n-choose-2 estimators
print('...')
cov.to.clonality <- function(target.matrix, reg.coefs) {
  print(c('c2c start', reg.coefs))
  linear.combo.matrix <- (target.matrix      * reg.coefs)
                      +  (pre.reg.matrix * (1 - reg.coefs))
  root.prec.matrix <- sqrtm(ginv(linear.combo.matrix))
  numerator <- rep(1, num.pairs) %*% root.prec.matrix %*% estimator.vec.forcov
  denominator <- t(rep(1, num.pairs)) %*% root.prec.matrix %*% rep(1, num.pairs)
  print('c2c end')
  return(abs(numerator / denominator)) # the clonality score
}


#regularize the n-choose-2 by n-choose-2 covariance matrix
mean.eps2 <- mean(epsilon.vec)
min.eps2 <- min(epsilon.vec)

regularization.method.names <- 
  c('unregularized', 'ue.zr.full', 
  'eq.zr.half', 'ue.zr.half', 'eq.eq.half', 'ue.eq.half',
  'ue.mn.half', 'ue.mn.full', 'ue.mn.js1')
num.reg.methods <- length(regularization.method.names)
regularized.estimates <- rep(NA, num.reg.methods)
names(regularized.estimates) <- regularization.method.names

# initialize default values for regularization
reg.coefs <- rep(0.5, num.reg.methods)
names(reg.coefs) <- regularization.method.names

target.matrices <- lapply(regularization.method.names, function(x){pre.reg.matrix})
names(target.matrices) <- regularization.method.names

reg.coefs['unregularized'] <- 1
target.matrices[['unregularized']] <- pre.reg.matrix

reg.coefs['ue.zr.full'] <- 1
target.matrices[['ue.zr.full']] <- diag(diag(pre.reg.matrix))

reg.coefs['eq.zr.half'] <- 0.5
target.matrices[['eq.zr.half']] <- diag(rep(2 * mean.eps2, nrow(pre.reg.matrix)))

reg.coefs['ue.zr.half'] <- 0.5
target.matrices[['ue.zr.half']] <- diag(diag(pre.reg.matrix))

reg.coefs['eq.eq.half'] <- 0.5
target.matrices[['eq.eq.half']][pre.reg.matrix > 0] <- mean.eps2
diag(target.matrices[['eq.eq.half']]) <- 2 * mean.eps2

reg.coefs['ue.eq.half'] <- 0.5
target.matrices[['ue.eq.half']][pre.reg.matrix > 0] <- mean.eps2
diag(target.matrices[['ue.eq.half']]) <- 2 * mean.eps2

ue.mn.coef.denom <- sum((epsilon.vec - min.eps2) ^ 2) + (1e-14)
ue.mn.coef.numer <- sum((epsilon.vec - mean(epsilon.vec)) ^ 2) + (2e-14)
reg.coefs['ue.mn.half'] <- 0.5
reg.coefs['ue.mn.full'] <- 1
reg.coefs['ue.mn.js1'] <- ue.mn.coef.numer / ue.mn.coef.denom

ue.mn.matrix <- pre.reg.matrix
ue.mn.matrix[pre.reg.matrix > 0] <- min.eps2
diag(ue.mn.matrix) <- diag(pre.reg.matrix)
target.matrices[['ue.mn.half']] <- ue.mn.matrix
target.matrices[['ue.mn.full']] <- ue.mn.matrix
target.matrices[['ue.mn.js1']]  <- ue.mn.matrix

for (curr.reg in regularization.method.names)
{
  regularized.estimates[curr.reg] <-
    cov.to.clonality(target.matrices[[curr.reg]], reg.coefs[curr.reg])
}
#print(target.matrices)
#print(reg.coefs)

curr.clonality.score.estimate <- as.numeric(regularized.estimates['ue.zr.half']) 
# fpc: fixed point iterations. Usually 1 is best

fpc.iter.estimates <- append(fpc.iter.estimates, curr.clonality.score.estimate)
} # for (curr.iter.number in 1: num.iterations)

if (estimate.abundances) {
  return.results <- list(
    simple.precision.clonality = simple.precision.clonality, 
    lymphclon.clonality = regularized.estimates['ue.zr.half'],
    estimated.abundances = (replicates %*% inv.eps.vec) / sum(inv.eps.vec),
    estimated.squared.errs = epsilon.vec,
    estimated.precisions = inv.eps.vec,
    variance.method = variance.method,
    fpc.iter.estimates = fpc.iter.estimates,
    regularized.estimates = regularized.estimates,
    internal.parameters = internal.parameters)
} else {
  return.results <- list(
    simple.precision.clonality = simple.precision.clonality, 
    lymphclon.clonality = regularized.estimates['ue.zr.half'],
    estimated.abundances = 
      'The estimate.abundances parameter was set to false when infer.clonality was called.',
    estimated.squared.errs = epsilon.vec,
    estimated.precisions = inv.eps.vec,
    variance.method = variance.method,
    fpc.iter.estimates = fpc.iter.estimates,
    regularized.estimates = regularized.estimates,
    internal.parameters = 
      'The estimate.abundances parameter was set to false when infer.clonality was called. internal.parameters suppressed for brevity')
}

return (return.results)
}
