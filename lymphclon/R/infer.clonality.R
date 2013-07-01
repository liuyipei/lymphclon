#!/usr/bin/Rscript
library(expm)
library(MASS)
library(Matrix)
library(VGAM)
library(corpcor)

infer.clonality <- function(
  read.count.matrix, 
  variance.method = 'fpc.add',
  estimate.abundances = F,
  num.iterations = 1,
  internal.parameters = list()
  ) {

take.pos <- function(x)
{
  x[x < 0] <- 0
  return(x)
}

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
  simple.precision.clonality <- sum(simple.precision.weights * rep.grahm.matrix) / 
    sum(simple.precision.weights)
}

if (length(internal.parameters$num.clones.est) > 0) {
  num.clones.est <- internal.parameters$num.clones.est
} else {
  int.chao <- function (x) { # taken from fossil package
    so <- length(x[x > 0])
    s1 <- length(x[x == 1])
    s2 <- length(x[x == 2])
    if ((s1 - s2)^2 == (s1 + s2)^2) 
        return(so + s1 * (s1 - 1)/((s2 + 1) * 2))
    else return(so + s1^2/(s2 * 2))
  }

  num.clones.est <- int.chao(apply(replicates > 0, 1, sum))
}
#print(num.clones.est)

if (length(internal.parameters$use.squared.err.est) > 0) {
  use.squared.err.est <- as.matrix(internal.parameters$use.squared.err.est)
} else {
  use.squared.err.est <- c()
}

if (length(internal.parameters$use.replicate.var.est) > 0) {
  use.replicate.var.est <- as.matrix(internal.parameters$use.replicate.var.est)
} else {
  use.replicate.var.est <- c()
}

internal.parameters <- list(
  replicates = replicates,
  rep.grahm.matrix = rep.grahm.matrix,
  simple.precision.clonality = simple.precision.clonality,
  num.clones.est = num.clones.est, 
  use.squared.err.est = use.squared.err.est,
  use.replicate.var.est = use.replicate.var.est
)

curr.clonality.score.estimate <- simple.precision.clonality
fpc.iter.estimates <- c()

for (curr.iter.number in 1:num.iterations) {

replicates.cov.off.diagonal.value <-
  take.pos(curr.clonality.score.estimate - (1 / num.clones.est)) 

replicates.cov.off.diagonals <- 
    matrix(replicates.cov.off.diagonal.value, num.replicates, num.replicates)
    - diag(rep(curr.clonality.score.estimate, num.replicates))

if (variance.method %in% c('fpc.add'))
{ # diagonal is the clonality score plus the abs difference of the self-inner products from it
  replicates.cov.diagonals <- 
    ifelse(
      diag(rep.grahm.matrix) > replicates.cov.off.diagonal.value,
      diag(rep.grahm.matrix), 
      2 * replicates.cov.off.diagonal.value - diag(rep.grahm.matrix))
    + rep((1 / num.clones.est), num.replicates) 
    # regularize by smallest possible clonality, given number of clones seen  
    # print(replicates.cov.off.diagonals)
    replicates.cov <- diag(replicates.cov.diagonals) + replicates.cov.off.diagonals      
} else if (variance.method %in% c('fpc.max')) 
{ # diagonal is the max of clonality score, or the self-inner products
  replicates.cov.diagonals <- 
    ifelse(
      diag(rep.grahm.matrix) > replicates.cov.off.diagonal.value,
      diag(rep.grahm.matrix), 
      replicates.cov.off.diagonal.value)
    + rep((1 / num.clones.est), num.replicates)  
  # print(replicates.cov.off.diagonals)
  # regularize by smallest possible clonality, given number of clones seen
  replicates.cov <- diag(replicates.cov.diagonals) + replicates.cov.off.diagonals      
} else if (variance.method %in% c('mle.cov')) {
  replicates.cov <- cov(replicates) * n
} else if (variance.method %in% c('usr.var')) {
  replicates.cov.diagonals <- use.replicate.var.est
  replicates.cov <- diag(replicates.cov.diagonals) + replicates.cov.off.diagonals      
}

# usr.rer: internal.parameters$use.squared.err.est specifies conditional variances of replicates
# usr.var internal.parameters$use.replicate.var.est specifies variances of replicates
# fpc.add: fixed point covariance, based on an unbiased clonality and self inner products
# mle.cov: use maximum likelihood estimate
# corpcor: corpcor covvariance

Lambda.matrix <- matrix(data = NA, nrow = num.replicates, ncol = num.replicates)
ptinv.Lambda.matrix <- matrix(data = NA, nrow = num.replicates, ncol = num.replicates)
if (variance.method == 'usr.rer') {
  epsilon.vec <- use.squared.err.est
  inv.eps.vec <- 1 / epsilon.vec
  diag(Lambda.matrix) <- inv.eps.vec
  diag(ptinv.Lambda.matrix) <- epsilon.vec
} else if (variance.method %in% c('fpc.add', 'fpc.max', 'mle.cov', 'corpcor', 'usr.var')) {
  if (variance.method %in% c('fpc.add', 'mle.cov', 'fpc.max', 'usr.var')) {
    inv.eps.vec <- diag(ginv(replicates.cov))
  } else { 
    # variance.method == 'corpcor'
    capture.output(inv.eps.vec <- diag(invcov.shrink(replicates)))
  }

  epsilon.vec <- 1 / inv.eps.vec # the estimated conditional errors associated with each replicate
  diag(Lambda.matrix) <- inv.eps.vec
  diag(ptinv.Lambda.matrix) <- epsilon.vec
} else {
  write(sprintf('unknown variance method: %s\n', variance.method), stderr())
  write('list of valid methods: usr.rer, fpc.add, fpc.max, corpcor', stderr())
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

unreg.matrix <- full.cov

# use the covariance matrix to compute the mvg MLE estimate, given the n-choose-2 estimators
cov.to.clonality <- function(target.matrix, reg.coefs) {
  target.term <- (target.matrix * reg.coefs)
  unreg.term <- (unreg.matrix * (1 - reg.coefs))
  linear.combo.matrix <- target.term + unreg.term
  root.prec.matrix <- sqrtm(ginv(linear.combo.matrix))
  numerator <- rep(1, num.pairs) %*% root.prec.matrix %*% estimator.vec.forcov
  denominator <- t(rep(1, num.pairs)) %*% root.prec.matrix %*% rep(1, num.pairs)
  return(abs(numerator / denominator)) # the clonality score
} # cov.to.clonality

# regularize the n-choose-2 by n-choose-2 covariance matrix
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

# populate the list with default values to be overwritten
target.matrices <- lapply(regularization.method.names, function(x){unreg.matrix}) 
names(target.matrices) <- regularization.method.names

reg.coefs['unregularized'] <- 0 # doesn't matter
target.matrices[['unregularized']] <- unreg.matrix

reg.coefs['ue.zr.full'] <- 1
target.matrices[['ue.zr.full']] <- diag(diag(unreg.matrix))

reg.coefs['eq.zr.half'] <- 0.5
target.matrices[['eq.zr.half']] <- diag(rep(2 * mean.eps2, nrow(unreg.matrix)))

reg.coefs['ue.zr.half'] <- 0.5
target.matrices[['ue.zr.half']] <- diag(diag(unreg.matrix))

reg.coefs['eq.eq.half'] <- 0.5
target.matrices[['eq.eq.half']][unreg.matrix > 0] <- mean.eps2
diag(target.matrices[['eq.eq.half']]) <- 2 * mean.eps2

reg.coefs['ue.eq.half'] <- 0.5
target.matrices[['ue.eq.half']][unreg.matrix > 0] <- mean.eps2
diag(target.matrices[['ue.eq.half']]) <- 2 * mean.eps2

ue.mn.coef.denom <- sum((epsilon.vec - min.eps2) ^ 2) + (1e-14)
ue.mn.coef.numer <- sum((epsilon.vec - mean(epsilon.vec)) ^ 2) + (2e-14)
reg.coefs['ue.mn.half'] <- 0.5
reg.coefs['ue.mn.full'] <- 1
reg.coefs['ue.mn.js1'] <- ue.mn.coef.numer / ue.mn.coef.denom

ue.mn.matrix <- unreg.matrix
ue.mn.matrix[unreg.matrix > 0] <- min.eps2
diag(ue.mn.matrix) <- diag(unreg.matrix)
target.matrices[['ue.mn.half']] <- ue.mn.matrix
target.matrices[['ue.mn.full']] <- ue.mn.matrix
target.matrices[['ue.mn.js1']]  <- ue.mn.matrix

for (curr.reg in regularization.method.names)
{
  regularized.estimates[curr.reg] <-
    cov.to.clonality(target.matrices[[curr.reg]], reg.coefs[curr.reg])
}

curr.clonality.score.estimate <- as.numeric(regularized.estimates['ue.zr.half']) 
# fpc: fixed point iterations. Usually 1 is best

fpc.iter.estimates <- append(fpc.iter.estimates, curr.clonality.score.estimate)
} # for (curr.iter.number in 1: num.iterations)

# use jackknife to estimate variance
compute.variances.d1jkn <- F
if (length(internal.parameters$compute.variances.d1jkn > 0)) {
  compute.variances.d1jkn <- internal.parameters$compute.variances.d1jkn
}

simple.clonality.variance <- NA
regularized.clonality.variance <- NA
if (compute.variances.d1jkn) {
subset.simple2.sum  <- 0
subset.reg2.sum     <- 0
for (i in c(1:num.replicates)) {
  curr.subset.internal.parameters <- list(
    replicates = replicates[, -i],
    rep.grahm.matrix = rep.grahm.matrix[-i, -i])
  curr.subset.clonality <- infer.clonality(
      read.count.matrix = read.count.matrix[, -i], 
      variance.method = variance.method,
      internal.parameters = curr.subset.internal.parameters)
  subset.simple2.sum <- subset.simple2.sum +
    (curr.subset.clonality$simple.precision.clonality - simple.precision.clonality) ^ 2
  subset.reg2.sum <- subset.reg2.sum +
    (curr.subset.clonality$regularized.estimates - regularized.estimates) ^ 2
} # for (i in c(1:num.replicates))
n.var.coef <- (num.replicates - 1) / num.replicates
simple.clonality.variance <- n.var.coef * subset.simple2.sum
regularized.clonality.variance <- n.var.coef * subset.reg2.sum
} # if (compute.variances.d1jkn)

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
    internal.parameters = internal.parameters,
    simple.clonality.variance = simple.clonality.variance,
    regularized.clonality.variance = regularized.clonality.variance
    )
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
    internal.parameters = internal.parameters,
    simple.clonality.variance = simple.clonality.variance,
    regularized.clonality.variance = regularized.clonality.variance
    )
}

return (return.results)
}
