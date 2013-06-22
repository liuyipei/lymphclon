#!/usr/bin/Rscript

if (all(c('infer.clonality.R', 'simulate.clonality.data.R') %in% list.files())) {
  write('using local source files', stderr())
  source('simulate.clonality.data.R')
  source('infer.clonality.R')
} else {
  library(lymphclon) 
}

curr.args <- commandArgs(trailingOnly = FALSE)
settings.index <- as.integer(curr.args[length(curr.args)])
#print(curr.args)

#sets.of.8.vec <- 1:4
sets.of.8.vec <- 1
clonal.power.vec <- - (0:20) / 20 
#num.cells.scaling.vec <- c(0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10)
num.cells.scaling.vec <- 1

raw.settings.table <- expand.grid(sets.of.8.vec, clonal.power.vec, num.cells.scaling.vec)
# at least 2 must match, the other parameter is floating
default.settings <- c(1, -0.7, 1.0)
valid.settings <- apply(raw.settings.table, 1, function(x){2 <= sum(1e-10 > abs(x - default.settings))}) 
settings.table <- raw.settings.table[valid.settings, ]

curr.settings <- as.numeric(settings.table[settings.index, ])
sets.of.8 <- curr.settings[1]
clonal.power <- curr.settings[2]
num.cells.scaling <- curr.settings[3]

base.num.cells.taken.vector <- round(c(rep(1000, 15), rep(10000, 3)) * num.cells.scaling)
num.cells.taken.vector <- rep(base.num.cells.taken.vector, sets.of.8)

#print(c(curr.settings, num.cells.taken.vector, clonal.power))
#clones <- 2e7 # 20m
clones <- 5e5 # 500k
#clones <- 5e2 # 500 

num.iterations <- 4
sim.data <- simulate.clonality.data(
  n = clones, 
  num.cells.taken.vector = num.cells.taken.vector,
  clonal.distribution.power = clonal.power)
#print ('data simulation complete')  

x <- sim.data$read.count.matrix
answer.opt1 <- infer.clonality(read.count.matrix = x, 
  estimate.abundances = T, variance.method = 'usr.1',
  use.squared.err.est = sim.data$replicate.squared.errs)
answer.fpc1 <- infer.clonality(read.count.matrix = x, 
  estimate.abundances = T, variance.method = 'fpc.1', num.iterations = num.iterations)
answer.mle1 <- infer.clonality(read.count.matrix = x, 
  estimate.abundances = T, variance.method = 'mle.1')
answer.cpc1 <- infer.clonality(read.count.matrix = x, 
  estimate.abundances = T, variance.method = 'corpcor.1')

simple.mean.abundances <- apply((x %*% diag(1/apply(x,2,sum)) ), 1,mean)
opt1.mean.abundances    <- as.numeric(answer.opt1$estimated.abundances)
fpc1.mean.abundances    <- as.numeric(answer.fpc1$estimated.abundances)
mle1.mean.abundances    <- as.numeric(answer.mle1$estimated.abundances)
cpc1.mean.abundances    <- as.numeric(answer.cpc1$estimated.abundances)

meta.cols <- c('power', 'replicates', 'clones', 'cells.scaling')
meta.values <- c(clonal.power,
  length(base.num.cells.taken.vector) * sets.of.8,
  clones, num.cells.scaling)
names(meta.values) <- meta.cols

baseline.cols <- c('true', 'bln', 'bln.abe2') # abe2: abundance error 2-norm
baseline.values <- c(
  sim.data$true.clonality,
  answer.fpc1$simple.precision.clonality,
  sum((simple.mean.abundances - sim.data$true.clone.prob)^2))
names(baseline.values) <- baseline.cols

regularization.method.names <- 
  c('unregularized', 'ue.zr.full', 
  'eq.zr.half', 'ue.zr.half', 'eq.eq.half', 'ue.eq.half',
  'ue.mn.half', 'ue.mn.full', 'ue.mn.js1')
var.meth.names <- c('opt1', 'fpc1', 'mle1', 'cpc1')

experiment.cols <- Reduce(c, 
  lapply(X = var.meth.names,
    FUN = function(curr.var.meth.name)
      {c(paste(curr.var.meth.name, regularization.method.names, 'err2', sep = '.'),
      paste(curr.var.meth.name, 'abe2', sep = '.'))}
))

experiment.values <- c(
  answer.opt1$regularized.estimates,
  sum((opt1.mean.abundances - sim.data$true.clone.prob)^2),

  answer.fpc1$regularized.estimates,
  sum((fpc1.mean.abundances - sim.data$true.clone.prob)^2),
  
  answer.mle1$regularized.estimates,
  sum((mle1.mean.abundances - sim.data$true.clone.prob)^2),

  answer.cpc1$regularized.estimates,
  sum((cpc1.mean.abundances - sim.data$true.clone.prob)^2)
  )
names(experiment.values) <- experiment.cols

all.cols <- c(meta.cols, baseline.cols, experiment.cols)
all.values <- c(meta.values, baseline.values, experiment.values)

description <- paste(all.cols, collapse = '/')
values.str <- paste(sprintf("%.8e", all.values), collapse = ',')
print(paste(c(description, values.str), collapse = ','))

