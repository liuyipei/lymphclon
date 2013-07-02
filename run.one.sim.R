#!/usr/bin/Rscript

if (all(c('infer.clonality.R', 'simulate.clonality.data.R') %in% list.files())) {
  write('using local source files', stderr())
  source('simulate.clonality.data.R')
  source('infer.clonality.R')
} else {
  library(lymphclon) 
}

curr.args <- commandArgs(trailingOnly = FALSE)
settings.index <- as.integer(curr.args[length(curr.args)]) # final argument
number.of.repeated.simulations <- as.integer(curr.args[length(curr.args) - 1]) # second last argument
if (!is.integer(number.of.repeated.simulations))
{
  number.of.repeated.simulations <- 1
}
if (!is.integer(settings.index))
{
  settings.index <- 3
}

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

base.num.cells.taken.vector <- round(c(rep(1e3, 6), rep(1e4, 2)) * num.cells.scaling)
num.cells.taken.vector <- rep(base.num.cells.taken.vector, sets.of.8)

#print(c(curr.settings, num.cells.taken.vector, clonal.power))
#clones <- 2e7 # 20m
#clones <- 5e5 # 500k
#clones <- 2e5 # 200k
clones <- 5e2 # 500 

for (i in c(1:number.of.repeated.simulations)) {

num.iterations <- 1 # refers to parameters for the fixed-point-clonality setting (fpc)

sim.data <- simulate.clonality.data(
  n = clones, 
  num.cells.taken.vector = num.cells.taken.vector,
  clonal.distribution.power = clonal.power)
#print ('data simulation complete')  

x <- sim.data$read.count.matrix
answer.opt.rer <- infer.clonality(read.count.matrix = x, 
  estimate.abundances = T, variance.method = 'usr.rer',
  internal.parameters = 
    list(use.squared.err.est = sim.data$replicate.squared.errs,
         compute.variances.d1jkn = T)
  )
internal.parameters <- answer.opt.rer$internal.parameters
internal.parameters[['use.replicate.var.est']] <- sim.data$replicate.squared.errs + sim.data$true.clonality

answer.opt.cov <- infer.clonality(read.count.matrix = x, 
  estimate.abundances = T, variance.method = 'usr.var',
  internal.parameters = internal.parameters)


answer.fpc.add <- infer.clonality(read.count.matrix = x, 
  estimate.abundances = T, variance.method = 'fpc.add', num.iterations = num.iterations,
  internal.parameters = internal.parameters)
answer.fpc.max <- infer.clonality(read.count.matrix = x, 
  estimate.abundances = T, variance.method = 'fpc.max', num.iterations = num.iterations,
  internal.parameters = internal.parameters)
answer.mle.cov <- infer.clonality(read.count.matrix = x, 
  estimate.abundances = T, variance.method = 'mle.cov',
  internal.parameters = internal.parameters)
answer.corpcor <- infer.clonality(read.count.matrix = x, 
  estimate.abundances = T, variance.method = 'corpcor',
  internal.parameters = internal.parameters)

simple.mean.abundances <- apply((x %*% diag(1/apply(x,2,sum)) ), 1,mean)

meta.cols <- c('power', 'replicates', 'clones', 'cells.scaling', 'chao2')
meta.values <- c(clonal.power,
  length(base.num.cells.taken.vector) * sets.of.8,
  clones, num.cells.scaling,
  internal.parameters$num.clones.est)
names(meta.values) <- meta.cols

baseline.cols <- c('true', 'bln', 'bln.abe2') # abe2: abundance error 2-norm
baseline.values <- c(
  sim.data$true.clonality,
  answer.fpc.add$simple.precision.clonality,
  sum((simple.mean.abundances - sim.data$true.clone.prob)^2))
names(baseline.values) <- baseline.cols

regularization.method.names <- 
  c('unregularized', 'ue.zr.full', 
  'eq.zr.half', 'ue.zr.half', 'eq.eq.half', 'ue.eq.half',
  'ue.mn.half', 'ue.mn.full', 'ue.mn.js1', 'mix1')
var.meth.names <- c('opt.rer', 'opt.cov', 'fpc.add', 'fpc.max', 'mle.cov', 'corpcor')

experiment.cols <- Reduce(c, 
  lapply(X = var.meth.names,
    FUN = function(curr.var.meth.name)
      {c(paste(curr.var.meth.name, regularization.method.names, sep = '.'),
         paste(curr.var.meth.name, 'abe2', sep = '.'))}
))


get.experiment.values.from.answer <- function(answer)
{
return(c(
  answer$regularized.estimates, 
  answer$mixture.clonality,
  sum((as.numeric(answer$estimated.abundances) - sim.data$true.clone.prob) ^ 2)
  ))
}

experiment.values <- c(
  get.experiment.values.from.answer(answer.opt.rer),
  get.experiment.values.from.answer(answer.opt.cov),
  get.experiment.values.from.answer(answer.fpc.add),
  get.experiment.values.from.answer(answer.fpc.max),
  get.experiment.values.from.answer(answer.mle.cov),
  get.experiment.values.from.answer(answer.corpcor)
  )
names(experiment.values) <- experiment.cols

all.cols <- c(meta.cols, baseline.cols, experiment.cols)
all.values <- c(meta.values, baseline.values, experiment.values)

#description <- paste(all.cols, collapse = '/')
description <- 'omitted'
values.str <- paste(sprintf("%.8e", all.values), collapse = ',')
print(paste(c(description, values.str), collapse = ','))

}
