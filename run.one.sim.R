#!/usr/bin/Rscript
#library(lymphclon)

source('simulate.clonality.data.R')
source('infer.clonality.R')

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

num.iterations <- 10
sim.data <- simulate.clonality.data(
  n = clones, 
  num.cells.taken.vector = num.cells.taken.vector,
  clonal.distribution.power = clonal.power)
#print ('data simulation complete')  

x <- sim.data$read.count.matrix
answer.opt1 <- infer.clonality(read.count.matrix = x, 
  estimate.abundances = T, variance.method = 'usr.1',
  use.squared.err.est = sim.data$replicate.squared.errs)
answer.opt2 <- infer.clonality(read.count.matrix = x, 
  estimate.abundances = T, variance.method = 'usr.2',
  use.squared.err.est = sim.data$grahm.errs)

answer.mle1 <- infer.clonality(read.count.matrix = x, 
  estimate.abundances = T, variance.method = 'mle.1', num.iterations = num.iterations)
answer.mle2 <- infer.clonality(read.count.matrix = x, 
  estimate.abundances = T, variance.method = 'mle.2', num.iterations = num.iterations)
answer.loo1 <- infer.clonality(read.count.matrix = x, 
  estimate.abundances = T, variance.method = 'loo.1', num.iterations = num.iterations)
answer.loo2 <- infer.clonality(read.count.matrix = x, 
  estimate.abundances = T, variance.method = 'loo.2', num.iterations = num.iterations)

answer.cpc1 <- infer.clonality(read.count.matrix = x, 
  estimate.abundances = T, variance.method = 'corpcor.1')
answer.cpc2 <- infer.clonality(read.count.matrix = x, 
  estimate.abundances = T, variance.method = 'corpcor.2')

#row <- as.matrix(t(
# c(answer.loo['simple.precision.clonality'], 
# answer.loo['rb.iter.estimates'],
# true.clonality <- sim.data$true.clonality)
#))

simple.mean.abundances <- apply((x %*% diag(1/apply(x,2,sum)) ), 1,mean)
opt1.mean.abundances    <- as.numeric(answer.opt1$estimated.abundances)
opt2.mean.abundances    <- as.numeric(answer.opt2$estimated.abundances)
mle1.mean.abundances    <- as.numeric(answer.mle1$estimated.abundances)
mle2.mean.abundances    <- as.numeric(answer.mle2$estimated.abundances)
loo1.mean.abundances    <- as.numeric(answer.loo1$estimated.abundances)
loo2.mean.abundances    <- as.numeric(answer.loo2$estimated.abundances)
cpc1.mean.abundances    <- as.numeric(answer.cpc1$estimated.abundances)
cpc2.mean.abundances    <- as.numeric(answer.cpc2$estimated.abundances)

meta.cols <- c('power', 'replicates', 'clones', 'cells.scaling')
meta.values <- c(clonal.power,
  length(base.num.cells.taken.vector) * sets.of.8,
  clones, num.cells.scaling)
names(meta.values) <- meta.cols

baseline.cols <- c('true', 'bln', 'bln.abe2') # abe2: abundance error 2-norm
baseline.values <- c(
  sim.data$true.clonality,
  answer.mle1$simple.precision.clonality,
  sum((simple.mean.abundances - sim.data$true.clone.prob)^2))
names(baseline.values) <- baseline.cols

experiment.cols <- 
  c('opt1', 'opt1.abe2', 'opt2', 'opt2.abe2', 
    paste('mle1', c(1:num.iterations), sep = '.'), 'mle1.abe2', 
    paste('mle2', c(1:num.iterations), sep = '.'), 'mle2.abe2', 
    paste('loo1', c(1:num.iterations), sep = '.'), 'loo1.abe2', 
    paste('loo2', c(1:num.iterations), sep = '.'), 'loo2.abe2', 
    'cpc1', 'cpc1.abe2', 'cpc2', 'cpc2.abe2')
experiment.values <- c(
  answer.opt1$rb.iter.estimates,
  sum((opt1.mean.abundances - sim.data$true.clone.prob)^2),
  answer.opt2$rb.iter.estimates,
  sum((opt2.mean.abundances - sim.data$true.clone.prob)^2),

  answer.mle1$rb.iter.estimates,
  sum((mle1.mean.abundances - sim.data$true.clone.prob)^2),
  answer.mle2$rb.iter.estimates,
  sum((mle2.mean.abundances - sim.data$true.clone.prob)^2),
  
  answer.loo1$rb.iter.estimates,
  sum((loo1.mean.abundances - sim.data$true.clone.prob)^2),
  answer.loo2$rb.iter.estimates,
  sum((loo2.mean.abundances - sim.data$true.clone.prob)^2),

  answer.cpc1$rb.iter.estimates,
  sum((cpc1.mean.abundances - sim.data$true.clone.prob)^2),
  answer.cpc2$rb.iter.estimates,
  sum((cpc2.mean.abundances - sim.data$true.clone.prob)^2)
  )
names(experiment.values) <- experiment.cols

all.cols <- c(meta.cols, baseline.cols, experiment.cols)
all.values <- c(meta.values, baseline.values, experiment.values)

description <- paste(all.cols, collapse = '/')
values.str <- paste(sprintf("%.8e", all.values), collapse = ',')
print(paste(c(description, values.str), collapse = ','))

