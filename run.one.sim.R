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
clones <- 5e2 # 500 

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
answer.opt2 <- infer.clonality(read.count.matrix = x, 
  estimate.abundances = T, variance.method = 'usr.2',
  use.squared.err.est = sim.data$grahm.errs)

answer.fpc1 <- infer.clonality(read.count.matrix = x, 
  estimate.abundances = T, variance.method = 'fpc.1', num.iterations = num.iterations)
answer.fpc2 <- infer.clonality(read.count.matrix = x, 
  estimate.abundances = T, variance.method = 'fpc.2', num.iterations = num.iterations)
answer.mle1 <- infer.clonality(read.count.matrix = x, 
  estimate.abundances = T, variance.method = 'mle.1')
answer.mle2 <- infer.clonality(read.count.matrix = x, 
  estimate.abundances = T, variance.method = 'mle.2')

answer.cpc1 <- infer.clonality(read.count.matrix = x, 
  estimate.abundances = T, variance.method = 'corpcor.1')
answer.cpc2 <- infer.clonality(read.count.matrix = x, 
  estimate.abundances = T, variance.method = 'corpcor.2')

answer.idt1 <- infer.clonality(read.count.matrix = x, 
  estimate.abundances = T, variance.method = 'usr.1',
  use.squared.err.est = rep(mean(sim.data$replicate.squared.errs), length(num.cells.taken.vector)))
answer.idt2 <- infer.clonality(read.count.matrix = x, 
  estimate.abundances = T, variance.method = 'usr.2',
  use.squared.err.est = diag(rep(mean(sim.data$replicate.squared.errs), length(num.cells.taken.vector))))

answer.reg1 <- infer.clonality(read.count.matrix = x, 
  estimate.abundances = T, variance.method = 'usr.1',
  use.squared.err.est = sim.data$replicate.squared.errs,
  regularization.method = 'ue.zr.full')
answer.reg2 <- infer.clonality(read.count.matrix = x, 
  estimate.abundances = T, variance.method = 'fpc.1', 
  num.iterations = num.iterations,
  regularization.method = 'ue.zr.full')

answer.reg3 <- infer.clonality(read.count.matrix = x, 
  estimate.abundances = T, variance.method = 'fpc.1', 
  regularization.method = 'eq.zr.half')
answer.reg4 <- infer.clonality(read.count.matrix = x, 
  estimate.abundances = T, variance.method = 'fpc.1', 
  regularization.method = 'ue.zr.half')
answer.reg5 <- infer.clonality(read.count.matrix = x, 
  estimate.abundances = T, variance.method = 'fpc.1', 
  regularization.method = 'eq.eq.half')
answer.reg6 <- infer.clonality(read.count.matrix = x, 
  estimate.abundances = T, variance.method = 'fpc.1', 
  regularization.method = 'ue.eq.half')

answer.rgn1 <- infer.clonality(read.count.matrix = x, 
  estimate.abundances = T, variance.method = 'fpc.1', 
  regularization.method = 'ue.mn.half')
answer.rgn2 <- infer.clonality(read.count.matrix = x, 
  estimate.abundances = T, variance.method = 'fpc.1', 
  regularization.method = 'ue.mn.full')
answer.rgn3 <- infer.clonality(read.count.matrix = x, 
  estimate.abundances = T, variance.method = 'fpc.1', 
  regularization.method = 'ue.mn.js1')


simple.mean.abundances <- apply((x %*% diag(1/apply(x,2,sum)) ), 1,mean)
opt1.mean.abundances    <- as.numeric(answer.opt1$estimated.abundances)
opt2.mean.abundances    <- as.numeric(answer.opt2$estimated.abundances)
fpc1.mean.abundances    <- as.numeric(answer.fpc1$estimated.abundances)
fpc2.mean.abundances    <- as.numeric(answer.fpc2$estimated.abundances)
mle1.mean.abundances    <- as.numeric(answer.mle1$estimated.abundances)
mle2.mean.abundances    <- as.numeric(answer.mle2$estimated.abundances)
cpc1.mean.abundances    <- as.numeric(answer.cpc1$estimated.abundances)
cpc2.mean.abundances    <- as.numeric(answer.cpc2$estimated.abundances)
idt1.mean.abundances    <- as.numeric(answer.idt1$estimated.abundances)
idt2.mean.abundances    <- as.numeric(answer.idt2$estimated.abundances)
reg1.mean.abundances    <- as.numeric(answer.reg1$estimated.abundances)
reg2.mean.abundances    <- as.numeric(answer.reg2$estimated.abundances)
reg3.mean.abundances    <- as.numeric(answer.reg3$estimated.abundances)
reg4.mean.abundances    <- as.numeric(answer.reg4$estimated.abundances)
reg5.mean.abundances    <- as.numeric(answer.reg5$estimated.abundances)
reg6.mean.abundances    <- as.numeric(answer.reg6$estimated.abundances)
rgn1.mean.abundances    <- as.numeric(answer.rgn1$estimated.abundances)
rgn2.mean.abundances    <- as.numeric(answer.rgn2$estimated.abundances)
rgn3.mean.abundances    <- as.numeric(answer.rgn3$estimated.abundances)

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

experiment.cols <- 
  c('opt1', 'opt1.abe2', 'opt2', 'opt2.abe2', 
    paste('fpc1', c(1:num.iterations), sep = '.'), 'fpc1.abe2', 
    paste('fpc2', c(1:num.iterations), sep = '.'), 'fpc2.abe2', 
    'mle1', 'mle1.abe2', 'mle2', 'mle2.abe2',
    'cpc1', 'cpc1.abe2', 'cpc2', 'cpc2.abe2',
    'idt1', 'idt1.abe2', 'idt2', 'idt2.abe2',
    'reg1', 'reg1.abe2', 
    paste('reg2', c(1:num.iterations), sep = '.'), 'reg2.abe2', # fpc.1
    'reg3', 'reg3.abe2', 'reg4', 'reg4.abe2', 
    'reg5', 'reg5.abe2', 'reg6', 'reg6.abe2',
    'rgn1', 'rgn1.abe2', 'rgn2', 'rgn2.abe2',
    'rgn3', 'rgn3.abe2')

experiment.values <- c(
  answer.opt1$rb.iter.estimates,
  sum((opt1.mean.abundances - sim.data$true.clone.prob)^2),
  answer.opt2$rb.iter.estimates,
  sum((opt2.mean.abundances - sim.data$true.clone.prob)^2),

  answer.fpc1$rb.iter.estimates,
  sum((fpc1.mean.abundances - sim.data$true.clone.prob)^2),
  answer.fpc2$rb.iter.estimates,
  sum((fpc2.mean.abundances - sim.data$true.clone.prob)^2),
  
  answer.mle1$rb.iter.estimates,
  sum((mle1.mean.abundances - sim.data$true.clone.prob)^2),
  answer.mle2$rb.iter.estimates,
  sum((mle2.mean.abundances - sim.data$true.clone.prob)^2),

  answer.cpc1$rb.iter.estimates,
  sum((cpc1.mean.abundances - sim.data$true.clone.prob)^2),
  answer.cpc2$rb.iter.estimates,
  sum((cpc2.mean.abundances - sim.data$true.clone.prob)^2),

  answer.idt1$rb.iter.estimates,
  sum((idt1.mean.abundances - sim.data$true.clone.prob)^2),
  answer.idt2$rb.iter.estimates,
  sum((idt2.mean.abundances - sim.data$true.clone.prob)^2),

  answer.reg1$rb.iter.estimates,
  sum((reg1.mean.abundances - sim.data$true.clone.prob)^2),
  answer.reg2$rb.iter.estimates,
  sum((reg2.mean.abundances - sim.data$true.clone.prob)^2),

  answer.reg3$rb.iter.estimates,
  sum((reg3.mean.abundances - sim.data$true.clone.prob)^2),
  answer.reg4$rb.iter.estimates,
  sum((reg4.mean.abundances - sim.data$true.clone.prob)^2),
  
  answer.reg5$rb.iter.estimates,
  sum((reg5.mean.abundances - sim.data$true.clone.prob)^2),
  answer.reg6$rb.iter.estimates,
  sum((reg6.mean.abundances - sim.data$true.clone.prob)^2),

  answer.rgn1$rb.iter.estimates,
  sum((rgn1.mean.abundances - sim.data$true.clone.prob)^2),
  answer.rgn2$rb.iter.estimates,
  sum((rgn2.mean.abundances - sim.data$true.clone.prob)^2),
  answer.rgn3$rb.iter.estimates,
  sum((rgn3.mean.abundances - sim.data$true.clone.prob)^2)
  )
names(experiment.values) <- experiment.cols

all.cols <- c(meta.cols, baseline.cols, experiment.cols)
all.values <- c(meta.values, baseline.values, experiment.values)

description <- paste(all.cols, collapse = '/')
values.str <- paste(sprintf("%.8e", all.values), collapse = ',')
print(paste(c(description, values.str), collapse = ','))

