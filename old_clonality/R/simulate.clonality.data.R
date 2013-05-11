#!/usr/bin/Rscript
library(expm)
#library(sqldf)
library(MASS)
library(Matrix)
library(VGAM)

simulate.clonality.data<-function(
n = 2e7, # number of clones -- reasonably 20 million
num.cells.taken.vector  = c(2e3,5e3,1e4,2e4,5e4,5e4))
{
unnorm.clone.prob  =  (1:n)^-0.7
true.clone.prob <- unnorm.clone.prob / sum(unnorm.clone.prob)
true.clonality<-sum(true.clone.prob*true.clone.prob)
num.replicates<- length(num.cells.taken.vector)
replicates<-matrix(0,n,num.replicates)

for (i in 1:num.replicates)
{

read.count.per.replicate <- 2e4

#num.cells.taken<-2e4
num.cells.taken <- num.cells.taken.vector[i]

num.obs.clones.frac <- num.cells.taken / n
sample.of.cells <- rmultinom(n=1, size =num.cells.taken,  prob = true.clone.prob)
readcount.sdlog <- 1
get.readcount.given.cellcount<-function(x){ifelse(x>0,sum(abs(rpareto(n=x, location = 1, shape =1))),0)}
raw.sample.of.reads <- sapply(X=sample.of.cells , FUN=get.readcount.given.cellcount)
sample.of.reads <- rpois(
  n=length(sample.of.cells),
  lambda=raw.sample.of.reads * (read.count.per.replicate / sum(raw.sample.of.reads)))
table(sample.of.reads)
# we won't normalize
replicates[,i] <- sample.of.reads # / sum(sample.of.reads)
}
simulated.data<-list(read.count.matrix=replicates, true.clone.prob = true.clone.prob, true.clonality = true.clonality)
return (simulated.data)
}
