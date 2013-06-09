source('simulate.clonality.data.R')
source('infer.clonality.R')

num.replicates <- 8
num.cells.taken.vector <- c(1e2, 1e2, 1e2, 1e2, 1e2, 1e2, 1e3, 1e3)
read.count.per.replicate.vector <- c(5e2, 5e2, 5e2, 5e2, 5e2, 5e2, 5e2, 5e2)

get.err.est <- function(x) {

simulate.clonality.data(n = 1000, num.cells.taken.vector = num.cells.taken.vector,
  read.count.per.replicate.vector = read.count.per.replicate.vector, 
  clonal.distribution.power = -0.7) -> sim.data

vm.vec <- c('mle.1', 'mle.2', 'loo.1', 'loo.2', 'corpcor.1', 'corpcor.2')
results <- rep(0, length(vm.vec))
for (i in 1:length(vm.vec))
{
  vm <- vm.vec[i]
  curr.results <- infer.clonality(sim.data$read.count.matrix, variance.method = vm)
  print(curr.results$variance.method)
  results[i] <- curr.results$rao.blackwell.mvg.clonality - sim.data$true.clonality
  names(results) <- vm.vec
}
return(c(curr.results$simple.precision.clonality, results))
}

compare.table <- sapply(1:50, get.err.est)
apply(compare.table^2, 1, sum)
