source('simulate.clonality.data.R')
source('infer.clonality.R')

num.replicates <- 8
num.cells.taken.vector <- c(1e2, 1e2, 1e2, 1e2, 1e2, 1e2, 1e3, 1e3)
read.count.per.replicate.vector <- c(5e2, 5e2, 5e2, 5e2, 5e2, 5e2, 5e2, 5e2)

simulate.clonality.data(n = 100, num.cells.taken.vector = num.cells.taken.vector,
  read.count.per.replicate.vector = read.count.per.replicate.vector, 
  clonal.distribution.power = -0.9) -> sim.data

infer.clonality(sim.data$read.count.matrix, variance.method = 'loo.1', num.iterations = 20)->answer

data.frame(
err = c(answer$simple.precision.clonality, answer$rb.iter.estimates) - sim.data$true.clonality,
delta = diff(c(0, answer$simple.precision.clonality, answer$rb.iter.estimates))
)
sim.data$true.clonality #0.04481423