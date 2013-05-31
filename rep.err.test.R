
source('simulate.clonality.data.R')
source('infer.clonality.R')

num.replicates <- 8
num.cells.taken.vector <- c(1e2, 1e2, 1e2, 1e2, 1e2, 1e2, 1e3, 1e3)
read.count.per.replicate.vector <- c(1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3)

get.err.est <- function(x) {

simulate.clonality.data(n = 10000, num.cells.taken.vector = num.cells.taken.vector,
  read.count.per.replicate.vector = read.count.per.replicate.vector) -> sim.data
results <- infer.clonality(sim.data$read.count.matrix, estimate.abundances = T)

correct.weights <- sim.data$replicate.squared.errs / sum(sim.data$replicate.squared.errs)
estimated.weights <- results$estimated.squared.errs / sum(results$estimated.squared.errs)

loo.weights <- rep(0, num.replicates)
for (i in (1:num.replicates)) {
  loo.weights[-i] <- infer.clonality(sim.data$read.count.matrix[, -i], estimate.abundances = T)$estimated.squared.errs
    + loo.weights[-i]
}
loo.weights <- loo.weights / sum(loo.weights)

replicate.size.weights <- data.frame(
  num.cells = num.cells.taken.vector,
  correct.weights,
  estimated.weights,
  loo.weights
)
}

master.err.list <- lapply(1:1000, get.err.est)
master.err.table <- Reduce(rbind, master.err.list)
master.err.table$j1 <- (num.replicates-1) * master.err.table$estimated.weights 
  - (num.replicates - 2) * master.err.table$loo.weights

library(car)
scatterplot(correct.weights~estimated.weights| num.cells, data = master.err.table)

err.table.100 <- master.err.table[master.err.table$num.cells == 100, ]
err.table.1k  <- master.err.table[master.err.table$num.cells == 1000, ]
library(glmnet)
reg.100 <- stats::glm(correct.weights~estimated.weights+loo.weights-1, data = err.table.100)
reg.1k <- stats::glm(correct.weights~estimated.weights+loo.weights-1, data = err.table.1k)


reg.1k <- stats::glm(correct.weights~estimated.weights-1, data = err.table.1k)
reg.100 <- stats::glm(correct.weights~estimated.weights-1, data = err.table.100)

reg.1k <- stats::glm(correct.weights~loo.weights-1, data = err.table.1k)
reg.100 <- stats::glm(correct.weights~loo.weights-1, data = err.table.100)

reg.1k <- stats::glm(correct.weights~j1-1, data = err.table.1k)
reg.100 <- stats::glm(correct.weights~j1-1, data = err.table.100)