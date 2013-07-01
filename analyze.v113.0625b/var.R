library(plyr)

num.iterations <- 4
meta.cols <- c('power', 'replicates', 'clones', 'cells.scaling', 'chao2')
bln.cols <- c('true', 'bln', 'bln.abe2') # abe2: abundance error 2-norm

regularization.method.names <- 
  c('unregularized', 'ue.zr.full', 
  'eq.zr.half', 'ue.zr.half', 'eq.eq.half', 'ue.eq.half',
  'ue.mn.half', 'ue.mn.full', 'ue.mn.js1')
var.meth.names <- c('opt.rer', 'opt.cov', 'fpc.add', 'fpc.max', 'mle.cov', 'corpcor')

experiment.cols <- Reduce(c, 
  lapply(X = var.meth.names,
    FUN = function(curr.var.meth.name)
      {c(paste(curr.var.meth.name, regularization.method.names, sep = '.'),
      paste(curr.var.meth.name, 'abe2', sep = '.'))}
))

all.cols <- c('labels', meta.cols, bln.cols, experiment.cols)

raw <- read.csv('var.cat', header=F)
colnames(raw) <- all.cols
num.table <- raw[apply(raw, 1, function(X){all(!is.na(X))}), -1]

var.settings.names <- c('bln', Reduce(c, 
  lapply(X = var.meth.names,
    FUN = function(curr.var.meth.name)
      {c(paste(curr.var.meth.name, regularization.method.names, sep = '.'))}
)))
var.settings.finaliter.names <- var.settings.names

for (curr.vsn in var.settings.names) 
{
  num.table[, paste(c(curr.vsn, 'err'), collapse = '.')] <- NULL
  num.table[, paste(c(curr.vsn, 'err'), collapse = '.')] <-
    num.table[, curr.vsn] - num.table$true
  num.table[, paste(c(curr.vsn, 'err2'), collapse = '.')] <-
    num.table[, paste(c(curr.vsn, 'err'), collapse = '.')] ^ 2
}

var.abe.names <- c('bln', var.meth.names)

for (curr.vsn in var.settings.names) 
{
  num.table[, paste(c(curr.vsn, 'err'), collapse = '.')] <- NULL
  num.table[, paste(c(curr.vsn, 'err'), collapse = '.')] <-
    num.table[, curr.vsn] - num.table$true
  num.table[, paste(c(curr.vsn, 'err2'), collapse = '.')] <-
    num.table[, paste(c(curr.vsn, 'err'), collapse = '.')] ^ 2
}

mean.err2.table <- ddply(num.table, .(power), 
.fun = function(curr.table){ 
  e2v<-apply(curr.table[, paste(var.settings.names, 'err2', sep = '.')], 2, var)
  names(e2v) <- paste(var.settings.names, 'e2v', sep = '.')
  return(c(
    apply(curr.table[, paste(var.settings.names, 'err', sep = '.')], 2, mean), 
    apply(curr.table[, paste(var.settings.names, 'err2', sep = '.')], 2, mean), 
    apply(curr.table[, paste(var.abe.names, 'abe2', sep = '.')], 2, mean),
    e2v
  ))
  })
  

err2.rat.table <- mean.err2.table$bln.err2 / mean.err2.table[, grep('err2', colnames(mean.err2.table))]
log.err2.rats <- apply(log(err2.rat.table[, grep('err2', colnames(err2.rat.table))]),2,mean)

abe2.rat.table <- mean.err2.table$bln.abe2 / mean.err2.table[, grep('abe2', colnames(mean.err2.table))]
e2v.table <- mean.err2.table[, grep('e2v', colnames(mean.err2.table))]
bias.rat.table <- 
  (mean.err2.table[, grep('err$', colnames(mean.err2.table))] ^ 2) /  
  mean.err2.table[, grep('err2$', colnames(mean.err2.table))]


methods.err2.compare <- data.frame(
geo = exp(apply(log(err2.rat.table),2, mean)),
min = apply(err2.rat.table,2, min),
med = apply(err2.rat.table,2, median)
)

write.csv(x=err2.rat.table, file = 'err2.csv')
library(pheatmap)
pheatmap(err2.rat.table, cluster_rows = F, display_numbers = T, fontsize_number=6)


meth.to.plot <- 'fpc.add.ue.zr.half'
abe.to.plot <- 'fpc.add'

err2.table <- data.frame(
method = c(rep('bln',nrow(num.table)), 
           rep(meth.to.plot,  nrow(num.table))),
err2 = c(num.table$bln.err2, num.table[, paste(c(meth.to.plot, 'err2'), collapse = '.')]),
power = c(num.table$power, num.table$power))

postscript('box.plot.eps')
par(mar=c(5,4,4,2)+2) #should alleviate the margin problem (default is +0.1)
boxplot(err2 ~ power, data = err2.table, 
        boxwex = 0.25, at = 1:21 + 0.2,
        subset = method == "bln", col = "yellow", log = 'y',
        main = "Simulated performance of baseline and our estimator",
        xlab = "Underlying Zipf Power",
        ylab = "Empirical Estimator Squared Error (Log Axis)")
boxplot(err2 ~ power, data = err2.table, add = TRUE,
        boxwex = 0.25, at = 1:21 - 0.2,
        subset = method == meth.to.plot, col = "orange")
legend(x = 'top', c(meth.to.plot, "Unbiased Baseline (Paramsweran 2013)"),
 fill = c("orange", "yellow"))
dev.off()

#plot(mean.err2.table$power, mean.err2.table$bt.err, col ='blue')
#points(mean.err2.table$power, mean.err2.table$bln.err, col = 'red')
#plot(log(mean.err2.table$bln.err), log(mean.err2.table$bt.err))

library(car)
#scatterplot(log(bln.err2) ~ power|sign(bln.err), data = mean.err2.table)
#scatterplot(log(bln.err^2) ~ power|sign(bln.err), data = mean.err2.table)
#scatterplot(log(bt.err^2) ~ power|sign(bt.err), data = mean.err2.table)

postscript('err2.point.plot.eps')
par(mar=c(5,4,4,2)+2) #should alleviate the margin problem (default is +0.1)
numbers.to.plot <- c(mean.err2.table$bln.err2, 
         mean.err2.table[, paste(c(meth.to.plot, 'err2'), collapse = '.')])
ylim <- c(min(numbers.to.plot), max(numbers.to.plot))
plot(mean.err2.table$power, mean.err2.table[, paste(c(meth.to.plot, 'err2'), collapse = '.')], 
  col = 'blue', log = 'y', pch = 1, lwd = 3,
  main = sprintf('Mean squred error compared over %d simulations', nrow(num.table)),
  xlab = 'Underlying Clonal Distribution Zipf Power', 
  ylab = 'Clonality Mean Squared Error (Log Axis)',
  ylim = ylim)
points(mean.err2.table$power, mean.err2.table$bln.err2, col = 'red', pch = 1, lwd = 3)
legend(x = 'top', c("Unbiased Baseline (Paramsweran 2013)", 
       sprintf('Low Variance Method [%s]', meth.to.plot)),
  col = c("red", "blue"), pch = 1)
dev.off()

postscript('abe2.point.plot.eps')
numbers.to.plot <- c(mean.err2.table$bln.abe2, 
         mean.err2.table[, paste(c(abe.to.plot, 'abe2'), collapse = '.')])
ylim <- c(0, 1.2*max(numbers.to.plot))
par(mar=c(5,4,4,2)+2) #should alleviate the margin problem (default is +0.1)
plot(mean.err2.table$power, mean.err2.table[, paste(c(abe.to.plot, 'abe2'), collapse = '.')], 
  col = 'blue', log = '', main = 'Abundance Squared Error Comparisons', pch = 2, ylim = ylim,
  xlab = 'Underlying Clonal Distribution Zipf Power', ylab = 'Abundance Squared Error (Log Axis)',
  )
points(mean.err2.table$power, mean.err2.table$bln.abe2, col = 'red', pch = 2)
legend(x = 'top', c(abe.to.plot, "Baseline"),
  col = c("blue", "red"), pch = 2)
dev.off()
