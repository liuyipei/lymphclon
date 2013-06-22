library(plyr)

num.iterations <- 4
meta.cols <- c('power', 'replicates', 'clones', 'cells.scaling')
bln.cols <- c('true', 'bln', 'bln.abe2') # abe2: abundance error 2-norm

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
    'rgn1', 'rgn1.abe2', 'rgn2', 'rgn2.abe2', 'rgn3', 'rgn3.abe2'
    )
all.cols <- c('labels', meta.cols, bln.cols, experiment.cols)

raw <- read.csv('var.cat', header=F)
colnames(raw) <- all.cols
num.table <- raw[apply(raw, 1, function(X){all(!is.na(X))}), -1]

var.settings.names <- c('bln', 'opt1', 'opt2',
    paste('fpc1', c(1:num.iterations), sep = '.'), 
    paste('fpc2', c(1:num.iterations), sep = '.'), 
    'mle1', 'mle2', 'cpc1', 'cpc2', 
    'idt1', 'idt2', 'reg1',
    paste('reg2', c(1:num.iterations), sep = '.'),
    'reg3', 'reg4', 'reg5','reg6',
    'rgn1', 'rgn2', 'rgn3')

var.settings.finaliter.names <- c('bln', 'opt1', 'opt2',
    paste('fpc1', num.iterations, sep = '.'), 
    paste('fpc2', num.iterations, sep = '.'), 
    'mle1', 'mle2', 'cpc1', 'cpc2',
    'idt1', 'idt2', 'reg1',
    paste('reg2', num.iterations, sep = '.'),
    'reg3', 'reg4', 'reg5','reg6',
    'rgn1', 'rgn2', 'rgn3')

for (curr.vsn in var.settings.names) 
{
  num.table[, paste(c(curr.vsn, 'err'), collapse = '.')] <- NULL
  num.table[, paste(c(curr.vsn, 'err'), collapse = '.')] <-
    num.table[, curr.vsn] - num.table$true
  num.table[, paste(c(curr.vsn, 'err2'), collapse = '.')] <-
    num.table[, paste(c(curr.vsn, 'err'), collapse = '.')] ^ 2
}

var.abe.names <- c('bln', 
    'opt1', 'opt2',
    'fpc1', 'fpc2',
    'mle1', 'mle2',
    'cpc1', 'cpc2',
    'idt1', 'idt2',
    'reg1', 'reg2',
    'reg3', 'reg4',
    'reg5', 'reg6',
    'rgn1', 'rgn2', 'rgn3'
    )


for (curr.vsn in var.settings.names) 
{
  num.table[, paste(c(curr.vsn, 'err'), collapse = '.')] <- NULL
  num.table[, paste(c(curr.vsn, 'err'), collapse = '.')] <-
    num.table[, curr.vsn] - num.table$true
  num.table[, paste(c(curr.vsn, 'err2'), collapse = '.')] <-
    num.table[, paste(c(curr.vsn, 'err'), collapse = '.')] ^ 2
}


mean.err2.table <- ddply(num.table, .(power), 
.fun = function(curr.table){ c(
  apply(curr.table[, paste(var.settings.names, 'err', sep = '.')], 2, mean), 
  apply(curr.table[, paste(var.settings.names, 'err2', sep = '.')], 2, mean), 
  apply(curr.table[, paste(var.abe.names, 'abe2', sep = '.')], 2, mean),
  sapply(paste(var.settings.names, 'err', sep = '.'), 
    function(curr.errname){
      var.test(curr.table[, 'bln.err'], curr.table[, curr.errname])$statistic
    })
  )})
  

err2.rat.table <- mean.err2.table$bln.err2 / mean.err2.table[, grep('err2', colnames(mean.err2.table))]
log.err2.rats <- apply(log(err2.rat.table[, grep('err2', colnames(err2.rat.table))]),2,mean)
F.rats <- mean.err2.table[, grep('.F', colnames(mean.err2.table))]

abe2.rat.table <- mean.err2.table$bln.abe2 / mean.err2.table[, grep('abe2', colnames(mean.err2.table))]

bias.rat.table <- 
  (mean.err2.table[, grep('err$', colnames(mean.err2.table))] ^ 2) /  
  mean.err2.table[, grep('err2$', colnames(mean.err2.table))]

meth.to.plot <- 'reg4'
abe.to.plot <- 'reg4'

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
par(mar=c(5,4,4,2)+2) #should alleviate the margin problem (default is +0.1)
plot(mean.err2.table$power, mean.err2.table[, paste(c(abe.to.plot, 'abe2'), collapse = '.')], 
  col = 'blue', log = 'y', main = 'Abundance Squared Error Comparisons', pch = 2, ylim = c(1e-4,1e-1),
  xlab = 'Underlying Clonal Distribution Zipf Power', ylab = 'Abundance Squared Error (Log Axis)')
points(mean.err2.table$power, mean.err2.table$bln.abe2, col = 'red', pch = 2)
legend(x = 'top', c(abe.to.plot, "Baseline"),
  col = c("blue", "red"), pch = 2)
dev.off()


methods.err2.compare <- data.frame(
geo = exp(apply(log(err2.rat.table),2, mean)),
min = apply(err2.rat.table,2, min),
med = apply(err2.rat.table,2, median)
)
