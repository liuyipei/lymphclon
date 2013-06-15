library(plyr)

num.iterations <- 4
meta.cols <- c('power', 'replicates', 'clones', 'cells.scaling')
bln.cols <- c('true', 'bln', 'bln.abe2') # abe2: abundance error 2-norm


experiment.cols <- 
  c('opt1', 'opt1.abe2', 'opt2', 'opt2.abe2', 
    paste('fpc1', c(1:num.iterations), sep = '.'), 'fpc1.abe2', 
    paste('fpc2', c(1:num.iterations), sep = '.'), 'fpc2.abe2', 
    'loo1', 'loo1.abe2', 'loo2', 'loo2.abe2',
    'mle1', 'mle1.abe2', 'mle2', 'mle2.abe2',
    'cpc1', 'cpc1.abe2', 'cpc2', 'cpc2.abe2',
    'idt1', 'idt1.abe2', 'idt2', 'idt2.abe2',
    'reg1', 'reg1.abe2', 
    paste('reg2', c(1:num.iterations), sep = '.'), 'reg2.abe2' # fpc.1
    )
all.cols <- c('labels', meta.cols, bln.cols, experiment.cols)

raw <- read.csv('var.cat')
colnames(raw) <- all.cols
num.table <- raw[, -1]

var.settings.names <- c('bln', 'opt1', 'opt2',
    paste('fpc1', c(1:num.iterations), sep = '.'), 
    paste('fpc2', c(1:num.iterations), sep = '.'), 
    'loo1', 'loo2', 'mle1', 'mle2', 'cpc1', 'cpc2', 
    'idt1', 'idt2', 'reg1',
    paste('reg2', c(1:num.iterations), sep = '.'))

var.settings.finaliter.names <- c('bln', 'opt1', 'opt2',
    paste('fpc1', num.iterations, sep = '.'), 
    paste('fpc2', num.iterations, sep = '.'), 
    'loo1', 'loo2', 'mle1', 'mle2', 'cpc1', 'cpc2',
    'idt1', 'idt2', 'reg1',
    paste('reg2', num.iterations, sep = '.'))

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
    'loo1', 'loo2',
    'mle1', 'mle2',
    'cpc1', 'cpc2',
    'idt1', 'idt2',
    'reg1', 'reg2')


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

meth.to.plot <- 'fpc1.4'
abe.to.plot <- 'fpc1'

err2.table <- data.frame(
method = c(rep('bln',nrow(num.table)), 
           rep(meth.to.plot,  nrow(num.table))),
err2 = c(num.table$bln.err2, num.table[, paste(c(meth.to.plot, 'err2'), collapse = '.')]),
power = c(num.table$power, num.table$power))

png('box.plot.png')
boxplot(err2 ~ power, data = err2.table, 
        boxwex = 0.25, at = 1:21 + 0.2,
        subset = method == "bln", col = "yellow", log = 'y',
        main = "Simulated performance of bln and our estimator (Liu2013)",
        xlab = "Underlying Zipf Power",
        ylab = "Empirical Estimator Squared Error (Log Axis)")
boxplot(err2 ~ power, data = err2.table, add = TRUE,
        boxwex = 0.25, at = 1:21 - 0.2,
        subset = method == meth.to.plot, col = "orange")
legend(x = 'top', c(meth.to.plot, "Baseline"),
 fill = c("orange", "yellow"))
dev.off()

#plot(mean.err2.table$power, mean.err2.table$bt.err, col ='blue')
#points(mean.err2.table$power, mean.err2.table$bln.err, col = 'red')
#plot(log(mean.err2.table$bln.err), log(mean.err2.table$bt.err))

library(car)
#scatterplot(log(bln.err2) ~ power|sign(bln.err), data = mean.err2.table)
#scatterplot(log(bln.err^2) ~ power|sign(bln.err), data = mean.err2.table)
#scatterplot(log(bt.err^2) ~ power|sign(bt.err), data = mean.err2.table)

png('err2.point.plot.png')
plot(mean.err2.table$power, mean.err2.table[, paste(c(meth.to.plot, 'err2'), collapse = '.')], 
  col = 'blue', log = 'y', main = 'Mean Squared Comparisons', pch = 1)
points(mean.err2.table$power, mean.err2.table$bln.err2, col = 'red', pch = 1)
legend(x = 'top', c(meth.to.plot, "Baseline"),
  col = c("blue", "red"), pch = 1)
dev.off()

png('abe2.point.plot.png')
plot(mean.err2.table$power, mean.err2.table[, paste(c(abe.to.plot, 'abe2'), collapse = '.')], 
  col = 'blue', log = 'y', main = 'Abundance Squared Error Comparisons', pch = 2, ylim = c(1e-4,1e-1))
points(mean.err2.table$power, mean.err2.table$bln.abe2, col = 'red', pch = 2)
legend(x = 'top', c(abe.to.plot, "Baseline"),
  col = c("blue", "red"), pch = 2)
dev.off()
