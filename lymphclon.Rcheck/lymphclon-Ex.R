pkgname <- "lymphclon"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('lymphclon')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("infer.clonality")
### * infer.clonality

flush(stderr()); flush(stdout())

### Name: infer.clonality
### Title: infer.clonality (part of lymphclon package)
### Aliases: infer.clonality
### Keywords: ~kwd1 ~kwd2

### ** Examples

my.data <- simulate.clonality.data(n=2e3) 
# n ~ 2e7 is more appropriate for a realistic B cell repertoire
my.lymphclon.results <- infer.clonality(my.data$read.count.matrix)
# a consistently improved estimate of clonality (the squared 2-norm of the underlying multinomial probabilistic distrubution vector)
my.lymphclon.results$lymphclon.clonality 



cleanEx()
nameEx("lymphclon-package")
### * lymphclon-package

flush(stderr()); flush(stdout())

### Name: lymphclon-package
### Title: Estimates the clonality score from replicate of abundances data
### Aliases: lymphclon-package lymphclon
### Keywords: diversity, clonality score, clonality

### ** Examples

my.data <- simulate.clonality.data(n=2e3) 
# n ~ 2e7 is more appropriate for a realistic B cell repertoire
my.lymphclon.results <- infer.clonality(my.data$read.count.matrix)
# a consistently improved estimate of clonality (the squared 2-norm of the underlying multinomial probabilistic distrubution vector)



cleanEx()
nameEx("simulate.clonality.data")
### * simulate.clonality.data

flush(stderr()); flush(stdout())

### Name: simulate.clonality.data
### Title: simulate.clonality.data (part of lymphclon package)
### Aliases: simulate.clonality.data
### Keywords: ~kwd1 ~kwd2

### ** Examples

my.data <- simulate.clonality.data(n=2e3) 
# n ~ 2e7 is more appropriate for a realistic B cell repertoire
my.lymphclon.results <- infer.clonality(my.data$read.count.matrix)
# a consistently improved estimate of clonality (the squared 2-norm of the underlying multinomial probabilistic distrubution vector)



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
