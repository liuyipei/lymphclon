pkgname <- "lymphclon"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('lymphclon')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("generate.clonal.data")
### * generate.clonal.data

flush(stderr()); flush(stdout())

### Name: generate.clonal.data
### Title: generate.clonal.data (part of lymphclon package)
### Aliases: generate.clonal.data
### Keywords: ~kwd1 ~kwd2

### ** Examples

my.data <- generate.clonal.data(n=2e3) 
# n ~ 2e7 is more appropriate for a realistic B cell repertoire
my.lymphclon.results <- infer.clonality(my.data$read.count.matrix)
# a consistently improved estimate of clonality (the squared 
# 2-norm of the underlying multinomial distribution)



cleanEx()
nameEx("infer.clonality")
### * infer.clonality

flush(stderr()); flush(stdout())

### Name: infer.clonality
### Title: infer.clonality (part of lymphclon package)
### Aliases: infer.clonality
### Keywords: ~kwd1 ~kwd2

### ** Examples

my.data <- generate.clonal.data(n=2e3) 
# n ~ 2e7 is more appropriate for a realistic B cell repertoire
my.lymphclon.results <- infer.clonality(my.data$read.count.matrix)
# a consistently improved estimate of clonality (the squared 
# 2-norm of the underlying multinomial distribution)
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

my.data <- generate.clonal.data(n=2e3) 
# n ~ 2e7 is more appropriate for a realistic B cell repertoire
my.lymphclon.results <- infer.clonality(my.data$read.count.matrix)
# a consistently improved estimate of clonality (the squared 
# 2-norm of the underlying multinomial distribution)



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
