\name{lymphclon-package}
\alias{lymphclon-package}
\alias{lymphclon}
\docType{package}
\title{
Estimates the clonality score from replicate of abundances data
}
\description{
There are an enormous number of clones(species) distributed in a highly unequal distribution. The experimenters collect a small number of very coarse samplings, and perform noisy measurements on the samplings, leading to abundance estimations for the indiviudal biological replicates. The goal of this package is to estimate the probability that two random cells belong to the same clone. Clonality estimation arises as a naturally interesting question when trying to understand the diversity of B cell populations, T cell populations, microbial populations, cancer cell subclones, and artificial biologicl libraries, and also broader ecological settings.

This package provides two primary functions; one computes the clonality score estimate (probability that two random individuals belong to the same clone), given replicate of abundances; it also infers the underlying clonal distribution, which may be of interest in diagnostic and scientific settings. 

The other function generates reasonable simulation data, for evaluation purposes. The default behavior of this function generates a very large number of classes, which is suitable for B and T cell settings. It can be easily repurposed for other applications. The simulation assumes random sampling, amplifying, and sequencing of individual cells; this is implemented by introducing sparse sampling of individual cells, followed by introduction of log-normal error, and finally unbiased "rounding" by use of a poisson distribution.

}
\details{
\tabular{ll}{
Package: \tab lymphclon\cr
Type: \tab Package\cr
Version: \tab 1.2.0\cr
Date: \tab 2013-07-15\cr
License: \tab LGPL-2 \cr
}
}
\author{
Author and Maintainer: Yi Liu <liuyipei@stanford.edu>
}
\references{
Accurate Estimation of Clonal Coincidences and Abundances from Biological Replicates: lymphclon. Manuscript submitted.
}
\keyword{ diversity, clonality score, clonality }
\examples{
my.data <- generate.clonal.data(n=2e3) 
# n ~ 2e7 is more appropriate for a realistic B cell repertoire
my.lymphclon.results <- infer.clonality(my.data$read.count.matrix)
# a consistently improved estimate of clonality (the squared 
# 2-norm of the underlying multinomial distribution)
}
