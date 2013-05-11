\name{lymphclon-package}
\alias{lymphclon-package}
\alias{lymphclon}
\docType{package}
\title{
Estimates the clonality score from replicate of abundances data
}
\description{
There are an enormous number of clones(species) distributed in a highly unequal distribution. The experimenters collect a small number of very coarse samplings, and perform noisy measurements on the samplings, leading to abundance estimations for the indiviudal biological replicates. The goal of this package is to estimate the probability that two random cells belong to the same clone. 

It provides two functions; one computes the estimate, given replicate of abundances; the other function generates reasonable simulation data (for B or T cell lymphocyte setting), for evaluation purposes.

This package was written for clonality quantification in the lymphocyte repertoires, as measured by blood-draw and PCR/sequencing; however it can be directly reused for broader ecological applications.
}
\details{
\tabular{ll}{
Package: \tab lymphclon\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-04-24\cr
License: \tab LGPL-2 \cr
}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
This package provides two functions; one computes the estimate, given replicate of abundances; the other function generates reasonable simulation data (for B or T cell lymphocyte setting), for evaluation purposes.
infer.clonality: infers the clonality score, given replicate of abundances data
simulate.clonality.data: generates simulated data, primarily for evaluation purposes
}
\author{
Author and Maintainer: Yi Liu <liuyipei@stanford.edu>
}
\references{
~~ Literature or other references for background information ~~
}
\keyword{ diversity, clonality score, clonality }
\seealso{
\code{\link[<vegan>:<vegan>-package]{<vegan>}}
}
\examples{
my.data<-simulate.clonality.data(n=2e3) 
# n ~ 2e7 is more appropriate for a realistic B cell repertoire
infer.clonality(my.data$read.count.matrix)
}
