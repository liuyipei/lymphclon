\name{infer.clonality}
\alias{infer.clonality}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
infer.clonality (part of lymphclon package)
}
\description{
Clonality score, a useful metric used in immunology, refers to the probability that two random lymphocyte receptor chain reads drawn with replacement (makes no difference in the immunology context) from an individual corresponded to the same clone, within some given repertoire of either B cells or T cells.  This package implements ane estimator which understands the multi-replicate-with-PCR structure of thes sequencing experiments.
}
\usage{
infer.clonality(
  read.count.matrix, 
  variance.method = 'fpc.max',
  estimate.abundances = F,
  num.iterations = 1,
  internal.parameters = list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{read.count.matrix}{
%%     ~~Describe \code{read.count.matrix} here~~
A matrix of read frequences, where each row corresponds to a distinct clone, and each column corresponds to a particular biological (rather than technical) replicate. All biological replicates should be drawn from the same person. Reads from technical replicates of the same underlying templates should be merged into a single column.
}
\item{variance.method}{
The method is a code defaulting to "fpc.max". This code determines how the covariance between the replicates -- the replicates -- an n by n matrix, is estimated.  Other possible codes are "fpc.add", "mle.cov", and "corpcor". Empirically, from simulations ,"fpc.max" performs the best. "mle.cov" uses the empirical covariance between the replicates directly. "fpc.max" and "fpc.add" recognize that the individual off diagonal entries of this matrix should all be the same, and that they can be estimated with a simple unbiased clonality score estimate (we use the simple estimate here), but the diagonal entries differ: "fpc.max". The diagonal entries should be larger than the off diagonal entries. When the empirical variance of a replicate is lower than the off diagonal entries, "fpc.max" uses the diagonal entries' value; in contrast, "fpc.add" penalizes the replicate by adding the difference onto the off diagonal entries' values. Empirically, "fpc.max" performs best, but the margin from "fpc.add" is very very small. "corpcor" refers to a very strongly regularized method of covariance estimation from the corpcor package. This is not recommended because the shrinkage is so strong that the useful signal is lost, and performance falls to levels produced by the simple clonality estimate. 
}
\item{estimate.abundances}{
A boolean value defaulting to false. If set to true, then the return value will be in the form of a list, and include an additional item named "estiamted.abundances". The estimated abundances are computed as the (conditional) precision weighted average of the read.count.matrix.
}
\item{num.iterations}{
An integer specifying the number of iterations in estimation; applicable only if the variance method is set to one of the loo or mle methods. Note that the ue.zr.half setting provides the most consistent improvement against the naive simple clonality score, in terms of MSE reduction.
}
\item{internal.parameters}{
A named list of internal parameters used primarily to speed up evaluation of MSE performance across a large number of simulations. Another use is to provide user-specified covariances internally into the method, for debugging and evaluation purposes; the typical user should ignore this parameter.
The entries are named "replicates", "rep.grahm.matrix", "simple.precision.clonality", and "use.squared.err.est".
"replicates" refer to a column normalized (to 1) matrix of read.count.matrix.
"rep.grahm.matrix" refers to t(replicates) %*% replicates.
"simple.precision.clonality" refers to the unweighted averaging of all n-choose-2 estimators of clonality.
"use.squared.err.est" refers to a user-provided set of squared error estimates corresponding to each of the replicates; this parameter is active iff  the variance method is set to 'usr.1'.
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{estimated.abundances}{
If the estimate.abundances parameter was set to True, then return a probability vector indicating the fractional contributions of each individual clone to the overall multinomial repertoire.
}
\item{d1jkn.covariance}{
Lymphclon performs a number of regularization schemes to the between comparison positive definite covariance matrix, with n-choose-2 rows and columns, where n is the number of replicates. These regularizatino schemes are jackknifed to determine their variances and covariances, so that they can be averaged. This is the covariance matrix, as determined by a leave-one-replicate out jackknife. This is return value is populated only if there are at least 4 replicates.
}
\item{estimated.squared.errs}{
A numerical measure of the estimated squared 2-norm error of each given replicate.
}
\item{estimated.precisions}{
A numerical measure of the estimated precisions of each given replicate.
}
\item{variance.method}{
The method code provided in the input, determining the method by which the replicate-level covariances are computed.
}
\item{fpc.iter.estimates}{
Specifies the intermediate clonality score estimates across iterations, when num.iterations is greater than 1. Note that, based on empirical simulations, it is best to use a num.iterations value of 1. Additional iterations occasionally, but usually do not help in improving the mean squared error of the clonality measures.
}
\item{simple.precision.clonality}{
This is a base line estimator of the clonality score based on simple modeling. Briefly, for each pairing of reads possible, where the two reads arise from different biological replicates -- this pairing is considered as an observation from a single share undelrying bernoulli distribution with parameter equal to the clonality score. This estimator computes maximum likelihood estimate of the underlying clonality score. This estimator can be thought of as a baseline estimate. Unlike the Gini-Simpson-Estimator, this baseline does not suffer from any convexity-based bias.
}
\item{regularized.estimates}{
Under a variety of regularization settings on the n-choose-2 by n-choose-2 covariance matrix for the n-choose-2 individually unbiased clonality estimators, return the respectively weighted clonality measures for each regularization setting. 
}
\item{mixture.estimates}{
There are several regularization schemes possible on the n-choose-2 by n-choose-2 covariance matrix between the n-choose pair-replicate comparisons. This estimate of clonality performs a jackknife to estimate the covariance bewteen estimates corresponding to each regularization scheme, and then performs a MLE averaging based on the covariance and the individual estimates. This is the best estimate of the clonality score provided by Lymphclon; however, it requires at least 4 replicates. There are several ways by which the selected regularized estimates can be averaged: by their jackknife covariance matrix, by their jackknife variances (ie. diagonal covariance), and simple equal weights.
}
\item{mixture.clonality}{
This defaults to the matrix weighted mixture estimate, unless it falls outside the range of the contributing regularized estimates; in which case, this value defaults to the scalar jackknife precision weighted average.
}
\item{lymphclon.clonality}{
This is the recommended clonality return value to use. This defaults to mixture.clonality;  if there are only 3 replicates, then defaults to the estimate provided by the "ue.zr.half" (diagonal elements are maintained, but off diagonals are divided by 2) regularization setting on the n-choose-2 by n-choose-2 between-comparison covariance matrix.
}
}
\author{
Yi Liu (liuyipei@stanford.edu / liu.yi.pei@gmail.com)
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
my.data <- generate.clonal.data(n=2e3) 
# n ~ 2e7 is more appropriate for a realistic B cell repertoire
my.lymphclon.results <- infer.clonality(my.data$read.count.matrix)
# a consistently improved estimate of clonality (the squared 
# 2-norm of the underlying multinomial distribution)
my.lymphclon.results$lymphclon.clonality 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
