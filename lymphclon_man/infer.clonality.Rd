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
  variance.method = 'fpc.1',
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
The method is a code defaulting to "loo.2".
usr.1: argument use.squared.err.est specifies 1st order variances
mle.1: use maximum likelihood estimate of 1st order variances
loo.1: use leave-one-out estimate of 1st order variances
corpcor.1: use shrinkage-based estimate of 1st order variances (see package corpcor)
}
\item{estimate.abundances}{
A boolean value defaulting to false. If set to true, then the return value will be in the form of a list, and include an additional item named "estiamted.abundances". The estimated abundances are computed as the (conditional) precision weighted average of the read.count.matrix.
}
\item{num.iterations}{
An integer specifying the number of iterations in estimation; applicable only if the variance method is set to one of the loo or mle methods. Note that the ue.zr.half setting provides the most consistent improvement against the naive simple clonality score, in terms of MSE reduction.
}
\item{internal.parameters}{
A named list of internal parameters used primarily to speed up evaluation of MSE performance across a large number of simulations.
The entries are named "replicates", "rep.grahm.matrix", "simple.precision.clonality", and "use.squared.err.est".
"replicates" refer to a column normalized (to 1) matrix of read.count.matrix.
"rep.grahm.matrix" refers to t(replicates) %*% replicates.
"simple.precision.clonality" refers to the unweighted averaging of all n-choose-2 estimators of clonality.
"use.squared.err.est" refers to a user-provided set of squared error estimates corresponding to each of the replicates; this parameter is active iff  the variance method is set to 'usr.1'.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{simple.precision.clonality}{This is a base line estimator of the clonality score based on simple modeling. Briefly, for each pairing of reads possible, where the two reads arise from different biological replicates -- this pairing is considered as an observation from a single share undelrying bernoulli distribution with parameter equal to the clonality score. This estimator computes maximum likelihood estimate of the underlying clonality score. This estimator can be thought of as a baseline estimate. Unlike the Gini-Simpson-Estimator, this baseline does not suffer from any convexity-based bias.}
\item{lymphclon.clonality}{This is an estimator is a weighted average of n-choose-2 individually unbiased estimators that arise from comparing the provided n biological replicates pairwise. The weighting is determined by the covariance between these estimators. This has much lower variance than the simple.precision.clonality. This is the value returned by the ue.zr.half regularization setting}
\item{estimated.abundances}{
If the estimate.abundances parameter was set to True, then return a probability vector indicating the fractional contributions of each individual clone to the overall multinomial repertoire.
}
\item{esimated.squared.errs}{
A numerical measure of the estimated squared 2-norm error of each given replicate.
}
\item{estimated.precisions}{
A numerical measure of the estimated precisions of each given replicate.
}
\item{variance.method}{
The method code provided in the input, determining the method by which the replicate-level covariances are computed.
}
\item{fpc.iter.estimates}{
Specifies the intermediate clonality score estimates across iterations, when num.iterations is greater than 1. Note that, in empirical simulations, it is best to use a num.iterations value of 1. Additional iterations occasionally, but usually do not help in improving the mean squared error of the clonality measures.
}
\item{regularized.estimates}{
Under a variety of regularization settings on the n-choose-2 by n-choose-2 covariance matrix for the n-choose-2 individually unbiased clonality estimators, return the respectively weighted clonality measures for each regularization setting. 
}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
my.data <- simulate.clonality.data(n=2e3) 
# n ~ 2e7 is more appropriate for a realistic B cell repertoire
my.lymphclon.results <- infer.clonality(my.data$read.count.matrix)
# a consistently improved estimate of clonality (the squared 2-norm of the underlying multinomial probabilistic distrubution vector)
my.lymphclon.results$lymphclon.clonality 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
