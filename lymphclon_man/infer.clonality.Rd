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
  variance.method = 'loo.2',
  estimate.abundances = F,
  loo.squared.err.est = c(),
  use.squared.err.est = c(),
  num.iterations = 1)
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
usr.2: argument use.squared.err.est specifies 2nd order variances
mle.1: use maximum likelihood estimate of 1st order variances
mle.2: use maximum likelihood estimate of 2nd order variances
loo.1: use leave-one-out estimate of 1st order variances
loo.2: use leave-one-out estimate of 2nd order variances
corpcor.1: use shrinkage-based estimate of 1st order variances (see package corpcor)
corpcor.2: use shrinkage-based estimate of 2nd order variances (see package corpcor)
}
\item{estimate.abundances}{
A boolean value defaulting to false. If set to true, then the return value will be in the form of a list, and include an additional item named "estiamted.abundances". The estimated abundances are computed as the (conditional) precision weighted average of the read.count.matrix.
}
\item{loo.squared.err.est}{
A boolean value defaulting to true. When estimating the squared errors and precisions associated with each replicate, we can use the maximum likelihood estimate (i.e. set loo.squared.err.est to false), or use an averaged leave-one-out estimate (i.e. set loo.squared.err.est to true). The leave-one-out estimate is slightly biased, and slightly slower to compute, but much less variance. Especially if there are at least three replicates, this value probably should be set to true.
}
\item{use.squared.err.est}{
A vector defaulting to the empty vector. If a non-empty vector is provided, then this forces infer.clonality to use these values as the conditional variances of the respective replicates. If non-empty, then this vector should be positive, and of length equal to the number of replicates.
}
\item{num.iterations}{
An integer specifying the number of iterations in estimation; applicable only if the variance method is set to one of the loo or mle methods.
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
\item{rao.blackwell.mvg.clonality}{This is an estimator is a weighted average of n-choose-2 individually unbiased estimators that arise from comparing the provided n biological replicates pairwise. The weighting is determined by the covariance between these estimators. This has much lower variance than the simple.precision.clonality.}
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
my.data<-simulate.clonality.data(n=2e3) 
# n ~ 2e7 is more appropriate for a realistic B cell repertoire
infer.clonality(my.data$read.count.matrix)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
