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
  regularization.clones = matrix(0, 0, ncol(read.count.matrix)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{read.count.matrix}{
%%     ~~Describe \code{read.count.matrix} here~~
A matrix of read frequences, where each row corresponds to a distinct clone, and each column corresponds to a particular biological (rather than technical) replicate. All biological replicates should be drawn from the same person. Reads from technical replicates of the same underlying templates should be merged into a single column.
}
  \item{regularization.clones}{
%%     ~~Describe \code{regularization.clones} here~~
A matrix of regularzation (i.e. "virtual") data to be appended to the bottom of the read.count.matrix. For example, this could be the vertical concatentation of an indentiy matrix the size of the number of replicates, with a row of 1s. Deafault value is no regularization (i.e. empty matrix). This is recommended if the data sampling is exceptionally sparse.
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
