\name{generate.clonal.data}
\alias{generate.clonal.data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
generate.clonal.data (part of lymphclon package)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function generates simulated data, for evaluation purposes. We start with an underlying multinomial population with entries proportional to rank^-power distribution, where power is fixed. Next, we draw multinomially from this distribution, a fixed number of cells, to generate each desired replicate. Then, this distribution is subject to log-normal error, and subsequently scaled up to the expected number of reads. To round into integers, the expected number of reads for each clone is finally pushed through a poisson process to generate integer read counts. The poisson "rounding" process is why the resulting read counts are not exactly as specified.
}
\usage{
generate.clonal.data(
  n = 2e+07, 
  num.cells.taken.vector = c(2000, 5000, 10000, 20000, 50000, 50000), 
  read.count.per.replicate.vector = rep(20000, length(num.cells.taken.vector)), 
  clonal.distribution.power = -sqrt(2))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
The true number of distinct clones in the underlying assemblage
}
  \item{num.cells.taken.vector}{
A vector specifying the number of cells taken in each independent biological replicate
}
  \item{read.count.per.replicate.vector}{
A vector of the same length as num.cells.taken.vector, specifying the number of reads generated from each biological replicate, of the same corresponding indices
}
  \item{clonal.distribution.power}{
The true underlying clonal multinomial distribution is proportional to (1:n)^-clonal.distribution.power
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{read.count.matrix}{This is a matrix of simulated counts, with rows corresponding to clones (classes, or species), and columns corresponding to biological replicates}
\item{true.clone.prob}{This is the underlying simulated assemblage multinomial distribution used to generate read.count.matrix}
\item{true.clonality}{This is the true clonality score of the underlying simulated assemblage}
}
\author{
%%  ~~who you are~~
Yi Liu (liuyipei@stanford.edu / liu.yi.pei@gmail.com)
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
my.data <- generate.clonal.data(n=2e3) 
# n ~ 2e7 is more appropriate for a realistic B cell repertoire
my.lymphclon.results <- infer.clonality(my.data$read.count.matrix)
# a consistently improved estimate of clonality (the squared 
# 2-norm of the underlying multinomial distribution)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
