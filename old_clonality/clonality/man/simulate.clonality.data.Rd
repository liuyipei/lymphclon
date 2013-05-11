\name{simulate.clonality.data}
\alias{simulate.clonality.data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Simulates Clonality Data
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function generates simulated data, for evaluation purposes. We start with an underlying multinomial population with entries proportional to rank^-power distribution, where power is fixed. Next, we draw multinomially from this distribution, a fixed number of cells, to generate each desired replicate. Then, this distribution is subject to log-normal error, and subsequently scaled up to the expected number of reads. To round into integers, the expected number of reads for each clone is finally pushed through a poisson process to generate integer read counts. The poisson "rounding" process is why the resulting read counts are not exactly as specified.
}
\usage{
simulate.clonality.data(n = 2e+07, num.cells.taken.vector = c(2000, 5000, 10000, 20000, 50000, 50000))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{num.cells.taken.vector}{
%%     ~~Describe \code{num.cells.taken.vector} here~~
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (n = 2e+07, num.cells.taken.vector = c(2000, 5000, 10000, 
    20000, 50000, 50000)) 
{
    unnorm.clone.prob = (1:n)^-0.7
    true.clone.prob <- unnorm.clone.prob/sum(unnorm.clone.prob)
    true.clonality <- sum(true.clone.prob * true.clone.prob)
    num.replicates <- length(num.cells.taken.vector)
    replicates <- matrix(0, n, num.replicates)
    for (i in 1:num.replicates) {
        read.count.per.replicate <- 20000
        num.cells.taken <- num.cells.taken.vector[i]
        num.obs.clones.frac <- num.cells.taken/n
        sample.of.cells <- rmultinom(n = 1, size = num.cells.taken, 
            prob = true.clone.prob)
        readcount.sdlog <- 1
        get.readcount.given.cellcount <- function(x) {
            ifelse(x > 0, sum(abs(rpareto(n = x, location = 1, 
                shape = 1))), 0)
        }
        raw.sample.of.reads <- sapply(X = sample.of.cells, FUN = get.readcount.given.cellcount)
        sample.of.reads <- rpois(n = length(sample.of.cells), 
            lambda = raw.sample.of.reads * (read.count.per.replicate/sum(raw.sample.of.reads)))
        table(sample.of.reads)
        replicates[, i] <- sample.of.reads
    }
    simulated.data <- list(read.count.matrix = replicates, true.clone.prob = true.clone.prob, 
        true.clonality = true.clonality)
    return(simulated.data)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
