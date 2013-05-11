\name{infer.clonality}
\alias{infer.clonality}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Clonality Score Estimation
}
\description{
    Clonality score, a useful metric used in immunology, refers to the probability that two random lymphocyte receptor chain reads drawn with replacement (makes no difference in the immunology context) from an individual corresponded to the same clone, within some given repertoire of either B cells or T cells.  This package implements ane estimator which understands the multi-replicate-with-PCR structure of thes sequencing experiments.

    Names of analogous scores other area include: Simpson's Index (Ecology), Coincidence Index (Cryptogrpahy), Herfindahl-Hirschman Index (Economics), 1 - Gini Index (Economics)
}
\usage{
infer.clonality(read.count.matrix, regularization.clones = matrix(0, 0, ncol(read.count.matrix)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{read.count.matrix}{
%%     ~~Describe \code{read.count.matrix} here~~
   A matrix of read frequences, where each row corresponds to a distinct clone, and each column corresponds to a particular biological (rather than technical) replicate. All biological replicates should be drawn from the same person. Reads from technical replicates of the same underlying templates should be merged into a single column.

   It is also ok to provide within-replicate contributions of each clone, though it would create inconveniences for you if you wish to regularize in a meaningful manner.
}
  \item{regularization.clones}{
%%     ~~Describe \code{regularization.clones} here~~
A matrix of regularzation (i.e. "virtual") data to be appended to the bottom of the read.count.matrix. For example, this could be the vertical concatentation of an indentiy matrix the size of the number of replicates, with a row of 1s. Deafault value is no regularization (i.e. empty matrix). This is recommended especially if there are practically no common collisions between replicates.
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
  my.data<-simulate.clonality.data(2e6) # 2 million is on the very low end for B cells
  infer.clonality(my.data$read.count.matrix)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (read.count.matrix, regularization.clones = matrix(0, 
    0, ncol(read.count.matrix))) 
{
    replicates <- read.count.matrix
    n <- nrow(replicates)
    num.replicates <- ncol(replicates)
    for (i in 1:num.replicates) {
        replicates[, i] <- read.count.matrix[, i]/sum(read.count.matrix[, 
            i])
    }
    clonality.matrix <- t(replicates) \%*\% replicates
    num.estimators <- num.replicates * (num.replicates - 1)/2
    reads.per.replicate <- as.numeric(apply(read.count.matrix, 
        2, sum))
    simple.precision.weights <- matrix(reads.per.replicate, nrow = 6) \%*\% 
        matrix(reads.per.replicate, ncol = 6)
    simple.precision.weights <- lower.tri(simple.precision.weights) * 
        simple.precision.weights
    simple.precision.clonality <- sum(simple.precision.weights * 
        clonality.matrix)/sum(simple.precision.weights)
    estimators.rownums <- diag(c(1:num.replicates)) \%*\% matrix(1, 
        num.replicates, num.replicates)
    estimators.colnums <- t(estimators.rownums)
    lowtri.indx <- as.vector(lower.tri(clonality.matrix))
    estimators.vec <- as.vector(clonality.matrix)[lowtri.indx]
    rownums.vec <- as.vector(estimators.rownums)[lowtri.indx]
    colnums.vec <- as.vector(estimators.colnums)[lowtri.indx]
    estimator.table.raw <- data.frame(theta = estimators.vec, 
        rn = rownums.vec, cn = colnums.vec)
    etr2 <- merge(estimator.table.raw, estimator.table.raw, by = c())
    indx <- etr2$cn.x == etr2$rn.y
    estimator.table.pairs1 <- data.frame(cn1 = etr2$cn.x[indx], 
        rn1 = etr2$rn.x[indx], cn2 = etr2$cn.y[indx], rn2 = etr2$rn.y[indx], 
        shared = etr2$cn.x[indx], rn = etr2$rn.x[indx], cn = etr2$cn.y[indx])
    indx <- (etr2$rn.x == etr2$rn.y) & (etr2$cn.x < etr2$cn.y)
    estimator.table.pairs2 <- data.frame(cn1 = etr2$cn.x[indx], 
        rn1 = etr2$rn.x[indx], cn2 = etr2$cn.y[indx], rn2 = etr2$rn.y[indx], 
        shared = etr2$rn.x[indx], rn = etr2$cn.y[indx], cn = etr2$cn.x[indx])
    indx <- (etr2$cn.x == etr2$cn.y) & (etr2$rn.x < etr2$rn.y)
    estimator.table.pairs3 <- data.frame(cn1 = etr2$cn.x[indx], 
        rn1 = etr2$rn.x[indx], cn2 = etr2$cn.y[indx], rn2 = etr2$rn.y[indx], 
        shared = etr2$cn.x[indx], rn = etr2$rn.y[indx], cn = etr2$rn.x[indx])
    estimator.table.pairs <- rbind(estimator.table.pairs1, estimator.table.pairs2, 
        estimator.table.pairs3)
    estimator.table.pairs$covariance <- diag(clonality.matrix)[estimator.table.pairs$shared]
    row.col.to.indx <- function(r, c) {
        return(((r - 1) * (r - 2))/2 + c)
    }
    indx.to.row <- rep(0, num.replicates * num.replicates)
    indx.to.col <- rep(0, num.replicates * num.replicates)
    curr.indx <- 0
    for (r in 2:num.replicates) {
        for (c in 1:(r - 1)) {
            curr.indx <- row.col.to.indx(r, c)
            indx.to.row[curr.indx] <- r
            indx.to.col[curr.indx] <- c
        }
    }
    estimator.table.pairs$cov.r <- row.col.to.indx(estimator.table.pairs$rn1, 
        estimator.table.pairs$cn1)
    estimator.table.pairs$cov.c <- row.col.to.indx(estimator.table.pairs$rn2, 
        estimator.table.pairs$cn2)
    num.pairs <- num.replicates * (num.replicates - 1)/2
    cov.matrix.halfdone <- sparseMatrix(i = estimator.table.pairs$cov.r, 
        j = estimator.table.pairs$cov.c, x = estimator.table.pairs$covariance)
    cov.matrix <- as.matrix(cov.matrix.halfdone + t(cov.matrix.halfdone))
    estimator.vec.forcov <- rep(0, num.pairs)
    for (r in 2:num.replicates) {
        for (c in 1:(r - 1)) {
            curr.indx <- row.col.to.indx(r, c)
            cov.matrix[curr.indx, curr.indx] <- clonality.matrix[r, 
                r] + clonality.matrix[c, c]
            estimator.vec.forcov[curr.indx] <- clonality.matrix[r, 
                c]
        }
    }
    root.prec.matrix <- sqrtm(ginv(cov.matrix))
    numerator <- rep(1, num.pairs) \%*\% root.prec.matrix \%*\% estimator.vec.forcov
    denominator <- t(rep(1, num.pairs)) \%*\% root.prec.matrix \%*\% 
        rep(1, num.pairs)
    rao.blackwell.mvg.clonality <- abs(numerator/denominator)
    return(c(simple.precision.clonality = simple.precision.clonality, 
        rao.blackwell.mvg.clonality = rao.blackwell.mvg.clonality))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
