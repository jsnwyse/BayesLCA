\name{Zscore}
\alias{Zscore}
\alias{Zscore.internal}
\title{
Evaluating Class Membership of Binary Data
}
\description{
For a fitted model of class \code{\link{blca}}, and binary data \code{X}, the probability of class membership for each data point is provided.
}
\usage{
Zscore(X, fit = NULL, itemprob = NULL, classprob = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A binary data matrix. \code{X} must have the same number of columns as the data that \code{fit} was applied to.
}
  \item{fit}{
An object of class \code{\link{blca}}.
}
  \item{itemprob}{
A matrix of item probabilities, conditional on class membership.
}
  \item{classprob}{
A vector denoting class membership probability.
}
}
\details{
Calculation of the probability of class membership for a data point relies on two parameters, class membership and item probability. These may be supplied directly to Zscore, or alternatively, a \code{\link{blca}} object containing both parameters can be used instead. 
}
\value{
A matrix of equal rows to \code{X} and with G, the number of classes, columns, where each row is a score denoting the probability of class membership. Each row should therefore sum to 1.
}
\author{
Arthur White
}
\note{
\code{Zscore.internal} has the same functionality as \code{Zscore}, but is only intended for internal use.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
set.seed(1)
type1 <- c(0.8, 0.8, 0.05, 0.2)
type2 <- c(0.2, 0.2, 0.05, 0.8)
x<- rlca(250, rbind(type1,type2), c(0.5,0.5))

fit <- blca.em(x, 2)
fit$Z ## Unique data types
Zscore(x, fit=fit) ## Whole data set
Zscore(c(0, 1, 1, 0), fit=fit) ## Not in data set
Zscore(x, itemprob=rbind(type1,type2), classprob=c(0.5,0.5))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ blca }
