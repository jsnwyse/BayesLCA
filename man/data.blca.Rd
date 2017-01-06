\name{data.blca}
\alias{data.blca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Conveniently Format Data for Bayesian Latent Class 
}
\description{
Conveniently format data for use with \code{\link{blca}}.
}
\usage{
data.blca(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A data matrix intended for latent class analysis. See details.
}
}
\details{
The data may take of one of two forms, either as a binary matrix, or as a matrix consisting of unique binary rows, with a column of counts. In either case \code{\link{data.blca}} will convert X into a list, with binary matrix and count vector explicitly identified.
}
\value{
A list of class data.blca, containing

  \item{counts.n }{A vector of counts of each unique data entry.}
  \item{data }{A matrix consisting of each unique data entry.}

}

\author{
Arthur White
}
\note{
This function is used internally by \code{\link{blca}}, so its use is not necessary, though it will speed up computation time to supply the model with data of this form if repeated use of a function is required.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{blca}}
}
\examples{
type1 <- c(0.8, 0.8, 0.2, 0.2)
type2 <- c(0.2, 0.2, 0.8, 0.8)
x<- rlca(1000, rbind(type1,type2), c(0.6,0.4)) ##Only 16 unique observations possible

data.blca(x)
}
\keyword{ blca }
\keyword{ data.blca }% __ONLY ONE__ keyword per line
