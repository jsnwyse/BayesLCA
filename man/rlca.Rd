\name{rlca}
\alias{rlca}

\title{
Randomly Generate Binary Data with Underlying Latent Classes
}
\description{
A function which randomly generates data with respect to some underlying latent class. Data may be generated either by specifying item and class probabilities or by utilising an object previously fitted to data.
}
\usage{
rlca(n, itemprob = 0.5, classprob = 1, fit = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
Number of data points to be generated.
}
  \item{itemprob}{
The item probabilities, conditional on class membership. Defaults to 0.5.
}
  \item{classprob}{
The class probabilities. Defaults to 1, i.e., a one class model.
}
  \item{fit}{
An object of class blca. If fit is supplied, data is generated using the class and item probabilities obtained. Defaults to NULL.
}
}


\author{
Arthur White
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{data.blca}}
}
\examples{
type1 <- c(0.8, 0.8, 0.2, 0.2)
type2 <- c(0.2, 0.2, 0.8, 0.8)
x<- rlca(1000, rbind(type1,type2), c(0.6,0.4))

fit <- blca.em(x, 2)

x2<- rlca(1000, fit=fit)
fit2<- blca.em(x2,2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{blca }
\keyword{ random }% __ONLY ONE__ keyword per line
