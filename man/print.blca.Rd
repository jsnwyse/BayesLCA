\name{print.blca}
\alias{print.blca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bayesian Latent Class Analysis
}
\description{
Print a \code{\link{blca}} object.
}
\usage{
\method{print}{blca}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class \code{\link{blca}}.
}
  \item{\dots}{
Additional arguments to be passed onto lower-level functions at a later stage of development.
}
}
\details{
Prints parameter maximum \emph{a posteriori} (map) and standard deviation estimates. The latter are sometimes unavailable for \code{\link{blca.em}} objects.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
The blca object itself.
}

\author{
Arthur White
}



\examples{

data(Alzheimer)
fit1<- blca(Alzheimer, 2, method="em")
class(fit1)
print(fit1)
fit2<- blca(Alzheimer, 2, method="em", sd=TRUE)
print(fit2) ## Standard Errors also printed

fit3<- blca(Alzheimer, 2, method="vb")
print(fit3) ## Standard Errors as standard
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ blca}
\keyword{ print }% __ONLY ONE__ keyword per line
