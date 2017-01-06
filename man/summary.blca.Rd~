\name{summary.blca}
\alias{summary.blca}
\alias{summary.blca.boot}
\alias{summary.blca.em}
\alias{summary.blca.gibbs}
\alias{summary.blca.vb}
\alias{print.summary.blca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bayesian Latent Class Analysis
}
\description{
Summary method for class "blca".
}
\usage{
\method{summary}{blca}(object, ...)
%\method{print}{summary.blca}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Object of class \code{\link{blca}}.
}
  \item{\dots}{
Additional arguments to be passed onto lower-level functions at a later stage of development.
}
}

\value{
A brief summary consisting of two parts: the prior values specified to the model,  and model diagnostics specific to the inference method used, such as information about the log-posterior (or lower bound in the case of \code{blca.vb}), as well the number of iterations the algorithm ran for, etc..
}

\author{
Arthur White
}

\examples{
data(Alzheimer)
summary(blca.em(Alzheimer, 2))
summary(blca.vb(Alzheimer, 2, alpha=2, beta=2, delta=0.5))

\dontrun{(fit.gibbs)<- blca.gibbs(Alzheimer, 2, delta=2)}
\dontrun{summary(fit.gibbs)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ summary}
\keyword{ blca }% __ONLY ONE__ keyword per line
