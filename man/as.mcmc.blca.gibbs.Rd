\name{as.mcmc.blca.gibbs}
\alias{as.mcmc.blca.gibbs}
\alias{blca2mcmc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Converts \code{blca.gibbs} Objects to type \code{mcmc}
}
\description{
Converts blca objects to mcmc objects.  This is only to be used with the Gibbs sampling method.
}
 \usage{
%-\method{as.mcmc}{blca.gibbs}(x)
\method{as.mcmc}{blca.gibbs}(x, \dots)
 blca2mcmc(x)
 }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class blca.gibbs. An error is returned if this is not the case.
}
\item{\dots}{Additional arguments to be passed to the \code{\link[coda]{mcmc}} function. }
}
\details{
Whenever a Gibbs sampler is employed, it is always a good idea to ensure that parameter samples are being obtained correctly - that burn-in has been achieved, and that appropriate mixing is taking place, for example.  \code{as.mcmc.blca.gibbs} converts an object of class \code{\link{blca}} to that of \code{\link[coda]{mcmc}} to avail of the diagnostic  checks available in other R packages, particularly those in the \bold{coda} package.
}
\value{
An \eqn{N \times G*(M+1)} matrix of class \code{\link[coda]{mcmc}}, where N is the number of data points, M the number of columns and G the number of classes. The first G columns (labelled ClassProb 1 , \ldots, ClassProb G) are class membership probability samples, the next G*M columns (labelled ItemProb 1 1 , ItemProb 1 2, \ldots, ItemProb G 1, \ldots, ItemProb G M) are item response probability samples.
}
\author{
Arthur White
}

\note{
This function replaces the function \code{mcmc2blca}, which appeared in the original version of the package, and which is retained as an internal function for backwards compatibility reasons.

}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{blca.gibbs}},  \code{\link[coda]{geweke.diag}}, \code{\link[coda]{raftery.diag}}
}
\examples{
data(Alzheimer)

\dontrun{ fit.gibbs <- blca.gibbs(Alzheimer, 2) }
\dontrun{ raftery.diag(as.mcmc(fit.gibbs)) }

\dontrun{ fit.gibbs <- blca.gibbs(Alzheimer, 2, K=100000, accept=0.1, burn.in=100) }
\dontrun{ plot(as.mcmc(fit.gibbs)) }

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ blca }
\keyword{ as.mcmc }% __ONLY ONE__ keyword per line
