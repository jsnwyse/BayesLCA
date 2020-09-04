\name{blca}
\alias{blca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bayesian Latent Class Analysis with one of several methods
}
\description{
Latent class analysis (LCA) attempts to find G hidden classes in binary data X. blca utilises one of: an EM algorithm, a variational Bayes approximation, Gibbs sampling or boot-strapping techniques to find maximum {\emph{ a posteriori}} (MAP), standard error and density estimates of the parameters.
}
\usage{
blca(X, G, ncat=NULL, method = c("em", "gibbs", "boot", "vb", "collapsed"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
The data matrix. This may take one of several forms, see \code{\link{data.blca}}.
}
\item{ncat}{
	The number of categories in each variable (column).
}
  \item{G}{
The number of classes to run lca for.
}
  \item{method}{
The method with which to perform lca on the data. Four methods are currently available, "em", "gibbs", "boot", "collapsed" or "vb". Defaults to "em", with a warning.
}
  \item{\dots}{
Additional arguments to be passed on, depending on the method. See additional help files for details.
}
}
\details{
The function calls to one of \code{\link{blca.em}}, \code{\link{blca.boot}}, \code{\link{blca.gibbs}}, \code{\link{blca.vb}}, depending on the method specified.
}
\value{

A list of class "blca" is returned. All methods return the following items:
\item{classprob}{The class probabilities.}
\item{itemprob}{The item probabilities, conditional on class membership.}
\item{Z}{Estimate of class membership for each unique datapoint.}
\item{prior}{A list containing the prior values specified for the model.}

See additional help files for  details.
}
\author{
Arthur White
}

\references{
Arthur White, Thomas Brendan Murphy (2014). BayesLCA: An R Package for Bayesian Latent Class Analysis." Journal of Statistical Software, 61(13), 1-28. URL: http://www.jstatsoft.org/v61/i13/.

White, A., Wyse, J. and Murphy, T. B. (2016). Bayesian variable selection for latent class analysis using a collapsed Gibbs sampler. Statistics and Computing, volume 26, 511-527.

}

\note{
Earlier versions of this function erroneously referred to posterior standard deviations as standard errors. This also extended to some of the variable names of the returned function, which are now returned with the corrected suffix \code{blca.em.sd} (for standard deviation). For backwards compatability reasons, the earlier suffix \code{.se} has been retained as a returned argument.
}

\seealso{
\code{\link{blca.em}}, \code{\link{blca.boot}}, \code{\link{blca.gibbs}}, \code{\link{blca.vb}}
}
\examples{
type1 <- c(0.8, 0.8, 0.2, 0.2)
type2 <- c(0.2, 0.2, 0.8, 0.8)
x<- rlca(1000, rbind(type1,type2), c(0.6,0.4))

set.seed(1)
fit <- blca(x, 2) ## EM algorithm used, warning returned
print(fit)	  ## No posterior standard deviations returned
summary(fit)

set.seed(1)
fit2 <- blca(x, 2, method="em", sd=TRUE) ##No warning - same fit
print(fit2) 				 ##Posterior standard deviations returned

set.seed(1)
##Variational Bayes approximation, with priors specified.
fit3 <- blca(x, 2, method="vb", delta=5, alpha=2, beta=1) 
print(fit3)	##Posterior standard deviations returned also.
par(mfrow=c(3,2))			  
plot(fit3, which=3:4)
par(mfrow=c(1,1))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ blca }
\keyword{ lca }% __ONLY ONE__ keyword per line
