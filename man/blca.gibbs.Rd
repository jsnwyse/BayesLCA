\name{blca.gibbs}
\alias{blca.gibbs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bayesian Latent Class Analysis via Gibbs Sampling
}
\description{
Latent class analysis (LCA) attempts to find G hidden classes in binary data X. blca.gibbs performs Gibbs sampling to sample from the parameters' true distribution.
}
\usage{
blca.gibbs(X, G, alpha = 1, beta = 1, delta = 1, 
	   start.vals = c("prior", "single", "across"), 
	   counts.n = NULL, iter = 5000, thin = 1, 
	   accept=thin, burn.in = 100, relabel = TRUE, 
           verbose = TRUE, verbose.update = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{X}{
The data matrix. This may take one of several forms, see \code{\link{data.blca}}.
}
  \item{G}{
The number of classes to run lca for.
}

  \item{alpha, beta}{
The prior values for the data conditional on group membership. These may take several forms: a single value, recycled across all groups and columns, a vector of length G or M (the number of columns in the data), or finally, a G x M matrix specifying each prior value separately. Defaults to 1, i.e, a uniform prior, for each value.
}
  \item{delta}{
Prior values for the mixture components in model.  Defaults to 1, i.e., a uniform prior.  May be single or vector valued (of length G).
}

 \item{start.vals}{
Denotes how class membership is to be assigned during the initial step of the algorithm. One of three character values may be chosen: "prior", which samples parameter values from prior distributions, "single", which randomly assigns data points exclusively to one class, or "across", which assigns class membership via  \code{\link{runif}}. Alternatively, class membership may be pre-specified, either as a vector of class membership, or as a matrix of probabilities. Defaults to "single".
}
\item{counts.n}{
If data patterns have already been counted, a data matrix consisting of each unique data pattern can be supplied to the function, in addition to a vector counts.n, which supplies the corresponding number of times each pattern occurs in the data. 
}
  \item{iter}{
The number of iterations to run the gibbs sampler for \bold{after} burn-in.
}
  \item{thin}{
The thinning rate for samples from the distribution, in order to achieve good mixing. Should take a value greater  >0 and <=1. Defaults to 1.
}
  \item{accept}{
Similarly to \code{accept}, specifies the thinning rate for samples from the distribution, in order to achieve good mixing, however, its use is discouraged. Should always agree with \code{sd}. Retained for backwards compatability reasons. See `Note'.
}
  \item{burn.in}{
Number of iterations to run the Gibbs sampler for before beginning to store values.
}
 \item{relabel}{
Logical, indicating whether a mechanism to prevent label-switching be used or not. Defaults to TRUE.
}
  \item{verbose}{
Logical valued. If TRUE, the current number of completed samples is printed at regular intervals.
}
  \item{verbose.update}{
If \code{verbose=TRUE}, \code{verbose.update} determines the periodicity with which updates are printed.
}
}
\details{
The library coda provide extensive tools to diagnose and visualise MCMC chains. The generic function \code{\link{as.mcmc.blca.gibbs}}, makes \code{blca.gibbs} objects compatible with functions such as \code{\link[coda]{summary.mcmc}} and \code{\link[coda]{raftery.diag}}.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list of class "blca.gibbs" is returned, containing:
\item{call}{The initial call passed to the function.}
\item{classprob}{The class probabilities.}
\item{itemprob}{The item probabilities, conditional on class membership.}
\item{classprob.sd}{Posterior standard deviation estimates of the class probabilities.}
\item{itemprob.sd}{Posterior standard deviation estimates of the item probabilities.}
\item{logpost}{The log-posterior of the estimated model.}
\item{Z}{Estimate of class membership for each unique datapoint.}
\item{samples}{A list containing Gibbs samples of the item and class probabilities and log-posterior.}
\item{DIC}{The Deviance Information Criterion for the estimated model.}
\item{BICM}{The Bayesian Information Criterion (Monte Carlo) for the estimated model.}
\item{AICM}{Akaike's Information Criterion (Monte Carlo) for the estimated model.}
\item{counts}{ The number of times each unique datapoint point occured.}
\item{prior}{A list containing the prior values specified for the model.}
\item{thin}{The acceptance rate for samples from the distribution.}
\item{burn.in}{The number of iterations the gibbs sampler was run before beginning to store values.}
\item{relabel}{Logical, indicating whether a mechanism to prevent label-switching was used.}
\item{labelstore}{The stored labels during the sampling run. If relabel=TRUE, these show how labels were permuted in an attempt to avoid label-switching in the model.}
}
\references{
Spiegelhalter DJ, Best NG, Carlin BP, Linde Avd (2002). ``Bayesian Measures of Model Complexity and Fit.'' Journal of the Royal Statistical Society. Series B (Statistical Methodology), 
64(4), pp. 583-639. ISSN 13697412. URL http://www.jstor.org/stable/3088806.

Raftery AE, Newton MA, Satagopan JM, Krivitsky PN (2007). ``Estimating the integrated 
likelihood via posterior simulation using the harmonic mean identity.'' In Bayesian Statistics, 
pp. 1-45. 

}
\author{
Arthur White
}
\note{
Earlier versions of this function erroneously referred to posterior standard deviations as standard errors. This also extended to arguments supplied to and returned by the function, some of which are now returned with the corrected corrected suffix \code{blca.em.sd} (for standard deviation). For backwards compatability reasons, the earlier suffix \code{.se} has been retained as a returned argument.
The argument \code{thin} replaces \code{accept}, which appeared in the earliest version of the package. This is to maintain consistency with other packages, such as \bold{rjags}.
}
\seealso{
\code{\link{blca}}, \code{\link{as.mcmc.blca.gibbs}}, \code{\link[coda]{raftery.diag}}
}
\examples{
## Generate a 4-dim. sample with 2 latent classes of 500 data points each.
## The probabilities for the 2 classes are given by type1 and type2.

type1 <- c(0.8, 0.8, 0.2, 0.2)
type2 <- c(0.2, 0.2, 0.8, 0.8)
x<- rlca(1000, rbind(type1,type2), c(0.6,0.4))

\dontrun{fit.gibbs<-blca.gibbs(x,2, iter=1000, burn.in=10)}
\dontrun{summary(fit.gibbs)}
\dontrun{plot(fit.gibbs)}
\dontrun{raftery.diag(as.mcmc(fit.gibbs))}


\dontrun{fit.gibbs<-blca.gibbs(x,2, iter=10000, burn.in=100, thin=0.5) }
\dontrun{plot(fit.gibbs, which=4)}
\dontrun{raftery.diag(as.mcmc(fit.gibbs))}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ blca }
\keyword{ gibbs }% __ONLY ONE__ keyword per line
