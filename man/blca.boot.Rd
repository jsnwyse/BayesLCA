\name{blca.boot}
\alias{blca.boot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bayesian Latent Class Analysis via an EM Algorithm and Using Empirical Bootstrapping
}
\description{
Latent class analysis (LCA) attempts to find G hidden classes in binary data X.  blca.boot repeatedly samples from X with replacement then utilises an EM algorithm to find maximum posterior (MAP) and standard error estimates of the parameters.
}
\usage{
blca.boot(X, G, alpha = 1, beta = 1, delta = rep(1, G), 
	  start.vals = c("single", "across"), counts.n = NULL, 
	  fit = NULL, iter = 50, B = 100, relabel = FALSE, 
          verbose = TRUE, verbose.update = 10, small = 1e-100)
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
The prior values for the data conditional on group membership. These may take several forms: a single value, recycled across all groups and columns, a vector of length G or M (the number of columns in the data), or finally, a \eqn{G \times M} matrix specifying each prior value separately. Defaults to 1, i.e, a uniform prior, for each value.
}
  \item{delta}{
Prior values for the mixture components in model.  Defaults to 1, i.e., a uniform prior.  May be single or vector valued (of length G).
}
 \item{start.vals}{
Denotes how class membership is to be assigned during the initial step of the algorithm. Two character values may be chosen, "single", which randomly assigns data points exclusively to one class, or "across", which assigns class membership via  \code{\link{runif}}. Alternatively, class membership may be pre-specified, either as a vector of class membership, or as a matrix of probabilities. Defaults to "single".
}
\item{counts.n}{
If data patterns have already been counted, a data matrix consisting of each unique data pattern can be supplied to the function, in addition to a vector counts.n, which supplies the corresponding number of times each pattern occurs in the data. 
}
 \item{fit}{
Previously fitted models may be supplied in order to approximate standard error and unbiased point estimates. fit should be an object of class "blca.em". Defaults to NULL if no object is supplied.
}
\item{iter}{
The maximum number of iterations that the algorithm runs over, for each bootstrapped sample.  Will stop earlier if the algorithm converges.
}
  \item{B}{
The number of bootstrap samples to run. Defaults to 100.
}
  \item{relabel}{
Logical valued. As the data is recursively sampled, it is possible that label-switching may occur with respect to parameter estimates.  If TRUE, parameter estimates are checked at each iteration, and relabeled if necessary.  Defaults to FALSE.
}
  \item{verbose}{
Logical valued. If TRUE, the current number of completed bootstrap samples is printed at regular intervals.
}
  \item{verbose.update}{
If \code{verbose=TRUE}, \code{verbose.update} determines the periodicity with which updates are printed.
}
  \item{small}{
To ensure numerical stability a small constant is added to certain parameter estimates. Defaults to 1e-100.
}
}
\details{
Bootstrapping methods can be used to estimate properties of a distribution's parameters, such as the standard error estimates, by constructing multiple resamples of an observed dataset, obtained by sampling with replacement from said dataset.  The multiple parameter estimates obtained from these resamples may then be analysed.  This method is implemented in blca.boot by first running blca.em over the full data set and then using the returned values of the item and class probabilities as the initial values when running the algorithm for each bootstrapped sample.  Alternatively, initial parameter estimates may be specified using the fit argument.

Note that if a previously fitted model is supplied, then the prior values with which the model was fitted will be used for the sampling run, regardless of the values supplied to the prior arguments.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list of class "blca.boot" is returned, containing:
\item{call}{The initial call passed to the function.}
\item{itemprob}{The item probabilities, conditional on class membership.}
\item{classprob}{The class probabilities.}
\item{Z}{Estimate of class membership for each unique datapoint.}
\item{itemprob.sd}{Posterior standard deviation estimates of the item probabilities.}
\item{classprob.sd}{Posterior standard deviation estimates of the class probabilities.}
\item{classprob.initial, itemprob.initial}{Initial parameter values for classprob and itemprob, used to run over each bootstrapped sample.}
\item{samples}{A list containing the parameter estimates for each bootstrapped sample.}
\item{logpost}{The log-posterior of the estimated model.}
\item{BIC}{The Bayesian Information Criterion for the estimated model.}
\item{AIC}{Akaike's Information Criterion for the estimated model.}
\item{label}{Logical value, indicating whether label switching has been checked for.}
\item{counts}{ The number of times each unique datapoint point occured.}
\item{prior}{A list containing the prior values specified for the model.}
}
\references{
   Wasserman, L, 22nd May 2007, \emph{All of Nonparametric Statistics}, Springer-Verlag.
}
\author{
Arthur White
}

\note{
Earlier versions of this function erroneously referred to posterior standard deviations as standard errors. This also extended to arguments supplied to and returned by the function, some of which are now returned with the corrected corrected suffix \code{blca.em.sd} (for standard deviation). For backwards compatability reasons, the earlier suffix \code{.se} has been retained as a returned argument.
}
\seealso{
\code{\link{blca.em}}, \code{\link{blca}}
}
\examples{
type1 <- c(0.8, 0.8, 0.2, 0.2)
type2 <- c(0.2, 0.2, 0.8, 0.8)
x <- rlca(1000, rbind(type1,type2), c(0.6,0.4))
fit.boot <- blca.boot(x, 2)
summary(fit.boot)

fit <- blca.em(x, 2, se=FALSE)
fit.boot <- blca.boot(x, 2, fit=fit)
fit.boot
plot(fit.boot, which=1:4)
}
\keyword{ bootstrap }
\keyword{ blca }% __ONLY ONE__ keyword per line
