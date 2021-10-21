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
blca.boot(X, G, formula=NULL, ncat=NULL, alpha=1, beta=1, delta=1, 
          start.vals= c("single","across"), counts.n=NULL, 
          fit=NULL, iter=2000,  B=100, verbose=TRUE, verbose.update=100,
          conv=1e-06, small=1e-10, MAP=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Data matrix or dataframe. This may take one of several forms, see \code{\link{data.blca}}.
}
  \item{G}{
The number of classes to run lca for.
}
  \item{formula}{
Optional formula declaring variables to be used from \code{X} e.g. ``\code{~ x1 + x2}''.
}
\item{ncat}{
	The number of categories in each variable (column).
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

  \item{verbose}{
Logical valued. If TRUE, the current number of completed bootstrap samples is printed at regular intervals.
}
  \item{verbose.update}{
If \code{verbose=TRUE}, \code{verbose.update} determines the periodicity with which updates are printed.
}

\item{conv}{Convergence threshold for runs of the EM algorithm within bootstrap samples.}

  \item{small}{
To ensure numerical stability a small constant is added to certain parameter estimates. Defaults to 1e-100.
}
\item{MAP}{
Logical indicating whether priors should be included in the analysis or not.
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
\item{method}{Method used to fit LCA to this dataset.}
\item{call}{The initial call passed to the function.}
\item{args}{List of arguments passed to the function.}
\item{G}{Number of classes.}
\item{classprob}{The class probabilities.}
\item{classprob.sd}{Posterior standard deviation or standard error.}
\item{itemprob}{The item probabilities, conditional on class membership. Either a matrix or a list dependending on whether there are only dichotomous variables or a mix of di and polychotomous variables.}
\item{itemprob.sd}{Posterior standard deviation or standard error.}
\item{itemprob.tidy}{Another representation of itemprob for use by plotting functions.}
\item{Z}{Estimate of class membership for each unique datapoint.}
\item{logpost}{The log-posterior evaluated at the estimated parameters i.e. the parameters from averaging over bootstrap replications.}
\item{loglik}{The maximum log-likelihood of the estimated model (only returned when \code{MAP=FALSE}). }
\item{AIC}{Akaike's Information Criterion for the estimated model.}
\item{BIC}{The Bayesian Information Criterion for the estimated model.}
\item{samples}{A list containing the parameter estimates for each bootstrap sample.}
\item{B}{Number of bootstrap samples.}
\item{prior}{A list containing the prior values specified for the model.}
\item{classprob.initial, itemprob.initial}{Initial parameter values for classprob and itemprob, used to run over each bootstrapped sample.}
\item{model.indicator}{Binary vector indicating the variables which were included for clustering.}
\item{MAP}{Logical indicating whether the MAP was targeted i.e. whether the maximum of the log posterior or log likelihood was used.}
\item{ncat}{The ncat vector from the call.}
}

\references{
   Wasserman, L, 22nd May 2007, \emph{All of Nonparametric Statistics}, Springer-Verlag.
}
\author{
Arthur White, Jason Wyse
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

fit <- blca.em(x, 2, sd=FALSE)
fit.boot <- blca.boot(x, 2, fit=fit)
fit.boot
plot(fit.boot, which=1:4)
}
\keyword{ bootstrap }
\keyword{ blca }% __ONLY ONE__ keyword per line
