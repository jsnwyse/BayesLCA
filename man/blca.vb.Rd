\name{blca.vb}
\alias{blca.vb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bayesian Latent Class Analysis  via a variational Bayes algorithm
}
\description{
Latent class analysis (LCA) attempts to find G hidden classes in binary data X. blca.vb uses a variational EM algorithm to find the distribution which best approximates the parameters' true distribution.
}
\usage{
blca.vb(X, G, alpha = 1, beta = 1, delta = 1, 
	start.vals = c("single", "across"), counts.n = NULL, 
        iter = 500, restarts = 1, verbose = TRUE, conv = 1e-06, 
	small = 1e-100)
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
Denotes how class membership is to be assigned during the initial step of the algorithm. Two character values may be chosen, "single", which randomly assigns data points exclusively to one class, or "across", which assigns class membership via  \code{\link{runif}}. Alternatively, class membership may be pre-specified, either as a vector of class membership, or as a matrix of probabilities. Defaults to "single".
}
\item{counts.n}{
If data patterns have already been counted, a data matrix consisting of each unique data pattern can be supplied to the function, in addition to a vector counts.n, which supplies the corresponding number of times each pattern occurs in the data. 
}
  \item{iter}{
The maximum number of iterations that the algorithm runs over.  Will stop earlier if the algorithm converges.
}
  \item{restarts}{
\code{restarts} determines how many times the algorithm is run with different starting values. Parameter estimates from the run which achieved the highest log-posterior are returned. If starting values are supplied, these are used for the first run, after which random starting points are used. Defaults to 1. 
}
  \item{verbose}{
Logical valued. If TRUE, the log-posterior from each run is printed.
}
  \item{conv}{
Convergence criteria, i.e., how small should the log-posterior increase become before the algorithm is deemed to have converged?  Set relative to the size of the data matrix.
}
  \item{small}{
To ensure numerical stability a small constant is added to certain parameter estimates. Defaults to 1e-100.
}
}
\details{
The variational Bayes method approximates the posterior using as a product of independent distributions. Parameters are then estimated for this approximate distribution using a variational EM algorithm. This method has a tendency to underestimate parameter's variance; as such the standard error and density estimates should be interpreted with caution. 

While it is worth starting the algorithm from multiple starting points, variational algorithms have less of a tendency to cpnverge at saddle point or sub-optimal local maxima.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list of class "blca.vb" is returned, containing:
\item{call}{The initial call passed to the function.}
\item{itemprob}{The item probabilities, conditional on class membership.}
\item{classprob}{The class probabilities.}
\item{itemprob.sd}{Posterior standard deviation estimates of the item probabilities.}
\item{classprob.sd}{Posterior standard deviation estimates of the class probabilities.}
\item{parameters}{A list containing posterior parameter values for item and class probabilities, which are assumed to follow beta and Dirichlet distributions respectively.}
\item{Z}{Estimate of class membership for each unique datapoint.}
\item{LB}{The lower bound estimate of the log-posterior of the estimated model.}
\item{lbstore }{The value of the lower bound estimate for each iteration.}
\item{iter}{The number of iterations required before convergence.}
\item{eps}{The amount that the lower bound increased at the final iteration  of the algorithm's run. }
\item{counts}{ The number of times each unique datapoint point occured.}
\item{prior}{A list containing the prior values specified for the model.}
}
\references{
Ormerod J, Wand M (2010). ``Explaining Variational Approximations.'' The American Statistician, 64(2), 140-153. 
}
\author{
Arthur White
}
\note{
Variational Bayes approximations, are known to often underestimate the standard errors of the parameters under investigation, so caution is advised when checking their values.

Earlier versions of this function erroneously referred to posterior standard deviations as standard errors. This also extended to arguments supplied to and returned by the function, some of which are now returned with the corrected corrected suffix \code{blca.em.sd} (for standard deviation). For backwards compatability reasons, the earlier suffix \code{.se} has been retained as a returned argument.
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{blca.em}}, \code{\link{blca.gibbs}}
}
\examples{
type1 <- c(0.8, 0.8, 0.2, 0.2)
type2 <- c(0.2, 0.2, 0.8, 0.8)
x<- rlca(1000, rbind(type1,type2), c(0.6,0.4))

fit <- blca.vb(x, 2)
print(fit)
summary(fit)
par(mfrow=c(3,2))
plot(fit)
par(mfrow=c(1,1))

data(Alzheimer)
sj <- blca.vb(Alzheimer, 10, delta=1/10)
sj$classprob    ##Empty Groups
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ variational }
\keyword{ blca }% __ONLY ONE__ keyword per line
