\name{blca.em}
\alias{blca.em}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bayesian Latent Class Analysis  via an EM Algorithm
}
\description{
Latent class analysis (LCA) attempts to find G hidden classes in binary data X. blca.em utilises an expectation-maximisation algorithm to find maximum \emph{a posteriori} (map) estimates of the parameters.
}
\usage{
blca.em(X, G, ncat=NULL,  alpha=1, beta=1, delta=1, 
	start.vals=c("single", "across"), counts.n=NULL,
	model.indicator=NULL, iter=1000, restarts=5, verbose=TRUE, 
	sd = FALSE, se=sd, conv=1e-06, small=1e-100, MAP=TRUE, pars.init=NULL, for.boot=FALSE )
}


%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
The data matrix. This may take one of several forms, see \code{\link{data.blca}}.
}

  \item{G}{
The number of classes to run lca for.
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
Denotes how class membership is to be assigned during the initial step of the algorithm. Two character values may be chosen, "single", which randomly assigns data points exclusively to one class, or "across", which assigns class membership via  \code{\link[stats]{runif}}. Alternatively, class membership may be pre-specified, either as a vector of class membership, or as a matrix of probabilities. Defaults to "single".
}

\item{counts.n}{
If data patterns have already been counted, a data matrix consisting of each unique data pattern can be supplied to the function, in addition to a vector counts.n, which supplies the corresponding number of times each pattern occurs in the data. 
}

\item{model.indicator}{
Binary vector of length the number of columns of X giving the variables to include when fitting the LCA model (a 1 means include).
}
  \item{iter}{
The maximum number of iterations that the algorithm runs over.  Will stop early if the algorithm is deemed to converge.
}
  \item{restarts}{
\code{restarts} determines how many times the algorithm is run with different starting values. Parameter estimates from the run which achieved the highest log-posterior are returned. If starting values are supplied, these are used for the first run, after which random starting points are used. Defaults to 5. 
}
  \item{verbose}{
Logical valued. If \code{TRUE}, the log-posterior from each run is printed.
}
  \item{sd}{
Specifies whether posterior standard deviation estimates should also be returned. If TRUE, calls to {\link{blca.em.sd}}.
}
  \item{se}{
Similarly to \code{sd}, specifies whether posterior standard deviation estimates should also be returned, however, its use is discouraged. Should always agree with \code{sd}. Retained for backwards compatability reasons. See `Note'.
}
  \item{conv}{
Convergence criteria, i.e., how small should the log-posterior increase become before the algorithm is deemed to have converged?  Set relative to the size of the data matrix.
}
  \item{small}{
To ensure numerical stability a small constant is added to certain parameter estimates. Defaults to 1e-100.
}
\item{MAP}{
Logical determining whether the maximum a posteriori estimate should be computed. If \code{TRUE} then the priors are included in the M-step update.
}
\item{pars.init}{Parameter for use by the fuction \code{\link{blca.boot}} to provide initialization to the EM algorithm.}
\item{for.boot}{
Logical argument used by the function \code{\link{blca.boot}}. Should be set to \code{TRUE} if the result of the EM call will be used as input to the function \code{\link{blca.boot}}.
}
}
\details{
Regardless of the form of the data supplied to \code{\link{blca.em}}, it is internally converted to be of the form \code{\link{data.blca}}.  In particular, this should be noted when supplying starting values: the object must be of either the same length or have the same number of rows as the number of \bold{unique} observations in the dataset, as opposed to the total number.

Posterior standard deviations and convergence checks are calculated using \code{\link{blca.em.sd}}. 
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list of class "blca.em" is returned, containing:
\item{call}{The initial call passed to the function.}
\item{itemprob}{The item probabilities, conditional on class membership. Either a matrix or a list dependending on whether there are only dichotomous variables or a mix of di and polychotomous variables.}
\item{itemprob.tidy}{Another representation of itemprob for use by plotting functions.}
\item{classprob}{The class probabilities.}
\item{itemprob.sd}{If returned, standard error estimates of the item probabilities. Either a matrix or a list dependending on whether there are only dichotomous variables or a mix of di and polychotomous variables.}
\item{classprob.sd}{If returned, standard error estimates of the class probabilities.}
\item{Z}{Estimate of class membership for each unique datapoint.}
\item{logpost}{The log-posterior of the estimated model.}
\item{BIC}{The Bayesian Information Criterion for the estimated model.}
\item{AIC}{Akaike's Information Criterion for the estimated model.}
\item{iter}{The number of iterations required before convergence.}
\item{poststore }{The value of the log-posterior for each iteration.}
\item{converged}{Logical indicating whether the algorithm converged.}
\item{eps}{The value for which the algorithm was deemed to have converged at.}
\item{conv}{The value of \code{conv} from the call.}
\item{lpstarts}{ The log-posterior achieved after each of the multiple starts of the algorithm. }
\item{convergence}{If posterior standard deviations are calculated, then the Hessian of the model is also checked to determine whether the algorithm has converged to at least a local maximum. The convergence status is calculated by an integer value: 1 denotes acceptable convergence, 2 denotes that it converged at a saddle point, 3 that the algorithm ended before it had satisfactorily converged  and 4 denotes that at least one parameter value converged at a boundary value (i.e., a 1 or 0). 0 denotes that the algorithm converged satisfactorily but that the Hessian has not been checked.}
\item{prior}{A list containing the prior values specified for the model.}
\item{sd}{A logical value indicating whether standard error estimates were returned. }
\item{model.indicator}{Binary vector indicating the variables which were included for clustering.}
\item{MAP}{Logical indicating whether the MAP was computed i.e. whether the maximum of the log posterior or log likelihood was calculated.}
\item{ncat}{The ncat vector from the call.}
\item{pars.init}{Argument mainly for use by \code{blca.boot} to initialize the EM algorithm using the parameters found from a previous EM run.}
\item{for.boot}{Logical indicating whether the EM run is to be used as an initializer for a bootstrap run.}
\item{boot.init}{If \code{for.boot == TRUE} then this is a list giving a convenient arrangement of maximized parameters to initialize the bootstrap.}
}
\references{
%% ~put references to the literature/web site here ~
Dempster AP, Laird NM, Rubin DB (1977). ``Maximum Likelihood from Incomplete Data via 
the EM Algorithm.'' Journal of the Royal Statistical Society. Series B (Methodological), \bold{39(1)}, 
pp. 1--38. ISSN 00359246. doi:10.2307/2984875. URL http://dx.doi.org/10.2307/2984875.
}
\author{
Arthur White
}

\note{
Earlier versions of this function erroneously referred to posterior standard deviations as standard errors. This also extended to arguments supplied to and returned by the function, some of which are now returned with the corrected corrected suffix \code{blca.em.sd} (for standard deviation). For backwards compatability reasons, the earlier suffix \code{.se} has been retained as a returned argument.
}

\seealso{
\code{\link{blca}},\code{\link{blca.em.sd}}, \code{\link{blca.boot}}, \code{\link{blca.vb}}, \code{\link{blca.collapsed}}
}
\examples{
type1 <- c(0.8, 0.8, 0.2, 0.2)
type2 <- c(0.2, 0.2, 0.8, 0.8)
x <- rlca(1000, rbind(type1,type2), c(0.6,0.4))

fit <- blca.em(x, 2)
print(fit)
fit <- blca.em(x, 2, sd=TRUE) ##Returns posterior standard deviations
summary(fit)
plot(fit)

## Different starting values
fit <- blca.em(x, 2, start.vals="across")
xx <- data.blca(x)
fit <- blca.em(xx, 2, start.vals=sample(1:2, sum(xx$counts.n) , replace=TRUE))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ em }
\keyword{ blca }% __ONLY ONE__ keyword per line
