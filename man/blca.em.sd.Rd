\name{blca.em.sd}
\alias{blca.em.sd}
\alias{blca.em.se}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Posterior Standard Deviation Estimates for Bayesian Latent Class Analysis  via an EM Algorithm
}
\description{
Returns posterior standard deviation estimates for point estimates returned by \code{\link{blca.em}}. These are obtained via asymptotic estmation of the Observed Information matrix. The Hessian of the log-posterior is also checked to determine whether point estimates occur at at least a local maximum.
}
\usage{
blca.em.sd(fit, x, counts.n = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
An object of class "blca.em". 
}
  \item{x}{
A binary matrix. An object of class \code{\link{data.blca}} may also be supplied. In this case the argument \code{counts.n} is ignored.
}
  \item{counts.n}{
A vector which supplies the corresponding number of times each pattern in X occurs in the data.
}
}
\details{
This function is primarily intended for use in conjunction with \code{\link{blca.em}}, and may be called directly by that function by setting \code{se=TRUE}. However it can in fact be used with any blca object. 
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list containing:
\item{itemprob }{Posterior standard deviation estimates of the item probabilities.}
\item{classprob }{Posterior standard deviation estimates of the class probabilities.}
\item{convergence}{An integer value denoting whether point estimates occur at at least a local maximum. 1 denotes acceptable convergence, 2 denotes that it converged at a saddle point, 3 that the algorithm ended before it converged and 4 denotes that at least one parameter value converged at a boundary value.}
}
\author{
Arthur White
}
\note{
The posterior standard deviation estimates are derived asymptotically, i.e., by inverting the information matrix of the parameters.  These values are known to be unreliable in cases where parameters estimates are close to  1 or 0, so caution is advised when checking their values.  Bootstrapping methods may provide better estimates.

Computationally, the method becomes becomes unstable for values close to 1 or 0.  If the distance of any of the supplied parameter values from 0 or 1 is  <1e-5, then posterior standard deviation estimates for these values are returned as 0.

Earlier versions of this function erroneously referred to posterior standard deviations as standard errors.  This also extended to the function name, which has now been corrected to \code{blca.em.sd} (for standard deviation).  For backwards compatability reasons, the earlier function \code{blca.em.se} has been retained as an in internal function.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{blca.em}}, \code{\link{blca.boot}}
}
\examples{
type1 <- c(0.8, 0.8, 0.2, 0.2)
type2 <- c(0.2, 0.2, 0.8, 0.8)
x<- rlca(1000, rbind(type1,type2), c(0.6,0.4))
dat<- data.blca(x)

set.seed(1)
fit1 <- blca.em(dat, 2, se=TRUE)
fit1$itemprob.sd
fit1$classprob.sd

set.seed(1)
fit2<- blca.em(dat, 2, se=FALSE)
fit2.sd<- blca.em.sd(fit2, dat)
fit2.sd$itemprob
fit2.sd$classprob
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ posterior }
\keyword{ standard }
\keyword{ deviation }
\keyword{ blca }% __ONLY ONE__ keyword per line
