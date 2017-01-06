\name{blca.collapsed.post.hoc.estimates}
\alias{blca.collapsed.post.hoc.estimates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bayesian Latent Class Analysis variable and group selection via collapsed sampling
}
\description{
Latent class analysis (LCA) attempts to find G hidden classes in binary data X. Use post-hoc methods to extract parameter estimates from a blca.collapsed object.
}
\usage{
blca.collapsed.post.hoc.estimates( x, variables, G )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{x}{
	An object of class \code{blca.collapsed}.
}
  \item{variables}{
	A vector giving the indexes of the variables to be used for the post-hoc analysis (these have to be fixed). 
}

  \item{G}{
	The number of latent classes to be used when carrying out the post-hoc analysis (this has to be fixed).
}
  

}
\details{

}
\value{
A list of class is returned, containing:
\item{call}{The initial call passed to the function.}
\item{classprob}{Vector giving the post-hoc estimates of the class probabilities.}
\item{sd.classprob}{Vector giving associated standard deviations of the class probabilities.}
\item{variables}{A vector variables for which the item probabilities and associated standard deviations were computed.}
\item{itemprob}{List with entry \code{i} the estimated item probabilites for \code{variables[i]}. }
\item{sd.itemprob}{List with entry \code{i} the estimated standard deviations of the item probabilities for \code{variables[i]}. }

}
\references{
White, A., Wyse, J. and Murphy, T. B. (2016). Bayesian variable selection for latent class analysis using a collapsed Gibbs sampler. Statistics and Computing, volume 26, 511-527.
}
\author{
Jason Wyse
}

\seealso{
\code{\link{blca}}, \code{\link{blca.collapsed}} 
}
\examples{
#something in here...
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ blca }
\keyword{ collapsed }% __ONLY ONE__ keyword per line
