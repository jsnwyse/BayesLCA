\name{blca.collapsed}
\alias{blca.collapsed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bayesian Latent Class Analysis variable and group selection via collapsed sampling
}
\description{
Latent class analysis (LCA) attempts to find G hidden classes in binary data X. blca.collapsed performs sampling from a collapsed posterior distribution to select the  number of hidden classes G and the most relevant clustering variables.
}
\usage{
blca.collapsed(X , G, alpha = 1, beta = 1, delta = 1, num.categories = NULL, 
iter = 6000, burn.in = 1000, thin = 1, fixed.G = FALSE, 
just.gibbs.updates = FALSE, n.gibbs = nrow(X), prior.G = 1, 
G.max =  30, variable.selection = FALSE, prob.inc = .5, 
hprior.model = FALSE, relabel = TRUE, silent = TRUE  )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{X}{
The data matrix. 
}
  \item{G}{
The number of classes to run lca for in the case of \code{fixed.G  = TRUE}, otherwise this gives the initial number of groups for the sampler.
}

  \item{alpha, beta}{
The prior values for the data conditional on group membership. Defaults to 1, i.e, a uniform prior, for each value.
}
  \item{delta}{
Prior values for the mixture components in model.  Defaults to 1, i.e., a uniform prior.  May be single or vector valued (of length G).
}

\item{num.categories}{
A vector giving the number of categories for each variable. If not specified, a test for binary data is carried out in which case 2 categories will be assumed.
}

  \item{iter}{
The number of iterations to run the collapsed gibbs sampler for \bold{after} burn-in.
}

  \item{burn.in}{
Number of iterations to run the Gibbs sampler for before beginning to store values.
}

  \item{thin}{
The thinning rate for samples from the distribution, in order to achieve good mixing. Should take a value greater  >0 and <=1. Defaults to 1.
}
 \item{fixed.G}{
Logical, if G is allowed to vary over the run of the algorithm.
}
  \item{just.gibbs.updates}{
Logical, indicating if just gibbs sampling alone should be used to update labels, or to allow other moves which can reallocated a number of items at a time.
}
  \item{n.gibbs}{
For large data, only update a random sample of \code{n.gibbs} item labels.
}
\item{prior.G}{
Uniform (value 0) or truncated Poisson(rate=1) prior on the number of latent classes G.
}
\item{G.max}{
The maximum number of groups to search for (default 30).
}
\item{variable.selection}{
Logical, inidicating if the sampler should also search for the most relevant clustering variables. 
}
\item{prob.inc}{
Prior probability of including a variable for clustering when inclusion is modelled as a Bernoulli distriubtion a priori.
}
\item{hprior.model}{
Logical, indicating if the prior probability of variable inclusion should be modelled using a hyperprior and sampled at each iteration.
}
\item{relabel}{
Logical, indicating whether a mechanism to correct for label-switching be used or not. Defaults to TRUE.
}
\item{write.to.file}{
Logical, indicating whether the results should be written to file in the current working directory (useful for large datasets).
} 
\item{filenames}{
If \code{write.to.file} is \code{TRUE} a list of file names.
} 
\item{silent}{
Logical, indicating whether to print a summary to screen when finished.
}

}
\details{

}
\value{
A list of class "blca.collapsed" is returned, containing:
\item{call}{The initial call passed to the function.}
\item{Z}{List giving the estimated class membership for each unique datapoint for each value of G visited by the sampler. Entries in the list are indexed 1,2,..., each a matrix with the corresponding number of components for entry 1 in \code{G.Z[1]}. Only computed if \code{relabel=TRUE}}
\item{G.Z}{Vector indicating the values of G corresponding to entries of the list \code{Z}}
\item{samples}{A list containing the output of the sampler at each iteration.}
\item{prior}{A list containing the prior values specified for the model.}
\item{iter}{The number of iterations the collapsed sampler was run for.}
\item{thin}{The sub-sampling gap used when storing samples from the distribution.}
\item{burn.in}{The number of iterations the collapsed sampler was run before beginning to store values.}
\item{fixed.G}{Logical, indicating whether G was static throughout the sampler run}
\item{just.gibbs.updates}{Logical, indicating whether only Gibbs updates were used to update data item labels.}
\item{n.gibbs}{Number of label updates applied during each sweep of the algorithm.}
\item{prior.G}{Indicates the type of prior taken on the number of latent classes G.}
\item{G.max}{The maximum number of  groups allowed to be  visited by the algorithm if \code{fixed.G} was \code{FALSE}}
\item{variable.selection}{Logical, indicating whether the sampler searched for relevant clustering variables.}
\item{hprior.model}{Logical, indicating whether \code{prob.inc} was sampled using a hyperprior.}
\item{relabel}{Logical, indicating whether a mechanism to correct for label-switching was used.}
\item{labelstore}{The stored labels during the sampling run. If relabel=TRUE, these have already been permuted using the label correcting algorithm.}
\item{X}{The matrix \code{X} that was passed when calling the algorithm}
}
\references{
White, A., Wyse, J. and Murphy, T. B. (2016). Bayesian variable selection for latent class analysis using a collapsed Gibbs sampler. Statistics and Computing, volume 26, 511-527.
}
\author{
Jason Wyse
}

\seealso{
\code{\link{blca}}, \code{\link{blca.gibbs}} 
}
\examples{
#something in here...
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ blca }
\keyword{ collapsed }% __ONLY ONE__ keyword per line
