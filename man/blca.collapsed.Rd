\name{blca.collapsed}
\alias{blca.collapsed}
\title{
Bayesian Latent Class Analysis variable and group selection via collapsed sampling
}
\description{
\code{blca.collapsed} performs sampling from a collapsed posterior distribution to select the  number of hidden classes G and the most relevant clustering variables if requested.
}
\usage{
blca.collapsed( X, G, formula=NULL, ncat=NULL, alpha=1, beta=1, delta=1, 
                start.vals=c("single"), counts.n=NULL, iter=5000, burn.in=500, thin=1, 
                G.sel=TRUE, var.sel=FALSE, post.hoc.run=TRUE, 
                control.post.hoc=list(iter=2000,burn.in=500,thin=1), 
                var.prob.thresh=0.75, n.gibbs=nrow(X), only.gibbs=TRUE, 
                G.max=30, G.prior=dpois(1:G.max, lambda=1)/sum(dpois(1:G.max, lambda=1)), 
                prob.inc=0.5, hprior.model=FALSE, relabel=TRUE, verbose=FALSE, 
                verbose.update=1000 )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{X}{
Data matrix or data frame. This may take one of several forms, see \code{\link{data.blca}}.
}
  \item{G}{
The number of classes to run LCA for in the case of \code{G.sel=FALSE}, otherwise this gives the initial number of groups for the sampler.
}
  \item{formula}{
Optional formula declaring variables to be used from \code{X} e.g. ``\code{~ x1 + x2}''.
}
  \item{ncat}{
	The number of categories in each variable (optional).
}
  \item{alpha, beta}{
The prior values for the data conditional on group membership. Defaults to 1, i.e, a uniform prior, for each value.
}
  \item{delta}{
Prior values for the mixture components in model.  Defaults to 1, i.e., a uniform prior.  May be single or vector valued (of length G).
}
 \item{start.vals}{
Denotes how class membership is to be assigned during the initial step of the algorithm. Initialization is always "single", which randomly assigns data points exclusively to one class.
}
  \item{counts.n}{
If data patterns have already been counted, a data matrix consisting of each unique data pattern can be supplied to the function, in addition to a vector counts.n, which supplies the corresponding number of times each pattern occurs in the data. 
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
 \item{G.sel}{
Logical, indicating whether a search over the number of components should be performed.
}
\item{var.sel}{
	Logical, indicating whether a variable selection search should be performed.
}
  \item{post.hoc.run}{
	Logical, indicating whether a post hoc run should be performed on the optimal model found from model search algorithm.
}
  \item{control.post.hoc}{
  List specifying settings for a post hoc Gibbs sampling run to get parameter estimates for most visited model if \code{G.sel=TRUE}.
}
\item{var.prob.thresh}{
 Threshold on posterior probability for inclusion in post hoc run.
}
  \item{n.gibbs}{
	The number of items to randomly sample in each sweep of the class allocation. Defaults to nrow(X), but for large datasets some random selection of rows may be preferred. 
}
  \item{only.gibbs}{
Logical, indicating if just gibbs sampling alone should be used to update labels, or to allow other moves which can reallocated a number of items at a time.
}
  \item{G.max}{
	The maximum number of components that should be searched for if \code{G.sel==TRUE}
}
  \item{G.prior}{
Prior on the number of components. Defaults to Poisson(1).
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
  \item{verbose}{
Logical, indicating whether a running update should be printed to the screen.
}
  \item{verbose.update}{
If \code{verbose==TRUE} then the interval between summaries are printed to screen. 
}

}
\details{
Model selection for both the number of groups and (optionally) variables in an LCA model. Only symmetric priors on both class weights and item probabilities are allowed when using this function. 
}
\value{
A list of class "blca.collapsed" is returned, containing:
\item{method}{The method used to fit LCA.}
\item{call}{The initial call passed to the function.}
\item{args}{List of arguments passed to the function.}
\item{classprob}{The class probabilities for the optimal model if \code{post.hoc.run == TRUE}. These are based on an auxiliary run of a Gibbs sampler.}
\item{classprob.sd}{Posterior standard deviation estimates of the class probabilities for the optimal model if \code{post.hoc.run == TRUE}. These are based on an auxiliary run of a Gibbs sampler.}
\item{itemprob}{The item probabilities, conditional on class membership for the optimal model if \code{post.hoc.run==TRUE}. These are based on an auxiliary run of a Gibbs sampler. Either a matrix or a list dependending on whether there are only dichotomous variables or a mix of di and polychotomous variables.}
\item{itemprob.sd}{Posterior standard deviation estimates of the item probabilities for the optimal model if \code{post.hoc.run == TRUE}. These are based on an auxiliary run of a Gibbs sampler. Either a matrix or a list dependending on whether there are only dichotomous variables or a mix of di and polychotomous variables.}
\item{itemprob.tidy}{Another representation of itemprob for use by plotting functions.}
\item{Z}{Estimate of class membership for each unique datapoint.}
\item{samples}{A list containing the iterates of the collapsed Gibbs sampler.}
\item{accrates}{A list giving the acceptance rates for each of the Metropolis-Hastings moves in the collapsed sampler.}
\item{iter}{The number of iterations the collapsed sampler was run for post burn-in.}
\item{burn.in}{The number of iterations the gibbs sampler was run before beginning to store values.}
\item{thin}{The thinning rate for samples from the collapsed sampler.}
\item{G.sel}{The value of \code{G.sel} from the call.}
\item{var.sel}{The value of \code{var.sel} from the call.}
\item{prob.inc}{The value of \code{prob.inc} from the call.}
\item{hprior.model}{The value of \code{hprior.model} from the call.}
\item{relabel}{The value of \code{relabel} from the call.}
\item{prior}{A list containing the prior values specified for the model.}
\item{model.indicator}{Binary vector indicating the variables which were included for clustering.}
\item{ncat}{The num.categories from the call.}
}

\references{
White, A., Wyse, J. and Murphy, T. B. (2016). Bayesian variable selection for latent class analysis using a collapsed Gibbs sampler. Statistics and Computing, volume 26, 511-527.
}
\author{
Jason Wyse, Arthur White
}

\seealso{
\code{\link{blca}}, \code{\link{blca.em}}, \code{\link{blca.vb}}, \code{\link{blca.gibbs},} 
}
\examples{
data(Alzheimer)
fit <- blca.collapsed(Alzheimer, 2, G.sel=TRUE, var.sel=TRUE, iter=2000, burn.in=500, thin=1 )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ blca }
\keyword{ collapsed }% __ONLY ONE__ keyword per line
