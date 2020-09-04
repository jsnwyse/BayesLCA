\name{plot.blca}
\alias{plot.blca}
\alias{plot.blca.boot}
\alias{plot.blca.em}
\alias{plot.blca.gibbs}
\alias{plot.blca.vb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot Parameter Summaries, Density Estimates and Model Diagnostics for Bayesian Latent Class Analysis 
}
\description{
Plot summarising item and class probability in the fitted LCA model.
}
\usage{
\method{plot}{blca}(x, main = "", variables = NULL, col1 = NULL, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class \code{\link{blca}}.
}
  \item{main}{
An overall title for the plot: see \code{\link[graphics]{title}}.
}

\item{variables}{
 If there are more than two levels in some variables, passing a vector of column numbers in data frame or variable names will plot only those variables.
}
  \item{col1}{
Specifies a list of colours to be used when plotting.\code{\link[graphics]{image}}. Uses \code{\link[grDevices]{heat.colors}} by default, but several other choices are available. See the help files of \code{\link[fields]{image.plot}}, \code{\link[graphics]{image}} and \code{\link[grDevices]{palette}} for details.
}
  \item{\dots}{
Further arguments to be passed onto the plotting devices. 
}
}
\details{
Creates a plot showing class size and item probabilities. For completely dichotomous items, this is a mosaic plot. For polychotomous items, side-by-side bars where the width of the bars is proportional to the class side and the height is coloured proportionally to the estimated item probabilities. 
}


\author{
Arthur White, Jason Wyse
}

\references{
Arthur White, Thomas Brendan Murphy (2014). BayesLCA: An R Package for Bayesian Latent Class Analysis." Journal of Statistical Software, 61(13), 1-28. URL: http://www.jstatsoft.org/v61/i13/.

}

\seealso{
\code{\link[fields]{image.plot}}, \code{\link[graphics]{mosaicplot}}.
}
\examples{
type1 <- c(0.8, 0.8, 0.2, 0.2)
type2 <- c(0.2, 0.2, 0.8, 0.8)
x<- rlca(1000, rbind(type1,type2), c(0.6,0.4))

fit <- blca.em(x, 2)
plot(fit) ## Parameter summary plot.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ blca }
\keyword{ plot }% __ONLY ONE__ keyword per line
