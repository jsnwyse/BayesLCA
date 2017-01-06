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
Five plots are selectable: a plot summarising item and class probability, a mosaic plot representing classification uncertainty, item probability density estimates, conditional class probability density estimates, and a diagnostics plot. The default setting is for the first four plots to be displayed, with the exception of plot.blca.em, which cannot produce density plots and so only produces the first two plots by default.

}
\usage{
\method{plot}{blca}(x, which = 1L, main = "", col1 = heat.colors(12), \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class \code{\link{blca}}.
}
  \item{which}{
Which plots to select. May be any subset of \code{1:5}, with some exceptions. See `Details'.  
}
  \item{main}{
An overall title for the plot: see \code{\link[graphics]{title}}.
}
  \item{col1}{
Specifies a list of colours to be used by the heat map plot used when \code{which = 1}.\code{\link[graphics]{image}}. Uses \code{\link[grDevices]{heat.colors}} by default, but several other choices are available. See the help files of \code{\link[fields]{image.plot}}, \code{\link[graphics]{image}} and \code{\link[grDevices]{palette}} for details.
}
  \item{\dots}{
Further arguments to be passed onto the plotting devices. When \code{which = 1}, the plotting device is \code{\link[fields]{image.plot}},  \code{\link[graphics]{mosaicplot}} in the case of \code{which=2}, and when \code{which=3:5}, \code{\link[graphics]{plot}}.
}
}
\details{
Not all plots are available for some object classes. If the object is of class \code{\link{blca.em}}, density plots (\code{which = 3:4}) are unavailable, and a warning is returned. Similarly, diagnostic plots (\code{which = 5}) for \code{\link{blca.boot}} objects are unavailable. 

The available diagnostic plots differ depending on the class of the object in question. For \code{\link{blca.em}} and \code{\link{blca.vb}} objects, the plot is intended as visual aid to check whether the respective algorithms have converged, i.e., that the log-posterior or lower bound have ceased increasing after successive iterations. The main aim of the diagnostic plot for \code{\link{blca.gibbs}} objects is to visually check diagnostic measures such as mixing and burn-in, and also to assess whether label-switching has occurred, or been corrected for satisfactorily.

Currently, the colors used in a plot can only be specified directly for \code{which = 1}. For classification uncertainty (\code{which = 2}) and density plots (\code{which = 3:4}), each group is colored by the \code{\link[grDevices]{palette}} function so that Group g takes color \code{palette()[g+1]}. For the default settings, Group 1 will then be colored red, Group 2 green, and so on.
}


\author{
Arthur White
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
plot(fit, which = 1:2) ## Parameter summary and classification uncertainty plots.

palette(rainbow(6)) ## Change color scheme
plot(fit, which = 2) 
palette("default") ## Restore default color scheme

fit2<- blca.vb(x,2)
par(mfrow = c(3,4))
plot(fit2, which = 3) ## Approximate density plots for item probability parameters.
par(mfrow = c(1,1))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ blca }
\keyword{ plot }% __ONLY ONE__ keyword per line
