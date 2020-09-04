\name{rlca}
\alias{rlca}

\title{
Randomly Generate Data from a Latent Class Analysis model
}
\description{
A function which randomly generates data with respect to some underlying latent class. Data may be generated either by specifying item and class probabilities or by utilising an object previously fitted to data.
}
\usage{
rlca(n, itemprob = 0.5, classprob = 1, fit = NULL, ncat=NULL, return.class=FALSE )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
Number of data points to be generated.
}
  \item{itemprob}{
The item probabilities, conditional on class membership. Defaults to 0.5. In the case of binary data, this argument can be passed as a matrix with row dimension G and column dimension equal to the number of variables. Each entry gives the probability of observing a success. In the case where any of the variables are encoded by more than two categories (or different numbers of categories), this argument should be passed as a list of matrices,  a matrix for each variable, with the matrix in the corresponding slot in the list having G rows and number of columns given by the number of categories in that variable.
}
  \item{classprob}{
The class probabilities. Defaults to 1, i.e., a one class model. 
}
  \item{fit}{
An object of class blca. If fit is supplied, data is generated using the class and item probabilities obtained. Defaults to NULL.
}
\item{ncat}{
	The number of categories in each variable (column).
}
\item{return.class}{
Logical, indicating whether the class labels and probabilities should be returned with the simulated data. 
}

}

\value{
If \code{return.class} is FALSE, then an n by \code{length(ncat)} matrix containing the data. If \code{return.class} is true, a slot with entry \code{X} equal to the data, a slot \code{class} with the labels and a slot \code{class.memb.prob} which gives the a matrix of the probabilities of having generated the rows of \code{X} from each component using the \code{classprob} and \code{itemprob} values.
}


\author{
Arthur White, Jason Wyse
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{data.blca}}
}
\examples{
type1 <- c(0.8, 0.8, 0.2, 0.2)
type2 <- c(0.2, 0.2, 0.8, 0.8)
x<- rlca(1000, rbind(type1,type2), c(0.6,0.4))

fit <- blca.em(x, 2)

x2<- rlca(1000, fit=fit)
fit2<- blca.em(x2,2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{blca }
\keyword{ random }% __ONLY ONE__ keyword per line
