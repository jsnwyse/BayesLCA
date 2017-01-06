\name{MAP}
\alias{MAP}
\alias{unMAP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Maximum \emph{a posteriori} (MAP) Classification
}
\description{
\code{MAP} obtains maximum \emph{a posteriori} (MAP) classifications. \code{unMAP} converts a classification vector into an indicator matrix.  
}
\usage{
MAP(mat, tie = c("random", "standard"))
unMAP(vec)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{
An \eqn{N \times G} matrix, typically \eqn{N} denotes observations from a dataset and \eqn{G} denotes the number of underlying groups in the data. Each row is expected to contain positive entries which sum to 1, but this isn't strictly necessary.
}
  \item{tie}{
May take one of two values, \code{"random"} or \code{"standard"}. Takes the value \code{"random"} by default. See 'Details'.
}
  \item{vec}{
An vector consisting of integer entries. \code{unMAP} is intended to be used with a vector whose entries are classifications of dataset observations.
}

}
\details{
For each row in \code{mat}, \code{MAP} assigns an indexing value identifying the entry in the row taking the highest value. In the case where multiple values in a row share a common largest value, \code{tie} determines how such a value is chosen. If \code{tie = "random"}, one of the suitable values is chosen at random; when \code{tie = "standard"}, the first such suitable value is selected, in common with other packages. Defaults to \code{"random"}. 

}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\code{MAP} returns a classification vector. \code{unMAP} returns a classification matrix, with each row indicating group membership by the column entry which is non-zero (and equal to one). 
}
\author{
Arthur White
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{Zscore}}
}
\examples{
##Simple example
s1<- sample(1:2, 10, replace=TRUE)
unMAP(s1)
MAP(unMAP(s1))

##More to the point
data(Alzheimer)
fit<- blca.em(Alzheimer, 2)
MAP(fit$Z) ## Best estimates of group membership. 

mat1<- matrix(1/3, nrow=10, ncol=3) ##demonstrating the use of "tie" argument
MAP(mat1, tie="random")
MAP(mat1, tie="standard")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ map }
\keyword{ unmap }% __ONLY ONE__ keyword per line
