blca.cat.to.indicators <- function( X, ncat )
{
  if( is.data.frame(X) ) stop("argument X must be a matrix.")
  
  W <- matrix(0, nrow=nrow(X), ncol=sum(ncat) )
  
  coff <- cumsum(ncat)
  coff <- c(0,coff[-length(ncat)])
  
  # add one to X to shift the zero indexing to one offset
  idx <- t( t(X+1) + coff )
  
  for( k in 1:nrow(idx) ) W[ k, idx[k,] ] <- 1
  
  return(W)
}

