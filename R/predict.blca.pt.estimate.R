predict.blca.pt.estimate <- function( object, newdata = NULL )
{
  ncat <- object$args$ncat
  if( is.null(newdata) ) newdata <- na.omit( object$args$X )
  
  # blca.check.missing has to be wrapped in blca.check.data...
  missing <- blca.check.missing( newdata )$miss
  if( missing )
  {
    warning("Missing values encountered in X: rows with NA have been removed", call.=FALSE)
    newdata <- na.omit(newdata)
  }
  
  X <- blca.check.data( newdata, NULL, ncat )$X
  
  # extract the parameters
  classprob <- object$classprob
  itemprob <- object$itemprob
  # only take non excluded columns (this needs to match classprob etc.)
  var.in <- which(object$model.indicator == 1)
  X <- X[,var.in]
  G <- length(classprob)
  N <- nrow(X)
  M <- length(var.in)
  if( !any(ncat > 2) )
  {
    # all categories are binary- reshape the result matrix
    itpr <- vector( ncol(itemprob), mode="list")
    for( k in seq_along(itpr) ) itpr[[k]] <- cbind(1-itemprob[,k],itemprob[,k])
    itemprob <- itpr
  }
  
  gf <- function(P,z)
  {
    log( apply(P,1,function(zz){zz[z+1]}) )
  }
  sarr <- array(0,dim=c(M,N,G))
  for( k in seq_along(var.in) ) sarr[k,,] <- gf( itemprob[[k]], X[,var.in[k]] )
  small <- 1e-10
  a <- t( apply( sarr, c(3,2), sum)  + log(classprob + small) )
  a <- exp(a)
  rsa <- rowSums(a)
  s <- a / rsa
  rownames(s) <- rownames(X)
  return( s ) 
  
}

