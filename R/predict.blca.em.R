predict.blca.em <- function( object, newdata = NULL )
{
  
  ncat <- object$args$ncat
  G <- object$args$G
  
  if( is.null(newdata) ) 
  {
    # if new data is null, return predicted classes for 
    #  data passed to fit the function
    X <- object$args$X
  }else{
    X <- blca.check.data( newdata, NULL, ncat )$X
  }
  
  W <- blca.cat.to.indicators( X, ncat )
  
  class.probs <- object$classprob
  # stack the variable probabilities into an sum(ncat) by G matrix
  if( !any( ncat > 2 ) )
  {
    M <- ncol( object$itemprob )
    var.probs <- as.vector( object$itemprob )
    var.probs <- rbind( var.probs, 1-var.probs )
    vp <- matrix( nrow= sum(ncat), ncol=G ) #var.probs[,1:G]
    for( k in 1:M ) vp[ ((k-1)*2 + 1):(2*k), ] <- var.probs[, ((k-1)*G+1):(k*G) ] 
    var.probs <- vp
  }else{
    var.probs <- matrix( unlist(object$itemprob), ncol=G, byrow=T )
  }
  
  log.var.probs <- log(var.probs)
  
  # compute unnormalised log probability of group membership
  log.pr <- class.probs + t( W %*% log.var.probs )
  
  max.val <- apply( log.pr, 2, max )
  log.pr <- sweep( log.pr, 2, max.val, "-")
  
  log.pr <- exp(log.pr)
  z <- apply( log.pr, 2, sum )
  log.pr <- sweep( log.pr, 2, z, "/")
  
  return( t(log.pr) ) 
  
}

