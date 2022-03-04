predict.blca <- function( object, newdata=NULL, ...)
{
  ncat <- object$args$ncat
  if( is.null(newdata) ) newdata <- object$args$X 
  df <- blca.prep.data( newdata, object$args$formula, NULL, ncat )
  newdata <- df$X
  
  # blca.check.missing has to be wrapped in blca.prep.data...
  if( !is.null(df$missing.idx) )
  {
    warning("Missing values encountered in X: rows with NA have been removed", call.=FALSE)
    newdata <- na.omit(newdata)
  }
  
  if(any( class(object)[1] == c("blca.em", "blca.vb") )) return( predict.blca.pt.estimate( object, newdata ))
  if(any( class(object)[1] == c("blca.gibbs", "blca.boot"))) return( predict.blca.avg.samples( object, newdata))
}

predict.blca.pt.estimate <- function( object, newdata  )
{
  # prediction function for em/vb
  # extract the parameters
  ncat <- object$args$ncat
  classprob <- object$classprob
  itemprob <- object$itemprob
  # only take non excluded columns (this needs to match classprob etc.)
  var.in <- which(object$model.indicator == 1) # this is surplus but left to avoid issues elsewhere- will resolve at later date
  X <- newdata[,var.in]
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
  sarr <- array(0,dim=c(M,N,G))
  for( k in seq_along(var.in) ) sarr[k,,] <- t( log( itemprob[[k]][, X[,var.in[k]] + 1 ] ) )
  small <- 1e-10
  a <- t( apply( sarr, c(3,2), sum)  + log(classprob + small) )
  a <- exp(a)
  rsa <- rowSums(a)
  s <- a / rsa
  rownames(s) <- rownames(X)
  return( s ) 
  
}

predict.blca.avg.samples <- function( object, newdata )
{
  # prediction function for bootstrap/gibbs 
  # extract the parameters
  ncat <- object$args$ncat
  iter <- nrow( object$samples$classprob )
  classprob <- t( object$samples$classprob )
  itemprob <- object$samples$itemprob
  # only take non excluded columns (this needs to match classprob etc.)
  var.in <- which(object$model.indicator == 1) # this is surplus but left to avoid issues elsewhere- will resolve at later date
  X <- newdata[,var.in]
  G <- object$G #length(classprob)
  N <- nrow(X)
  M <- length(var.in)
  sarr <- array(0,dim=c(M,N,G*iter))
  for( k in seq_along(var.in) ) sarr[k,,] <- t( log( itemprob[[k]][, X[,var.in[k]] + 1 ] ) ) #gf( itemprob[[k]], X[,var.in[k]] )
  small <- 1e-10
  a <- t( apply( sarr, c(3,2), sum)  + log( as.vector(classprob) + small) )
  s <- matrix( 0, nrow=N, ncol=G )
  for( b in 1:iter )
  {
    ofse <- (b-1)*G
    aa <- a[, (ofse+1):(ofse+G)]
    aa <- exp(aa)
    rsa <- rowSums(aa)
    aa <- aa / rsa
    s <- s + aa
  }
  # now avg over all the samples
  s <- s/iter
  rownames(s) <- rownames(X)
  return( s ) 
}



