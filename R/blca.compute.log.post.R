blca.compute.log.post <- function( X, ncat, classprob, itemprob, model.indicator, prior=NULL, small=1e-10, reshape=TRUE )
{
  N <- nrow(X)
  M <- ncol(X)
  G <- length( classprob )
  S <- matrix( 0, nrow=N, ncol=G )
  hparam <- c(1,1)
  # need something here for the binary case.--
  if( any(ncat>2) | !reshape )
  {
    variable.probs <- unlist( lapply(itemprob, t) )
  }else{
    variable.probs <- vector( ncol(itemprob), mode="list")
    for( k in seq_along(variable.probs) ) variable.probs[[k]] <- cbind(1-itemprob[,k],itemprob[,k]) 
    itemprob <- variable.probs
    variable.probs <- unlist( lapply(variable.probs, t) )
  }
    #ncat  <- unlist( lapply(itemprob,ncol) )
  log.like <- 0
  x.vec <- as.vector( t(X) ) 
  
  w.ll <- .C(	"BLCA_LOG_LIKE",											        as.integer(x.vec),
              as.integer(N),																as.integer(M),
              as.integer(ncat),															as.double(hparam),
              as.integer(G),															 				
              weights = as.double(classprob),							  variable.probs = as.double(variable.probs),			
              as.integer(model.indicator), 									value = as.double(log.like),
              PACKAGE = "BayesLCA" )
  
  llike <- w.ll$value
  
  
  if( !is.null(prior) )
  {
    
    out.prior <- blca.check.prior( prior$alpha, prior$beta, prior$delta, G, M, ncat )
    gamma.r <- out.prior$gamma
    delta <- out.prior$delta

    lprior <- 0
    
    for( g in 1:G )
    {
      for( j in 1:M )
      {
        if( j == 1 ) idxstr <- 0 else idxstr <- G * sum( ncat[1:(j-1)] )
        #take out the terms relevant to var j
        t <- gamma.r[ ( idxstr + (g-1) * ncat[j] +  1 ) : ( idxstr + g * ncat[j]  ) ]
        K <- lgamma( sum(t) ) - sum( lgamma(t) ) 
        lidx <- sum( model.indicator[1:j] )
        # needs to be modified for the two category case
        lpr <- ( t - 1 ) * log( itemprob[[lidx]][g,] + small )
        if( model.indicator[j] == 1 ) lprior <- lprior + K + sum(lpr)
      }
    }	
    
    lprior <- lprior + lgamma( sum(delta) ) - sum( lgamma(delta) ) + sum( (delta-1)*log( classprob + small ) )

  }else{
    lprior <- 0
  }
  
  return( llike + lprior )
  
}
