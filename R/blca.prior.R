blca.check.prior <- function( alpha, beta, delta, G, M, ncat )
{
  
  prior.init.type <- 1
  
  # delta is the prior on the group weights
  if( length(delta) == 1 ) delta<- rep(delta,G)
  if( length(delta) != G ) stop("delta prior is wrong length (length is not either 1 or G)")
  
  # alpha either acts only as the first category prior for binary, or the entire prior for multicategory case
  if( !is.matrix(alpha) && !any(ncat > 2) ){
    if( any(length(alpha)==c(1,G)) ){
      alpha <- matrix(alpha,G,M)
    }else{
      if(length(alpha)==M){
        # this is ok as there is a different alpha for each variable but the same across groups
        alpha <- matrix( alpha, G, M, byrow=TRUE)
      } else {
        if( length(alpha) == G*sum(ncat) && beta==1 ){
          if( !any( ncat > 2 ) ) warning("using only entries in alpha to assign prior for sampling")
          # pass alpha directly, this looks ok
          prior.init.type <- 2
        }else stop("item probability prior improperly specified")
      }
    }
  } 
  
  if( is.matrix(beta) && any( ncat > 2 ) ) stop("item probability prior improperly specified: please use alpha to specify the prior." )
  
  # beta either acts only as the first category prior for binary 
  if( !is.matrix(beta) ){
    if( any(length(beta)==c(1,G)) ){
      beta <- matrix(beta,G,M)
    }else{
      if(length(beta)==M){
        beta <- matrix( beta, G, M, byrow=TRUE )
      } else {
        if( length(beta) == G*sum(ncat) && alpha==1 ){
          stop("item probability prior improperly specified: for varying numbers of categories, use alpha to specify the prior")
        }else stop("item probability prior improperly specified")
      }
    }
  } 
  
  if( !any( ncat > 2 ) )
  {
    if( !all(dim(alpha) == c(G,M)) ) stop("alpha has not been specified with the correct dimensions")
    if( !all( dim(beta) == c(G,M) ) ) stop("beta has not been specified with the correct dimensions")
    # now restack the alpha and beta matrices into compatible format
    gamma <- matrix( nrow=2*G, ncol=M )	
    for( k in 1:G ) 
    {
      gamma[ 2*(k-1) + (1:2) ,  ] <- rbind( alpha[k,], beta[k,] )
    }
    prior.init.type <- 2
  }else{
    if( length(alpha) == 1 ){
      prior.init.type <- 1
      gamma <- rep( alpha, G*sum(ncat) )
    }else{
      if( length(alpha) == sum(ncat) ) 
      {
        gamma <- rep( alpha, G )
      }else if( length(alpha) == G*sum(ncat) ){
        gamma <- alpha
      }else{
        stop("alpha provided is not of a compatible length: please check and rerun")
      }
    }
  }
  
  return( list( delta=delta, gamma=gamma, prior.init.type=prior.init.type ) )
  
}

blca.check.prior.collapsed <- function( alpha, beta, delta, G.sel, G, G.max, M, ncat, only.gibbs )
{
  prior.init.type <- 1
  
  # delta is the prior on the group weights
  if( length(delta) == 1 ) delta<- rep(delta,G.max)
  if( length(delta) != G.max ) stop("delta prior is wrong length (i.e., != 1 or G.max)")
  if( length(delta) == G.max && any( delta != delta[1] ) && G.sel )
    stop("unsymmetric priors on group weights are not allowed when G.sel is set to TRUE: rerun with delta set to a scalar")
  
  # alpha either acts only as the first category prior for binary, or the entire prior for multicategory case
  if( !is.matrix(alpha) && !any( ncat > 2 ) )
  {
    if( any(length(alpha)==c(1,G.max)) )
    {
      if( length(alpha) == G.max && G.sel && any( alpha != alpha[1] ) ) stop("the prior setting for alpha is required to be the same over all groups when G.sel is set to TRUE: rerun with alpha set to a scalar") 
      alpha <- matrix(alpha,G.max,M)
      if( G.sel ) alpha <- matrix( alpha[1], G.max, M ) 
    }else{
      if(length(alpha)==M)
      {
        # this is ok as there is a different alpha for each variable but the same across groups
        alpha <- matrix( alpha, G.max, M, byrow=TRUE)
      } else {
        if( length(alpha) == G.max*sum(ncat) && beta==1 ){
          if( !any( ncat > 2 ) ) warning("using only entries in alpha to assign prior for sampling")
          if( G.sel  ) 
          {
            # here we need a thorough check to ensure that the specification for the prior weight on a category does not change over groups
            u <- rep(FALSE, M )
            for( m in 1:M )
            {
              if( m == 1 ) idx0 <- 0 else idx0 <- G.max * sum( ncat[1:(m-1)] )
              v <- alpha[ (idx0 + 1) : (idx0 + G.max*ncat[m]) ]
              if( any( v != v[1] ) ) u[m] <- TRUE
            }
            if( any(u) ) stop("the prior setting for alpha is required to be the same over all groups when G.sel is set to TRUE: rerun with alpha set to a scalar") 
          }
          # pass alpha directly, this looks ok
          prior.init.type <- 2
        }else stop("item probability prior improperly specified")
      }
    }
  } 
  
  if( is.matrix(beta) && any( ncat > 2 ) ) stop("item probability prior improperly specified: please use alpha to specify the prior" )
  
  # beta either acts only as the first category prior for binary 
  if( !is.matrix(beta)  && !any( ncat > 2 ) )
  {
    if( any(length(beta)==c(1,G.max)) )
    {
      if( length(beta) == G.max && G.sel && any(  beta != beta[1] ) ) 
        stop("unsymmetric priors with G.sel set to TRUE is not available for collapsed sampler: rerun with beta set to a scalar")
      beta <- matrix(beta,G.max,M)
      if( G.sel ) beta <- matrix( beta[1], G.max, M )
    }else{
      if(length(beta)==M)
      {
        if( G.sel && any( beta != beta[1] ) ) 
          stop("unsymmetric priors with G.sel set to TRUE is not available for collapsed sampler: rerun with beta set to a scalar") 
        beta <- matrix(beta,G.max,M, byrow=TRUE)
      } else {
        if( length(beta) == G.max*sum(ncat) && alpha==1 )
        {
          stop("item probability prior improperly specified: for varying numbers of categories use alpha to specify the prior")
        }else stop("item probability prior improperly specified")
      }
    }
  } 
  
  if( !any( ncat > 2 ) )
  {
    if( !all(dim(alpha) == c(G.max,M)) ) stop("alpha has not been specified with the correct dimensions")
    if( !all( dim(beta) == c(G.max,M) ) ) stop("beta has not been specified with the correct dimensions")
    # if the priors are unsymmetric in all instances, disable moves 1, 2 and 3
    if( any( alpha != alpha[1] ) || any( beta != beta[1] ) )
    {
      warning("Gibbs sampling is the only label update move carried out in blca.collapsed due to identifiability of priors")
      only.gibbs <- TRUE
    }
    # now restack the alpha and beta matrices into compatible format
    gamma <- matrix( nrow=2*G.max, ncol=M )	
    for( k in 1:G.max ) 
    {
      gamma[ 2*(k-1) + (1:2) ,  ] <- rbind( alpha[k,], beta[k,] )
    }
    prior.init.type <- 2
  }else{
    if( length(alpha) == 1 ){
      prior.init.type <- 1
      gamma <- rep( alpha, G.max*sum(ncat) )
    }else{
      if( length(alpha) == sum(ncat) ) 
      {
        gamma <- rep( alpha, G.max )
      }else if( length(alpha) == G.max*sum(ncat) ){
        gamma <- alpha
      }else{
        stop("alpha provided is not of a compatible length. Please check and rerun.")
      }
    }
  }
  
  return( list( delta=delta, gamma=gamma, prior.init.type=prior.init.type, only.gibbs=only.gibbs ) )
  
}
