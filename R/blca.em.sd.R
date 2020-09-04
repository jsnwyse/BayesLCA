blca.em.sd <- function( X, ncat, classprob, itemprob, model.indicator, prior=NULL, small=1e-10 )
{
  G <- length( classprob )
  M <- length( ncat )
  
  # convert the class and itemprobs to cts scale
  eta <- log( classprob / classprob[G] )
  eta <- eta[1:(G-1)]
  
  zeta <- NULL
  l <- 1
  for( j in 1:M )
  {
    if( model.indicator[j] )
    {
      P <- itemprob[[l]]/itemprob[[l]][,ncat[j]]
      P <- log( P[,1:(ncat[j]-1)] )
      zeta <- c( zeta, as.vector( t(P) ) )
      l <- l + 1
    }
  }
  
  par <- c( eta, zeta )
  
  # use hessian function from numDeriv package
  obsI.par <- -hessian(func = blca.em.sd.log.post, par, ncat=ncat, G=G, M=M, X=X, model.indicator=model.indicator, prior=prior )
  obsI.par.inv <- ginv( obsI.par ) #chol2inv( chol(obsI.par) )
  
  # convert this to the obsI for theta
  # transformation matrix
  W <- matrix(0,nrow= G + sum(ncat*model.indicator)*G, ncol= G-1 + sum(ncat*model.indicator-1)*G )
  
  W[1:G,1:(G-1)] <- rbind( diag(1,nrow=G-1) * classprob[-G], numeric(G-1) ) - classprob %*% t( classprob[-G] )
  
  o <- G
  r <- G-1
  
  l <- 1
  for( j in 1:M )
  {
    if( model.indicator[j] )
    {
      for(g in 1:G )
      {
        a <- itemprob[[l]][g,]
        b <- a[ -ncat[j] ]
        W[(o+1):(o+ncat[j]),(r+1):(r+ncat[j]-1)] <- rbind( diag(1,nrow=ncat[j]-1)*b, numeric(ncat[j]-1)) - a %*% t(b)
        o <- o + ncat[j]
        r <- r + ncat[j]-1
      }
      l <- l + 1
    }
  }
  
  obsI.theta.inv <- W %*% ( obsI.par.inv %*% t(W))
  sds <- sqrt( diag(obsI.theta.inv) )
  #if( any(sds < small) ) sds[ sds < small ] <- NA
  
  # convert to list
  classprob.sd <- sds[1:G]
  o <- G
  itemprob.sd <- vector( M, mode="list" )
  itsd <- matrix( nrow=G, ncol=sum(model.indicator))
  l <- 1
  for( j in 1:M )
  {
    if( model.indicator[j] )
    {
      itemprob.sd[[l]] <- matrix( sds[(o+1):(o+G*ncat[j])], ncol=ncat[j], byrow=TRUE )
      if( !any(ncat>2) ) itsd[,l] <- itemprob.sd[[l]][,1]
      rownames( itemprob.sd[[l]] ) <- paste( "Group", 1:G )
      colnames( itemprob.sd[[l]] ) <- paste("Cat",0:(ncat[j]-1) )
      o <- o + G*ncat[j]
      l <- l + 1
    }
  }
  
  if( !any(ncat>2) ) 
  {
    itemprob.sd <- itsd
    rownames( itemprob.sd ) <- paste( "Group", 1:G )
    colnames( itemprob.sd ) <- names( itemprob )
  }
  return( list( classprob.sd = classprob.sd, itemprob.sd = itemprob.sd ))
}

blca.em.sd.log.post <- function( par, ncat, G, M, X, model.indicator, prior )
{
  
  eta <- c( par[1:(G-1)], 0 )
  classprob <- exp( eta )
  classprob <- classprob / sum(classprob)
  
  logJ <- 0
  # add Jacobian terms if needed
  if( !is.null(prior) ) 
  {
    d <- G - 1 + sum( ncat * model.indicator - 1 ) * G
    W <- matrix( 0, nrow=d, ncol=d )
    a <- classprob[1:(G-1)]
    W[1:(G-1),1:(G-1)] <- diag(a) - a %*% t(a)
  }
  
  o <- G-1
  
  itemprob <- vector( sum(model.indicator), mode="list")
  l <- 1
  for( j in 1:M )
  {
    if( model.indicator[j] )
    {
      zeta <- par[(o+1):(o+G*(ncat[j]-1))]
      u <- matrix( zeta, ncol=ncat[j]-1, byrow=TRUE )
      u <- cbind( u, numeric(G) )
      u <- exp(u)
      s <- rowSums(u)
      itemprob[[l]] <- u/s
      if( !is.null(prior) )
      {
        for( g in 1:G )
        {
          idx <- (g-1)*ncat[j] + (o+1):(o+ncat[j]-1)
          a <- itemprob[[l]][g,]
          W[idx,idx] <- diag(a) - a %*% t(a)
        }
      }
      o <- o + G * (ncat[j]-1)
      l <- l + 1
    }
  }
  
  if( !is.null(prior) ) logJ <- determinant( W, logarithm = TRUE )
  
  # if there is a prior there needs to be a Jacobian correction
  
  return( blca.compute.log.post( X, ncat, classprob, itemprob, model.indicator=model.indicator, prior, small=1e-10, reshape=FALSE ) + logJ )
}
