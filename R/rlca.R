rlca <-
  function(n, itemprob=0.5, classprob=1, fit=NULL, ncat=NULL, return.class=FALSE ){ 
    
    if( any(classprob < 1E-10) ) stop("a class probability of zero was encountered, reduce the number of groups to reflect this")
    
    # when a fitted model is not given
    if( is.null(fit) ) 
    {
      # determine category structure from itemprob
      if( is.list( itemprob) )
      {
        # check
        Gv <- unlist( lapply(itemprob,nrow) )
        # check all have same number of rows
        if( any(Gv!=Gv[1]) ) stop("itemprob does not have the same row dimension (groups) across slots")
        if( Gv[1] != length(classprob) ) stop("itemprob row dimension (groups) does not match length of classprob")
        G <- Gv[1]
        rs <- unlist( lapply(itemprob,rowSums) )
        if( any( abs(rs-1) > 1E-8 ) ) stop("each row must sum to 1 in slots of itemprob")
        M <- length(itemprob)
        ncat <- unlist( lapply(itemprob,ncol) )
      }else{
        if( !is.matrix(itemprob) ) stop("itemprob should be either a list or a matrix")
        if( any(itemprob > 1) | any(itemprob < 0)) stop("entries of itemprob matrix are probabilities and should be between 0 and 1")
        G<- nrow(itemprob)
        M<- ncol(itemprob)
        ncat <- rep(2,M)
      }
    }else{
      # fitted model is given
      if( !any( class(fit) == "blca" ) ) stop("fit should be a BayesLCA fitted model object for random generation")
      itemprob <- fit$itemprob
      classprob <- fit$classprob
      G <- length(classprob)
      if( is.list(itemprob) )
      {
        # polychotomous case
        M <- length(itemprob)
        ncat <- unlist( lapply(itemprob,ncol) )
      }else{
        M <- ncol(itemprob)
        ncat <- rep(2,M)
      }
    }
    
    if( !any(ncat>2) )
    {
      # binary only 
      x <- matrix( nrow=n, ncol=M )
      z <- sample( 1:G, size=n, replace=TRUE, prob=classprob )
      P <- t( itemprob[ z, ] )
      p <- as.vector( P )
      x <- matrix( rbinom( n*M, size=1, prob=p ), nrow=n, byrow=T )
      if( return.class )
      {
        # a bit repetitive but... 
        tr.itemprob <- t(itemprob)
        gpr <- function(q)
        {
          a <- log( (tr.itemprob^q)*((1-tr.itemprob)^(1-q)) )
          a <- log(classprob) + colSums(a)
          a <- a - max(a)
          a <- exp(a)
          return(a/sum(a))
        }
        s <- t( apply( t(x), 2, gpr ) )
      }
      
    }else{
      
      x <- matrix( nrow=n, ncol=M )
      z <- sample( 1:G, size=n, replace=T, prob=classprob )
      # function to apply for faster sampling
      sf <- function(P,grp,ngrp)
      {
        sample( 0:(ncol(P)-1), size=ngrp, replace=T, P[grp,]  )
      }
      # generate the data conditional on the labels
      for( g in 1:G )
      {
        idx <- which( z == g )
        ng <- length( idx )
        x[idx,] <- matrix( unlist(lapply( itemprob, sf, grp=g, ngrp=ng)) , nrow=ng )
      }
      # if the class labels are being returned, also return their generation probability from model
      if( return.class )
      {
        gf <- function(P,z)
        {
          log( apply(P,1,function(zz){zz[z+1]}) )
        }
        sarr <- array(0,dim=c(M,n,G))
        for( m in 1:M ) sarr[m,,] <- gf( itemprob[[m]], x[,m] )
        a <- t( apply( sarr, c(3,2), sum)  + log(classprob) )
        a <- exp(a)
        rsa <- rowSums(a)
        s <- a / rsa
      }
    }
    
    colnames(x) <- paste0("X",1:M)
    
    if( return.class ) x <- list( X=x, class=z, class.prob.memb=s ) 
    class(x) <- "blca.rand"
    return(x)
    
  }
