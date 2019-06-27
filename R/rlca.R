rlca <-
function(n, itemprob=0.5, classprob=1, fit=NULL, ncat=NULL, return.class=FALSE ){ 

	if( is.null(ncat) || ( !is.null(ncat) && !any(ncat>2) ) )
	{
	
		#if ncat not given, assume binary
	
		if(is.null(fit))
		{
			itemprob<- as.matrix(itemprob)
			G<- nrow(itemprob)
			M<- ncol(itemprob)
			} else{
				itemprob<- fit$itemprob
				classprob<- fit$classprob
				G<- nrow(itemprob)
				M<- ncol(itemprob)
			}

	  x <- matrix( nrow=n, ncol=M )
	  
	  z <- sample( 1:G, size=n, replace=TRUE, prob=classprob )
	  
	  P <- t( itemprob[ z, ] )
	  
	  p <- as.vector( P )
	  
	  x <- matrix( rbinom( n*M, size=1, prob=p ), nrow=n, byrow=T )
	  
		#x <- matrix(runif(n*M), nrow=n)
		#classvec<- as.vector(rmultinom(1, n, prob=classprob))

		#ind<- c(0, cumsum(classvec))

		#for(g in 1:G) x[(ind[g]+1):ind[g+1],] <- t(t(x[(ind[g]+1):ind[g+1],]) < itemprob[g,])*1
	
	}else{
	
		if( is.null(fit) )
		{
			if( is.null(ncat) ) stop("Either fit or ncat must be specified")
			if( typeof(itemprob) != "list" ) stop("For varying numbers of categories (>2) specify itemprob as a list of G times ncat matrices.")
			G <- length(classprob)
			M <- length(ncat)
		}else{
			itemprob <- fit$itemprob
			classprob <- fit$classprob
			ncat <- fit$ncat
			G <- nrow( itemprob[[1]] )
			M <- length(ncat)
		}
		
		x <- matrix( nrow=n, ncol=M )
		z <- sample( 1:G, size=n, replace=T, prob=classprob )
		
		for( j in 1:M )
		{
			for( k in 1:n )
			{
				x[ k, j ] <- sample( 0:(ncat[j]-1), size=1, prob=itemprob[[j]][ z[k], ]  )
			}
		}
		
	}

  if( return.class ) x <- list( X=x, class=z ) 
  
	return(x)

}
