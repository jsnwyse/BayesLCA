blca.check.data <- function( X, counts.n, ncat )
{
  
  if( class(X) == "data.blca" || !is.null(counts.n) )
  {
    if( class(X) == "data.blca" ) 
    {
      z <- X$counts.n
      Y <- X$data 
    }else{ 
      z <- counts.n
      Y <- as.matrix(X)
    }
    
    U <- matrix( rep( Y[1,], z[1] ) , nrow=z[1], byrow=TRUE )
    for( k in 2:length(z) )
    {
      U <- rbind( U, matrix( rep( Y[k,], z[k] ) , nrow=z[k], byrow=TRUE ) )
    }
    X <- as.matrix(U)
    
    N <- nrow(X)
    M <- ncol(X)
    
    levnames <- vector( M, mode="list")
    for( k in 1:M ) levnames[[k]] <- paste0("Cat ",0:1)
    
  }else{
    
    # check here if X is a data frame
    
    N <- nrow(X)
    M <- ncol(X)
    
    Y <- matrix( nrow=N, ncol=M )
    colnames(Y) <- colnames(X) 
    
    levnames <- vector( M, mode="list" )
    
		if( is.data.frame(X) )
		{
		  
		  # do character variables first and cast as factors
		  ischr <- unlist( lapply( X, is.character) )
		  chrs <- which( ischr == TRUE )
		  
		  for( k in chrs ) X[,k] <- factor( X[,k] )
		  
			isfact <- unlist( lapply( X, is.factor ) )
			factors <- which( isfact == TRUE ) 
		
			for( k in factors )
			{
				#warning("X is a data frame: factor levels converted to integers")
				levs <- levels( X[,k] )
				levnames[[k]] <- levs
				for( j in 1:length(levs) ) 
				{
				  idx <- which( X[,k] == levs[j] )  
				  Y[ idx, k ] <- j-1
				}
			}
			
			notfactors <- which( isfact == FALSE )
			for( k in notfactors) 
			{
			  u <- unique( X[,k] )
			  Y[ , k ] <- as.numeric( X[ , k ] - min(u) ) # this will map the lowest back to zero
			  levnames[[k]] <- paste0("Cat ",sort(u))
			}
		}else{
			if( is.matrix(X) ) Y <- X
			for( k in 1:M ) levnames[[k]] <- paste0("Cat ",0:(length(unique(X[,k])) - 1) )
		}
	
		X <-Y
	}
	

	#check that matrix is binary and/or ncat is passed
	if( is.null(ncat) ){
	  # determine ncat from the levnames vector
	  ncat <- unlist( lapply( levnames, length ) )
	}else{
		if(length(ncat)!= M) stop("the length of ncat must be the same as the number of variables")
	  if( any( ncat != unlist( lapply(levnames,length) ) ) ) stop("if passing ncat, the number of levels within variables must correspond to ncat entries")
	}
	
	return( list( X=X, ncat=ncat, levnames = levnames ) )
}


