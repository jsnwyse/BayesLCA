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
		
	}else{
		
		# check here if X is a data frame
		
		N <- nrow(X)
		M <- ncol(X)
		
		Y <- matrix( nrow=N, ncol=M )
		colnames(Y) <- colnames(X) 
		
		
		if( is.data.frame(X) )
		{
			isfact <- unlist( lapply( X, is.factor ) )
			factors <- which( isfact == TRUE ) 
		
			for( k in factors )
			{
				warning("X is a data frame: factor levels converted to integers")
				levs <- levels( X[,k] )
				for( j in 1:length(levs) ) 
				{
				  idx <- which( X[,k] == levs[j] )  
				  Y[ idx, k ] <- j-1
				}
			}
			
			notfactors <- which( isfact == FALSE )
			for( k in notfactors) Y[ , k ] <- as.numeric( X[ , k ] )
		
		}else{
			if( is.matrix(X) ) Y <- X
		}
	
		X <-Y
	}
	

	#check that matrix is binary and/or n.categories is passed
	if( is.null(ncat) ){
		if( !all( X[(!is.na(X)) & X>0]==1) ){
			stop("A matrix other than binary must have non-null ncat vector" )
		} else {
			#matrix is binary
			ncat <- rep( 2, M )
		}
	}else{
		if(length(ncat)!= M) stop("The length of ncat must be the same as the number of records per observation\n")
	}
	
	t <- apply( X, 2, function(x) max( x, na.rm=TRUE ) )
	if( sum( t + 1 - ncat ) > 0 ) stop("Number of categories in X exceeds ncat please recode categories from 0, ..., num cat-1")
	
	t <- apply( X, 2, function(x) min( x, na.rm=TRUE ) )
	if( sum(t) > 0 ) warning("It appears some variables have empty category 0. Analysis will proceed based assuming this is the case. If not, check X input.")
	
	return( list( X=X, ncat=ncat ) )

}

