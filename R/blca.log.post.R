blca.compute.log.post <- function( X, ncat, classprob, itemprob, model.indicator, prior=NULL, small=1e-10 )
{
	N <- nrow(X)
	M <- ncol(X)
	G <- length( classprob )
	S <- matrix( 0, nrow=N, ncol=G )
	
	#print( classprob )
	#print( itemprob )
	
	wch <- which(model.indicator==1)
	for( k in 1:N)
	{
		for( g in 1:G )
		{
			S[k,g] <- log( classprob[g] + small ) 
			for( j in wch )
			{
				lidx <- sum( model.indicator[1:j] )
				if( any( ncat > 2) )
				{
					S[k,g] <- S[k,g] + log( itemprob[[lidx]][ g, X[k,j] ] + small )
				}else{
					S[k,g] <- S[k,g] + X[k,j] * log( itemprob[ g, lidx ] + small ) + (1-X[k,j]) * log( 1 - itemprob[g, lidx]  + small )
				}
			}
		}
	}	
	
	# now apply the log sum exp
	mx <- apply( S, 1, max )
	
	Snew <- S - mx
	
	Snew <- exp(Snew)
	
	ap <- apply( Snew, 1, sum )
	
	L <- mx + log( ap )
	
	llike <- sum( L )
	
	
	Z <- Snew / ap
	
	if( !is.null(prior) )
	{
		alpha <- prior$alpha
		beta <- prior$beta
		delta <- prior$delta
		if( !any( ncat > 2 ) )
		{
			A <- lgamma( alpha + beta ) - lgamma( alpha ) - lgamma( beta ) + ( alpha - 1 ) * log(itemprob + small) + (beta-1)*log(1-itemprob+small)
			B <- lgamma( sum(delta) ) - sum( lgamma(delta) ) + (delta-1)*log( classprob + small )
			lprior <- sum(A) + sum(B)
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
				}
			}
			
			lprior <- 0
			
			for( g in 1:G )
			{
				for( j in 1:M )
				{
					if( j == 1 ) idxstr <- 0 else idxstr <- G * sum( ncat[1:(j-1)] )
					#take out the terms relevant to var j
					t <- gamma[ ( idxstr + (g-1) * ncat[j] +  1 ) : ( idxstr + g* ncat[j]  ) ]
					K <- lgamma( sum(t) ) - sum( lgamma(t) ) 
					lidx <- sum( model.indicator[1:j] )
					G <- ( t - 1 ) * log( itemprob[[lidx]][g,] + small )
					if( model.indicator[j] == 1 ) lprior <- lprior + K + sum(G)
				}
			}	
		
			lprior <- lprior + lgamma( sum(delta) ) - sum( lgamma(delta) ) + sum( (delta-1)*log( classprob + small ) )
		} 
			
	}else{
		lprior <- 0
	}
	
	return(  list( llike =  llike + lprior, Z=Z ) )
	
}
