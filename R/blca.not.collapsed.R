# this is the Gibbs sampler and needs to be renamed

blca.not.collapsed <- function( X , G, alpha = 1, beta = 1, delta = 1, num.categories = NULL, iter = 6000, burn.in = 1000, thin = 1, relabel = TRUE, silent = TRUE )
{

	t1 <- proc.time()

	X <- as.matrix(X)

	N <- nrow(X)
	M <- ncol(X)
	
	#check that matrix is binary and/or n.categories is passed
	if( is.null(num.categories) ){
		if( !all( X[X>0]==1) ){
			stop("\t A matrix other than binary must have non-null n.categories vector" )
		} else {
			#matrix is binary
			num.categories <- rep( 2, M )
		}
	}

	t <- apply( X, 2, min )
	if( sum(t) > 0 ) stop("\t please recode categories from 0,...,ncat  to use blca.collapsed")
	t <- apply( X, 2, max )
	if( sum( t + 1 - num.categories ) > 0 ) stop("\n number of categories in X exceeds num.categories please recode categories from 0, ..., ncat.")
	
	## safety checks ##
	
	if(length(num.categories)!= M){
		stop("\t The length of n.categories must be the same as the number of records per observation.")
	}	

	if(iter < burn.in){
		warning("\t The number of burn in iterations is greater than the number of iterations-- this will be automatically adjusted to the functions parameters.")
		iter = iter + burn.in
	}

	if((iter-burn.in)%%thin != 0)
		stop("\t Please thin by an amount that divides perfectly into iter - burn.in. ")

	stored <- (iter - burn.in) / thin

	memberships = numeric(stored*N)
	weights = numeric(stored*G)
	variable.probs = numeric(stored*G*sum(num.categories))
	log.post = numeric( iter )
	
	model.indicator <- rep(1,M) #modify this later to allow for flexible model specification

	hparam <- c( delta, beta )

	w <- .C(	"BLCA_GIBBS_SAMPLER",														as.integer(X),
				as.integer(N),																as.integer(M),
				as.integer(num.categories),													as.double(hparam),
				as.integer(G),															 	as.integer(iter),											
				as.integer(burn.in),
				as.integer(thin),															memberships = as.integer(memberships),					
				weights = as.double(weights),												variable.probs = as.double(variable.probs),			
				log.posterior = as.double(log.post),										as.integer(model.indicator),			
				PACKAGE = "BayesLCA" )

	membership.mat <- matrix( w$memberships, nrow = stored, ncol=N, byrow=FALSE ) + 1
	
	#split these into a list as it will be easier to rearrange from label processing
	var.probs.l <- list()
	
	k <- 1
	for( i in 1:M )
	{
		if(i == 1){
			gap <- 0
		}else{
			gap <- stored * G * sum( num.categories[1:(i-1)] )
		}
		#variable probabilities stacked by group and then iteration
		if( model.indicator[i] )
		{
			var.probs.l[[k]] <- matrix( w$variable.probs[(gap+1):(gap + stored*G*num.categories[i])] , nrow = stored * G, ncol=num.categories[i], byrow=TRUE )
			k <- k+1
		}
	}	
	
	weights.mat <- matrix( w$weights , nrow=stored, ncol=G, byrow=TRUE )

	if(relabel) relabelled <- undo.label.switching( membership.mat,rep(G, stored) )
	
	#post processed weights and probabilities
	
	l <- 1
	for( k in 1:stored )
	{
		weights.mat[ k, ] <- weights.mat[ k, relabelled$permutation[k,1:G] ]
		for( j in 1:M  )
		{
			it <- (k-1)*G
			if( model.indicator[j] )
			{
				var.probs.l[[l]][ it + 1:G , ] = var.probs.l[[j]][ it + relabelled$permutation[ k, 1:G ] , ]
				l <- l+1
			}
		}
	}
	
	v.probs <- list()
	tt <- seq( 0, stored-1, by=G )
	l <- 1
	for( j in 1:M )
	{
		if( model.indicator[j] )
		{
			v.probs[[l]] <- list()
			v.probs[[l]]$mean  <- matrix( nrow = G, ncol = num.categories[j] ) 
			v.probs[[l]]$se <- matrix( nrow = G, ncol = num.categories[j] )
			for( g in 1:G )
			{
				b <- var.probs.l[[l]][ tt+g , ]
				v.probs[[l]]$mean[g,] <- apply( b, 2, mean )
				v.probs[[l]]$se[g,] <- apply( b, 2, sd )
			}
			l <- l+1
		}
	}
	
	#compile the list to return

	x <- list()
	
	x$call <- match.call()
	
	x$classprob <- apply( weights.mat, 2, mean )
	x$classprob.se <- apply( weights.mat, 2, sd )
	
	x$itemprob <- v.probs
	
	x$X <- X #this is needed for the post-hoc parameter step
	x$num.categories <- num.categories
	
	x$samples <- list()
	
	x$samples$logpost <- w$log.posterior
	x$samples$Giter <- rep( G, stored )
	
	if(relabel)
	{
		x$labelstore <- relabelled$relab
		x$labelstore.permutation <- relabelled$permutation
		x$G.Z <- relabelled$components
		x$Z <- relabelled$label.prob
	}
	
	#inputs
	
	x$iter <- iter
	x$burn.in <- burn.in
	x$thin <- thin
	x$relabel <- relabel
	x$prior <- list( alpha=alpha, beta=beta, delta=delta )
	
	t2 = proc.time()
	
	ttime = t2-t1
  
  	x$time = ttime[3]	
  	
  	class(x) <- c("blca.not.collapsed", "blca")
	
	return(x)
}
