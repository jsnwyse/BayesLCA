# this is the Gibbs sampler and needs to be renamed

blca.gibbs <- function( X , G, alpha = 1, beta = 1, delta = 1, num.categories = NULL, iter = 5000, burn.in = 100, thin = 1, relabel = TRUE, verbose = TRUE, verbose.update = 1000 )
{

	t1 <- proc.time()

	X <- as.matrix(X)

	N <- nrow(X)
	M <- ncol(X)
	
	#check that matrix is binary and/or n.categories is passed
	if( is.null(num.categories) ){
		if( !all( X[X>0]==1) ){
			stop("\t A matrix other than binary must have non-null n.categories vector\n" )
		} else {
			#matrix is binary
			num.categories <- rep( 2, M )
		}
	}

	t <- apply( X, 2, min )
	if( sum(t) > 0 ) stop("\t please recode categories from 0, ..., ncat-1  to use blca.gibbs\n")
	t <- apply( X, 2, max )
	if( sum( t + 1 - num.categories ) > 0 ) stop("\n number of categories in X exceeds num.categories please recode categories from 0, ..., ncat-1\n")
	
	## safety checks ##
	
	if(length(num.categories)!= M){
		stop("\t The length of num.categories must be the same as the number of records per observation\n")
	}	

	if( length(alpha) > 1 ) 
	{
		warning("\t alpha value must be a positive scalar; only first entry will be used \n")
		alpha <- alpha[1]
	}
	if( length(beta) > 1 ) 
	{
		warning("\t beta value must be a positive scalar; only first entry will be used \n")
		beta <- beta[1]
	}
	if( length(delta) > 1 ) 
	{
		warning("\t delta value must be a positive scalar; only first entry will be used \n")
		delta <- delta[1]
	}

	#if(iter < burn.in){
	#	warning("\t The number of burn in iterations is greater than the number of iterations-- this will be automatically adjusted to the functions parameters.")
	#	iter = iter + burn.in
	#}

	#if((iter-burn.in)%%thin != 0)
	#	stop("\t Please thin by an amount that divides perfectly into iter - burn.in. ")

	stored <- iter / thin

	memberships <- numeric(stored*N)
	weights <- numeric(stored*G)
	variable.probs <- numeric(stored*G*sum(num.categories))
	log.post <- numeric( stored )
	
	model.indicator <- rep(1,M) #modify this later to allow for flexible model specification

	hparam <- c( delta, beta )

	w <- .C(	"BLCA_GIBBS_SAMPLER",													as.integer(X),
				as.integer(N),																as.integer(M),
				as.integer(num.categories),											as.double(hparam),
				as.integer(G),															 	as.integer(iter),											
				as.integer(burn.in),
				as.integer(thin),															memberships = as.integer(memberships),					
				weights = as.double(weights),											variable.probs = as.double(variable.probs),			
				as.integer(model.indicator), 											log.posterior = as.double(log.post),
				as.integer( verbose ),													as.integer( verbose.update ),			
				PACKAGE = "BayesLCA" )
	
	if( verbose ) cat("\nFinished sampling...\n")

	membership.mat <- matrix( w$memberships, nrow = stored, ncol=N, byrow=FALSE ) + 1
	
	#split these into a list as it will be easier to rearrange from label processing
	var.probs.l <- list()
	
	l <- 1
	for( j in 1:M )
	{
		if(j == 1){
			gap <- 0
		}else{
			gap <- stored * G * sum( num.categories[1:(j-1)] )
		}
		#variable probabilities stacked by group and then iteration
		if( model.indicator[j] )
		{
			var.probs.l[[l]] <- matrix( w$variable.probs[(gap+1):(gap + stored*G*num.categories[j])] , nrow = stored * G, ncol=num.categories[j], byrow=TRUE )
			l <- l+1
		}
	}	
	
	weights.mat <- matrix( w$weights , nrow=stored, ncol=G, byrow=TRUE )

	if( relabel ) relabelled <- undo.label.switching( membership.mat,rep(G, stored) )
	
	#post processed weights and probabilities
	
	for( k in 1:stored )
	{
		weights.mat[ k, ] <- weights.mat[ k, relabelled$permutation[k,1:G] ]
		l <- 1
		for( j in 1:M  )
		{
			it <- (k-1)*G
			if( model.indicator[j] )
			{
				var.probs.l[[l]][ it + 1:G , ] = var.probs.l[[l]][ it + relabelled$permutation[k,1:G] , ]
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
			rownames( v.probs[[l]]$mean ) <- paste( "Group", 1:G )
			colnames( v.probs[[l]]$mean ) <- paste("Cat",0:(num.categories[j]-1) )
			rownames( v.probs[[l]]$se ) <- paste( "Group", 1:G )
			colnames( v.probs[[l]]$se ) <- paste("Cat",0:(num.categories[j]-1) )
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
	
	Dbar<- mean(w$log.posterior)
	S2<- var(w$log.posterior)
	
	x$DIC<- 2*(2*Dbar - x$samples$logpost)
	# need to check definition of BICM
	#x$BICM<- 2*(Dbar - S2*(log(sum(counts.n))-1))
	x$AICM<- 2*(Dbar - S2)
	
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
