blca.gibbs <- function( X, G, ncat=NULL,  alpha=1, beta=1, delta=1, start.vals=c("prior","single","across"), counts.n=NULL, model.indicator=NULL, iter=5000, burn.in=100, thin=1, accept=thin, relabel=TRUE, verbose=TRUE, verbose.update=1000 ) 
{

	t1 <- proc.time()

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
	}else{
		X <- as.matrix(X)
	}

	N <- nrow(X)
	M <- ncol(X)
	
	#check that matrix is binary and/or n.categories is passed
	if( is.null(ncat) ){
		if( !all( X[X>0]==1) ){
			stop("A matrix other than binary must have non-null ncat vector" )
		} else {
			#matrix is binary
			ncat <- rep( 2, M )
		}
	}

	t <- apply( X, 2, min )
	if( sum(t) > 0 ) stop("Please recode categories from 0, ..., ncat-1  to use blca.collapsed")
	t <- apply( X, 2, max )
	if( sum( t + 1 - ncat ) > 0 ) stop("Number of categories in X exceeds ncat please recode categories from 0, ..., num cat-1")
		

	## safety checks ##
	
	if(length(ncat)!= M){
		stop("The length of ncat must be the same as the number of records per observation\n")
	}	
	
	prior.init.type <- 1
	
	# delta is the prior on the group weights
	if( length(delta) == 1 ) delta<- rep(delta,G)
	if( length(delta) != G ) stop("delta prior is wrong length (i.e., != 1 or G)")
	
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
			 		if( !any( ncat > 2 ) ) warning("Using only entries in alpha to assign prior for sampling")
			 		# pass alpha directly, this looks ok
			 		prior.init.type <- 2
			 	}else stop("Item probability prior improperly specified.")
			}
		}
	} 
	
	if( is.matrix(beta) && any( ncat > 2 ) ) stop("Item probability prior improperly specified. Please use alpha to specify the prior." )
	
	# beta either acts only as the first category prior for binary 
	if( !is.matrix(beta) ){
	if( any(length(beta)==c(1,G)) ){
		beta <- matrix(beta,G,M)
	}else{
		if(length(beta)==M){
			beta <- matrix( beta, G, M, byrow=TRUE )
			} else {
			 	if( length(beta) == G*sum(ncat) && alpha==1 ){
			 		stop("Item probability prior improperly specified. For varying numbers of categories, use alpha to specify the prior.")
			 	}else stop("Item probability prior improperly specified.")
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
				stop("alpha provided is not of a compatible length. Please check and rerun.")
			}
		}
	}
	
	if( is.character(start.vals[1]) ) 
	{
		if( start.vals[1] == "single" ) init <- 0 
		if( start.vals[1] == "across" ) init <- 1
		if( start.vals[1] == "prior" ) init <- 2
	}else{
		init <- 0
	}
	

	# thinning parameter
	
	Thinby <- floor(1/thin)

	#if((iter-burn.in)%%thin != 0)
	#	stop("\t Please thin by an amount that divides perfectly into iter - burn.in. ")

	stored <- floor( iter / Thinby )

	memberships <- numeric(stored*N)
	weights <- numeric(stored*G)
	variable.probs <- numeric(stored*G*sum(ncat))
	log.post <- numeric( stored )
	
	if( is.null(model.indicator) ) model.indicator <- rep(1,M) #modify this later to allow for flexible model specification

	hparam <- c( delta[1], beta[1] )
	
	# initialization is done randomly to groups
	
	w <- .C(	"BLCA_GIBBS_SAMPLER",													as.integer(X),
				as.integer(N),																as.integer(M),
				as.integer(ncat),															as.double(hparam),
				as.integer(prior.init.type),
				as.double(delta),															as.double(gamma),
				as.integer(init),
				as.integer(G),															 	as.integer(iter),											
				as.integer(burn.in),
				as.integer(Thinby),														memberships = as.integer(memberships),					
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
			gap <- stored * G * sum( ncat[1:(j-1)] )
		}
		#variable probabilities stacked by group and then iteration
		if( model.indicator[j] )
		{
			var.probs.l[[l]] <- matrix( w$variable.probs[(gap+1):(gap + stored*G*ncat[j])] , nrow = stored * G, ncol=ncat[j], byrow=TRUE )
			l <- l+1
		}
	}	
	
	weights.mat <- matrix( w$weights , nrow=stored, ncol=G, byrow=TRUE )

	if( relabel ) relabelled <- undo.label.switching( membership.mat, rep(G, stored) )
	
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
				var.probs.l[[l]][ it + 1:G , ] = var.probs.l[[l]][ it + relabelled$permutation[ k, 1:G ] , ]
				l <- l+1
			}
		}
	}
	
	v.probs <- list( mean=list(), sd=list() )
	tt <- seq( 0, stored-1, by=G )
	l <- 1
	for( j in 1:M )
	{
		if( model.indicator[j] )
		{
			v.probs$mean[[l]] <- list()
			v.probs$sd[[l]] <- list()
			v.probs$mean[[l]]  <- matrix( nrow = G, ncol = ncat[j] ) 
			v.probs$sd[[l]] <- matrix( nrow = G, ncol = ncat[j] )
			for( g in 1:G )
			{
				b <- var.probs.l[[l]][ tt+g , ]
				v.probs$mean[[l]][g,] <- apply( b, 2, mean )
				v.probs$sd[[l]][g,] <- apply( b, 2, sd )
			}
			rownames( v.probs$mean[[l]] ) <- paste( "Group", 1:G )
			colnames( v.probs$mean[[l]] ) <- paste("Cat",0:(ncat[j]-1) )
			rownames( v.probs$sd[[l]] ) <- paste( "Group", 1:G )
			colnames( v.probs$sd[[l]] ) <- paste("Cat",0:(ncat[j]-1) )
			l <- l+1
		}
	}
	
	if( is.null(colnames(X)) )
	{
		names( v.probs$mean ) <- names( v.probs$sd ) <- paste("Variable", which( model.indicator==1 ) )
	}else{
		names( v.probs$mean ) <- names( v.probs$sd ) <- colnames(X)[  which( model.indicator == 1 ) ]
	}
	
	#compile the list to return

	x <- list()
	
	x$call <- match.call()
	
	x$classprob <- apply( weights.mat, 2, mean )	
	x$itemprob <- v.probs$mean
	
	x$classprob.sd <- apply( weights.mat, 2, sd )
	x$itemprob.sd <- v.probs$sd
	
	# look after itemprob.tidy here
	vec.itemprobs <- unlist(x$itemprob)
	vec.itemprobs.sd <- unlist(x$itemprob.sd)
	
	itemprobs.group.ind <- rep(paste("Group", 1:G), times = sum(ncat[ which(model.indicator == 1) ]))
	itemprobs.var.ind <- rep(names(x$itemprob), times = G * ncat[ which(model.indicator == 1) ])
	itemprobs.cat.ind <- paste("Cat", rep(as.numeric(unlist(apply(t(ncat[ which(model.indicator == 1) ]), 2, function(x) 0:(x-1)))), each = G))
	
	x$samples <- list()
	
	x$samples$logpost <- w$log.posterior
	x$samples$Giter <- rep( G, stored )
	x$samples$itemprob <- var.probs.l
	x$samples$classprob <- t( weights.mat )
	
	if(any(ncat > 2)){
		x$itemprob.tidy <- data.frame(itemprob = vec.itemprobs, group = itemprobs.group.ind, variable = itemprobs.var.ind, category = itemprobs.cat.ind)
	} else { 
	  # rearrange for backwards compatibility with plotting functions
	  x$itemprob <- matrix(vec.itemprobs[itemprobs.cat.ind == "Cat 1"], nrow = G, ncol = sum(model.indicator), dimnames = list(paste("Group", 1:G), names(x$itemprob)))
	  x$itemprob.sd <- matrix(vec.itemprobs.sd[itemprobs.cat.ind == "Cat 1"], nrow = G, ncol = sum(model.indicator), dimnames = list(paste("Group", 1:G), names(x$itemprob)))
	  M.in <- sum(model.indicator)
	  arr <- array( dim=c(stored,G,M.in) )
	  for( m in 1:M.in ) 
	  {
	  	arr[,,m] <- matrix( var.probs.l[[m]][,1], nrow=stored , ncol=G , byrow = TRUE )
	  }
	  dimnames(arr)[[3]] <- names(v.probs$mean)
	  x$samples$itemprob <- arr
	  }
	  
	#x$itemprob.tidy <- itemprob.tidy		
	
	Dbar<- mean(w$log.posterior)
	S2<- var(w$log.posterior)
	
	x$DIC<- 2*(2*Dbar - x$samples$logpost)
	# need to check definition of BICM (number of samples)
	x$BICM<- 2*(Dbar - S2*(log(N)-1))
	x$AICM<- 2*(Dbar - S2)

	if(relabel)
	{
		x$labelstore <- relabelled$relab
		x$labelstore.permutation <- relabelled$permutation
		x$G.Z <- relabelled$components
		x$Z <- relabelled$label.probs
	}
			
	#inputs
	x$iter <- iter
	x$burn.in <- burn.in
	x$thin <- thin
	x$relabel <- relabel
	x$prior <- list( alpha=alpha, beta=beta, delta=delta )
	x$ncat <- ncat
	x$model.indicator <- model.indicator
	
	t2 = proc.time()
	
	ttime = t2-t1
  
  	x$time = ttime[3]	
  	
  	if( any(ncat > 2) ) class(x) <- c("blca.gibbs", "blca.multicat", "blca") else class(x) <- c("blca.gibbs", "blca")
	
	return(x)
}
