blca.em <-
function( X, G, ncat=NULL, alpha=1, beta=1, delta=1, start.vals = c("single","across"), counts.n=NULL, model.indicator=NULL, iter=1000, restarts=5, verbose=TRUE, sd=FALSE, se=sd, conv=1e-6, small=1e-100, MAP=TRUE, pars.init=NULL, for.boot=FALSE )
{
	
	# convert X into a numeric matrix and check inputs
	D <- blca.check.data( X, counts.n, ncat )
	
	X <- D$X
	ncat <- D$ncat

	N<-nrow(X) 
	M<-ncol(X)
	
	if( is.null(model.indicator) )
	{
		model.indicator <- rep(1,M)
	}else if( length(model.indicator) != M ){
		stop("model.indicator must have length ncol(X)")
	}
	
	prior.init.type <- 1
	
	# delta is the prior on the group weights
	if( length(delta) == 1 ) delta<- rep(delta,G)
	if( length(delta) != G ) stop("delta prior is wrong length (i.e., != 1 or G)")
	
	# alpha either acts only as the first category prior for binary, or the entire prior for multicategory case
	if( !is.matrix(alpha) && !any( ncat > 2 ) ){
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
	
	M.in <- sum( model.indicator ) 
	if( prod( ncat[model.indicator==1] ) <= (M.in+1)*G) warning(paste("Model may be improperly specified. Maximum number of classes that should be run is ", floor(prod(ncat[model.indicator==1])/(M.in+1)), "."))
	
	counter<- 0
	llstore<- 0
	llcheck<- -Inf

	if(is.numeric(restarts)){
	  if(length(restarts)>1){
	    restarts<- restarts[1]
	    warning("restarts improperly specified - first value will be used, other values will be ignored")
	    }# else {stop("restarts improperly specified. Must be an integer of length 1.")}
	} else {stop("restarts improperly specified. Must be an integer of length 1.")}

	multistart.lp.store<- rep(0, restarts)
	if( sd != se ) se<- sd
	
	#storage 
	weights <- numeric(G)
	variable.probs <- numeric(G*sum(ncat))
	group.probs <- numeric(N*G)
	log.post <- numeric( iter )	
	iters <- 0
	eps <- 0
	converged <- 0

	hparam <- c( delta[1], beta[1] )
	
	log.object.max <- 0.
	
	log.post.max <- -Inf
	w.max <- NULL
	
	init.vals <- numeric(N*G)
	if( is.character(start.vals[1]) ) 
	{
		if( start.vals[1] == "single" ) init <- 0 else init <- 1
	}else{
		init <- 2
		if(is.matrix(start.vals) & all(dim(as.matrix(start.vals)) == c(N,G))){
		  init.vals<- as.vector(start.vals) 
		  }else{
		  if(is.numeric(start.vals) & length(as.numeric(start.vals))==N)
		  { 
		    init.vals<- as.vector(unMAP(start.vals))
		    }else{ 
		    stop("start.vals improperly specified. See help files for more details.")
		   }
		}
	}
	
	lpstarts <- numeric(restarts)
	
	for( r in 1:restarts )
	{	
		# call the EM algorithm
		
		w <- .C( 		"BLCA_EM_FIT", 							as.integer(X), 			
							as.integer(N),								as.integer(M), 
							as.integer(ncat),	 						as.double(hparam),
							as.integer(prior.init.type), 			as.double(delta),
							as.double(gamma),
							as.integer(G),								as.integer(init),
							as.double(init.vals), 					as.integer(iter),
							iters = as.integer(iter),				group.probs = as.double(group.probs),
							weights = as.double(weights),			variable.probs = as.double(variable.probs),
							as.integer(model.indicator),			log.post = as.double(log.post),
							as.integer(MAP),							as.double(conv),
							eps = as.double(eps),
							converged = as.integer(converged),	log.object.max = as.double(log.object.max),
							PACKAGE = "BayesLCA" )
	
		log.post.this <- w$log.post[ w$iters ]
		
		lpstarts[r] <- log.post.this
		
		#store the run that gives the highest value  of the log posterior
		#	for the runs that have converged
		
		new.max <- FALSE
		if( log.post.this > log.post.max && w$converged == TRUE )
		{
			w.max <- w 
			log.post.max <- log.post.this
			new.max <- TRUE
		}
		
		if( verbose ) 
		{
			if( new.max && r>1 )
			{
				cat( "\nNew maximum found... Restart number ",r,", logpost = ", round(log.post.this,2),"...", sep = "" )
			}else{ 
				cat( "\nRestart number ",r,", logpost = ", round(log.post.this,2),"...", sep = "" )
			}
		}
		
		if( w$converged == FALSE ) warning("\nRestart number ", r," failed to converge. Rerun with a higher iter value.")
	
	}	
	if( verbose ) cat("\n")
	
	#print(w.max)
	
	x <- list()
	x$call <- match.call()
	
	
	var.probs.l <- list()
	
	l <- 1
	for( j in 1:M )
	{
		if(j == 1){
			gap <- 0
		}else{
			gap <- G * sum( ncat[1:(j-1)] )
		}
		#variable probabilities stacked by group and then iteration
		if( model.indicator[j] )
		{
			var.probs.l[[l]] <- matrix( w.max$variable.probs[(gap+1):(gap + G*ncat[j])] , nrow =  G, ncol=ncat[j], byrow=TRUE )
			rownames( var.probs.l[[l]] ) <-  paste( "Group", 1:G )
			colnames( var.probs.l[[l]] ) <- paste("Cat",0:(ncat[j]-1) )
			l <- l+1
		}
	}
	
	if( is.null(colnames(X)) )
	{
		names( var.probs.l ) <- paste("Variable", which( model.indicator==1 ) )
	}else{
		names( var.probs.l ) <- colnames(X)[  which( model.indicator == 1 ) ]
	}
	
	x$itemprob <- var.probs.l	
	
	vec.itemprobs <- unlist(x$itemprob)
	
	itemprobs.group.ind <- rep(paste("Group", 1:G), times = sum(ncat[ which(model.indicator == 1) ]))
	itemprobs.var.ind <- rep(names(x$itemprob), times = G * ncat[ which(model.indicator == 1) ])
	itemprobs.cat.ind <- paste("Cat", rep(as.numeric(unlist(apply(t(ncat[ which(model.indicator == 1) ]), 2, function(x) 0:(x-1)))), each = G))
	
	
	
	if(any(ncat > 2)){
		x$itemprob.tidy <- data.frame(itemprob = vec.itemprobs, group = itemprobs.group.ind, variable = itemprobs.var.ind, category = itemprobs.cat.ind)
	} else { 
	  x$itemprob <- matrix(vec.itemprobs[itemprobs.cat.ind == "Cat 1"], nrow = G, ncol = M, dimnames = list(paste("Group", 1:G), names(x$itemprob)))
	  }
	  
	#x$itemprob.tidy <- itemprob.tidy	
	
	x$classprob <- w.max$weights
	
	# include the itemprob.sd here?- this  is still left to do	
	
	
	
	x$Z <- matrix( w.max$group.probs, nrow=N, ncol=G )
	colnames(x$Z) <- paste( "Group", 1:G )	
	

	
	x$logpost <- log.post.max
  
  	likl <- w.max$log.object.max
  	
  	npars <- G * sum( model.indicator * ( ncat - 1 ) ) + G - 1
	
	x$BIC<- 2*likl- npars*log(N) ## this penalty will not be correct for multicat
	x$AIC<- 2*likl - 2*npars
	x$iter<- w.max$iters
	x$poststore <- w.max$log.post[ 1:w.max$iters ]
	x$converged <- w.max$converged
	x$eps <- w.max$eps
	x$conv <- conv
	x$lpstarts <- lpstarts
	
	
	if( MAP ) x$prior <- list( alpha=alpha, beta=beta, delta=delta )
	x$model.indicator <- model.indicator
	x$MAP <- MAP
	
	x$ncat <- ncat
	
	x$for.boot <- for.boot
	if( sd ) x$for.boot <- TRUE
	if( for.boot )
	{
		x$boot.init = list( weights=w.max$weights, var.probs=w.max$variable.probs )
	}
	
	if( sd )
	{
		if( verbose ) cat("\nRunning bootstrap sampler to get variance estimates...")
		y <- blca.boot( X, fit=x, iter=100, B=200, verbose=FALSE  )	
		x$itemprob.sd <-y$itemprob.sd
		x$classprob.sd <- y$classprob.sd
		if( verbose ) cat("\nBootstrap run complete...\n")
	}
	

	if(FALSE){
		# still need to do this
		if( sd && !any(ncat > 2) ){
#			if(any(x$itemprob==0)){ warning("some item probability estimates are exactly zero. standard errors in this case are undefined.")}
#			if(any(x$classprob==0)){ warning("some class probability estimates are exactly zero. standard errors in this case are undefined.")}
			s.e.<- blca.em.sd(x,X,counts.n)
			x$itemprob.sd<- x$itemprob.se<- s.e.$itemprob
			x$classprob.sd<- x$classprob.se<- s.e.$classprob
			convergence<- s.e.$convergence
		} else {
			if( any(ncat > 2 ) ) warning("Estimates of variablility only available for binary problems.")
			convergence<- 0
		}
		
		#if(counter>iter){ 
		#  convergence<- 3 
		#  warning("maximum iteration reached - algorithm not deemed to have converged. rerunning the function with 'iter' set to a higher value is recommended.")
		#} #else{ convergence<- s.e.$convergence}
		if(convergence==2) {warning("some point estimates likely converged at saddle-point. at least some points will not be at a local maximum. \n rerunning the function with a larger number of restarts is recommended.")}
		if(convergence==4){ warning("some point estimates located at boundary (i.e., are 1 or 0). posterior standard deviations will be 0 for these values.")}
		x$convergence<- convergence
		} else if(FALSE){ 
		  x$convergence<- 1
		  x$classprob.sd<- x$classprob.se<- 0
		  x$itemprob.sd<- x$itemprob.se<- sqrt( ((Thetat*N1 + alpha)*( (1 - Thetat)*N1 + beta))/( (N1 + alpha + beta + small)^2 * (N1 + alpha + beta + 1 + small) ) )
	}
		#x$small<- small
		#if((se==TRUE)&&(is.null(s.e.$classprob))) se<- FALSE
		#x$sd<- x$se<- se
	if(any(ncat > 2)) class(x)<- c("blca.em", "blca.multicat","blca") else class(x)<-c("blca.em", "blca")

		return(x)
}
