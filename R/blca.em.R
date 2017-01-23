blca.em <-
function( X, G, alpha=1, beta=1, delta=1, num.categories=NULL, iter=500, restarts=5, verbose=TRUE, sd=FALSE, se=sd, conv=1e-6, small=1e-100, MAP=FALSE )
{
	
	X <- as.matrix(X)

	N<-nrow(X) 
	M<-ncol(X) 
	
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
	if( sum(t) > 0 ) stop("\t please recode categories from 0, ..., ncat-1  to use blca.em\n")
	t <- apply( X, 2, max )
	if( sum( t + 1 - num.categories ) > 0 ) stop("\n number of categories in X exceeds num.categories please recode categories from 0, ..., ncat-1\n")
		

	## safety checks ##
	
	if(length(num.categories)!= M){
		stop("\t The length of num.categories must be the same as the number of records per observation\n")
	}	
	
	if( length(alpha) > 1 ) stop("\t alpha value must be a positive scalar \n")
	if( length(beta) > 1 ) stop("\t beta value must be a positive scalar \n")
	if( length(delta) > 1 ) stop("\t delta value must be a positive scalar \n")

	if(2^M <= (M+1)*G) warning(paste("Model may be improperly specified. Maximum number of classes that should be run is", floor(2^M/(M+1)), "."))
	
	counter<-0
	llstore<-0
	llcheck<- -Inf

	if(is.numeric(restarts)){
	  if(length(restarts)>1){
	    restarts<- restarts[1]
	    warning("restarts improperly specified - first value will be used, other values will be ignored")
	    }# else {stop("restarts improperly specified. Must be an integer of length 1.")}
	} else {stop("restarts improperly specified. Must be an integer of length 1.")}

	multistart.lp.store<- rep(0, restarts)
	if(sd!=se) se<- sd
	
	#storage 
	weights <- numeric(G)
	variable.probs <- numeric(G*sum(num.categories))
	group.probs <- numeric(N*G)
	log.post <- numeric( iter )	
	iters <- 0
	converged <- 0

	model.indicator <- rep(1,M) #modify this later to allow for flexible model specification

	hparam <- c( delta, beta )
	
	log.post.max <- -Inf
	w.max <- list()
	
	for( r in 1:restarts )
	{	
		# call the EM algorithm
		
		w <- .C( 		"BLCA_EM_FIT", 							as.integer(X), 			
							as.integer(N),								as.integer(M), 
							as.integer(num.categories),	 		as.double(hparam),
							as.integer(G),								as.integer(iter),
							iters = as.integer(iter),				group.probs = as.double(group.probs),
							weights = as.double(weights),			variable.probs = as.double(variable.probs),
							as.integer(model.indicator),			log.post = as.double(log.post),
							as.integer(MAP),							as.double(conv),
							converged = as.integer(converged),
							PACKAGE = "BayesLCA" )
	
		log.post.this <- w$log.post[ w$iters ]
		
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
	
	}	
	cat("\n")
	
	x <- list()
	x$call <- match.call()
	
	x$converged <- w$converged
	x$itersconv <- w$iters
	
	x$classprob <- w.max$weights
	x$Z <- matrix( w.max$group.probs, nrow=N, ncol=G )
	colnames(x$Z) <- paste( "Group", 1:G )
	
	var.probs.l <- list()
	
	l <- 1
	for( j in 1:M )
	{
		if(j == 1){
			gap <- 0
		}else{
			gap <- G * sum( num.categories[1:(j-1)] )
		}
		#variable probabilities stacked by group and then iteration
		if( model.indicator[j] )
		{
			var.probs.l[[l]] <- matrix( w$variable.probs[(gap+1):(gap + G*num.categories[j])] , nrow =  G, ncol=num.categories[j], byrow=TRUE )
			rownames( var.probs.l[[l]] ) <-  paste( "Group", 1:G )
			colnames( var.probs.l[[l]] ) <- paste("Cat",0:(num.categories[j]-1) )
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
	
	itemprobs.group.ind <- rep(paste("Group", 1:G), times = sum(num.categories))
	
	itemprobs.var.ind <- rep(names(x$itemprob), times = G * num.categories)
	
	itemprobs.cat.ind <- paste("Cat", rep(as.numeric(unlist(apply(t(num.categories), 2, function(x) 0:(x-1)))), each = G))
	
	if(any(num.categories > 2)){
	itemprob.tidy <- data.frame(itemprob = vec.itemprobs, group = itemprobs.group.ind, variable = itemprobs.var.ind, category = itemprobs.cat.ind)
	} else { 
	  itemprob.tidy <- matrix(vec.itemprobs[itemprobs.cat.ind == "Cat 1"], nrow = G, ncol = M, dimnames = list(paste("Group", 1:G), names(x$itemprob)))
	  }
	x$itemprob.tidy <- itemprob.tidy	
	
	x$poststore <- w$log.post[ 1:w$iters ]
	x$logpost <- log.post.max
  
		
	#x$BIC<- 2*likl-(G*M + G-1)*log(N1)
	#x$AIC<- 2*likl - 2*(G*M + G-1)
	#x$iter<- length(rstore$llstore)
	
	x$prior<-NULL
	x$prior$alpha<- alpha
	x$prior$beta<- beta
	x$prior$delta<- delta

		if(FALSE){
		if(sd){
#			if(any(x$itemprob==0)){ warning("some item probability estimates are exactly zero. standard errors in this case are undefined.")}
#			if(any(x$classprob==0)){ warning("some class probability estimates are exactly zero. standard errors in this case are undefined.")}
			s.e.<- blca.em.sd(x,X,counts.n)
			x$itemprob.sd<- x$itemprob.se<- s.e.$itemprob
			x$classprob.sd<- x$classprob.se<- s.e.$classprob
			convergence<- s.e.$convergence
		} else convergence<- 0
		
		if(counter>iter){ 
		  convergence<- 3 
		  warning("maximum iteration reached - algorithm not deemed to have converged. rerunning the function with 'iter' set to a higher value is recommended.")
		} #else{ convergence<- s.e.$convergence}
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
	if(any(num.categories > 2)) class(x)<- c("blca.em", "blca.multicat","blca") else class(x)<-c("blca.em", "blca")

		return(x)
}
