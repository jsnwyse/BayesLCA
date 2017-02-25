blca.vb <-
function( X, G, alpha=1, beta=1, delta=1, num.categories=NULL, model.indicator = NULL, iter=500, restarts=5, verbose=TRUE, conv=1e-6, small=1e-100 )
{
	
	X <- as.matrix(X)

	N<-nrow(X) 
	M<-ncol(X) 
	
	if( is.null(model.indicator) )
		model.indicator <- rep(1,M)
	
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

	if(is.numeric(restarts)){
	  if(length(restarts)>1){
	    restarts<- restarts[1]
	    warning("restarts improperly specified - first value will be used, other values will be ignored")
	    }# else {stop("restarts improperly specified. Must be an integer of length 1.")}
	} else {stop("restarts improperly specified. Must be an integer of length 1.")}
	
	#storage 
	weights <- numeric(G)
	sd.weights <- numeric(G)
	variable.probs <- numeric(G*sum(num.categories))
	sd.variable.probs <- numeric(G*sum(num.categories))
	group.probs <- numeric(N*G)
	lb <- numeric( iter )	
	iters <- 0
	converged <- 0
	logpost <- 0

	hparam <- c( delta, beta )
	
	lb.max <- -Inf
	w.max <- list()
	
	for( r in 1:restarts )
	{	
		# call the VB algorithm
		
		w <- .C( 		"BLCA_VB_FIT", 										as.integer(X), 			
							as.integer(N),											as.integer(M), 
							as.integer(num.categories),	 					as.double(hparam),
							as.integer(G),											as.integer(iter),
							iters = as.integer(iter),							group.probs = as.double(group.probs),
							weights = as.double(weights),						se.weights = as.double( sd.weights ),
							variable.probs = as.double(variable.probs),	se.variable.probs = as.double( sd.variable.probs ),
							as.integer(model.indicator),						lb = as.double(lb),
							as.double(conv),										converged = as.integer(converged),
							logpost = as.double(logpost),
							PACKAGE = "BayesLCA" )
	
		lb.this <- w$lb[ w$iters ]
		
		#store the run that gives the highest value  of the log posterior
		#	for the runs that have converged
		new.max <- FALSE
		if( lb.this > lb.max && w$converged == TRUE )
		{
			w.max <- w 
			lb.max <- lb.this
			new.max <- TRUE
		}
		
		if( verbose ) 
		{
			if( new.max && r>1 )
			{
				cat( "\nNew maximum found... Restart number ",r,", logpost = ", round(w$logpost,2),"...", sep = "" )
			}else{ 
				cat( "\nRestart number ",r,", logpost = ", round(w$logpost,2),"...", sep = "" )
			}
		}
	
	}	
	cat("\n")
	
	x <- list()
	x$call <- match.call()
	
	x$converged <- w$converged
	x$itersconv <- w$iters
	
	x$classprob <- w.max$weights
	x$classprob.se <- w.max$se.weights
	x$Z <- matrix( w.max$group.probs, nrow=N, ncol=G )
	colnames(x$Z) <- paste( "Group", 1:G )
	
	var.probs.l <- list()
	se.var.probs.l <- list()
	
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
			se.var.probs.l[[l]] <- matrix( w$se.variable.probs[(gap+1):(gap + G*num.categories[j])] , nrow =  G, ncol=num.categories[j], byrow=TRUE )
			rownames( var.probs.l[[l]] ) <- rownames( se.var.probs.l[[l]] ) <-  paste( "Group", 1:G )
			colnames( var.probs.l[[l]] ) <- rownames( se.var.probs.l[[l]] ) <-  paste("Cat",0:(num.categories[j]-1) )
			l <- l+1
		}
	}
	
	if( is.null(colnames(X)) )
	{
		names( var.probs.l ) <- names( se.var.probs.l ) <- paste("Variable", which( model.indicator==1 ) )
	}else{
		names( var.probs.l ) <- names( se.var.probs.l ) <- colnames(X)[  which( model.indicator == 1 ) ]
	}
	
	x$itemprob <- var.probs.l	
	x$itemprob.se <- se.var.probs.l
	
	vec.itemprobs <- unlist( x$itemprob )
	
	itemprobs.group.ind <- rep( paste("Group", 1:G), times = sum(num.categories) )
	
	itemprovs.var.ind <- rep( names(x$itemprob), times = G*num.categories )
	
	itemprobs.cat.ind <- paste("Cat", rep(as.numeric(unlist(apply(t(num.categories), 2, function(x) 0:(x-1)))), each = G))
	
	if(any(num.categories > 2)){
		itemprob.tidy <- data.frame(itemprob = vec.itemprobs, group = 			itemprobs.group.ind, variable = itemprobs.var.ind, category = itemprobs.cat.ind)
	} else { 
	  	itemprob.tidy <- matrix(vec.itemprobs[itemprobs.cat.ind == "Cat 1"], nrow = G, ncol = M, dimnames = list(paste("Group", 1:G), names(x$itemprob)))
	  }
	
	x$itemprob.tidy <- itemprob.tidy	
	
	x$lbstore <- w$lb[ 1:w$iters ]
	x$lb <- lb.max
  
	#x$BIC<- 2*likl-(G*M + G-1)*log(N1)
	#x$AIC<- 2*likl - 2*(G*M + G-1)
	#x$iter<- length(rstore$llstore)
	
	x$prior<-NULL
	x$prior$alpha<- alpha
	x$prior$beta<- beta
	x$prior$delta<- delta

	if( any(num.categories>2) ) class(x)<-c("blca.vb", "blca.multicat" "blca") else class(x) <- c("blca.vb", "blca" )

	return(x)
}
