blca.collapsed <- function( X,  G, ncat=NULL , alpha=1, beta=1, delta=1, start.vals=c("single"), counts.n=NULL, iter=5000, burn.in=1000, thin=1, G.sel=FALSE, var.sel=FALSE, post.hoc.run=TRUE, n.gibbs=nrow(X), only.gibbs=TRUE, G.max=30, G.prior=dpois(1:G.max, lambda=1), model.indicator=NULL, prob.inc=.5, hprior.model=FALSE, relabel=TRUE, verbose=TRUE, verbose.update=1000 )
{

	t1 <- proc.time()
	
	# convert X into a numeric matrix and check inputs
	D <- blca.check.data( X, counts.n, ncat )
	
	X <- D$X
	ncat <- D$ncat

	N <- nrow(X)
	M <- ncol(X)
	
	## safety checks ##
	
	#if( length(alpha) > 1 ) stop("alpha value must be a positive scalar when using collapsed sampler \n")
	#if( length(beta) > 1 ) stop("beta value must be a positive scalar when using collapsed sampler \n")
	#if( length(delta) > 1 ) stop("delta value must be a positive scalar when using collapsed sampler \n")
	
	# priors - there are alot of checks here to ensure things run smoothly
	
	prior.init.type <- 1
	
	if( G.sel == FALSE ) G.max <- G
	
	# delta is the prior on the group weights
	if( length(delta) == 1 ) delta<- rep(delta,G.max)
	if( length(delta) != G.max ) stop("delta prior is wrong length (i.e., != 1 or G.max)")
	if( length(delta) == G.max && any( delta != delta[1] ) && G.sel )
		stop("Unsymmetric priors on group weights are not allowed when G.sel is set to TRUE. Rerun with delta set to a scalar.")
	
	# alpha either acts only as the first category prior for binary, or the entire prior for multicategory case
	if( !is.matrix(alpha) && !any( ncat > 2 ) ){
	if( any(length(alpha)==c(1,G.max)) ){
		if( length(alpha) == G.max && G.sel && any( alpha != alpha[1] ) ) stop("The prior setting for alpha is required to be the same over all groups when G.sel is set to TRUE. Rerun with alpha set to a scalar.") 
		alpha <- matrix(alpha,G.max,M)
		if( G.sel ) alpha <- matrix( alpha[1], G.max, M ) 
	}else{
		if(length(alpha)==M){
			# this is ok as there is a different alpha for each variable but the same across groups
			alpha <- matrix( alpha, G.max, M, byrow=TRUE)
			} else {
			 	if( length(alpha) == G.max*sum(ncat) && beta==1 ){
			 		if( !any( ncat > 2 ) ) warning("Using only entries in alpha to assign prior for sampling")
			 		if( G.sel  ) 
			 		{
			 			# here we need a thorough check to ensure that the specification for the prior weight on a category does not change over groups
			 			u <- rep(FALSE, M )
			 			for( m in 1:M )
			 			{
			 				if( m == 1 ) idx0 <- 0 else idx0 <- G.max * sum( ncat[1:(m-1)] )
			 				v <- alpha[ (idx0 + 1) : (idx0 + G.max*ncat[m]) ]
			 				if( any( v != v[1] ) ) u[m] <- TRUE
			 			}
			 			if( any(u) ) stop("The prior setting for alpha is required to be the same over all groups when G.sel is set to TRUE. Rerun with alpha set to a scalar. Rerun") 
			 		}
			 		# pass alpha directly, this looks ok
			 		prior.init.type <- 2
			 	}else stop("alpha: Item probability prior improperly specified.")
			}
		}
	} 
	
	if( is.matrix(beta) && any( ncat > 2 ) ) stop("Item probability prior improperly specified. Please use alpha to specify the prior." )
	
	# beta either acts only as the first category prior for binary 
	if( !is.matrix(beta)  && !any( ncat > 2 ) ){
	if( any(length(beta)==c(1,G.max)) ){
		if( length(beta) == G.max && G.sel && any(  beta != beta[1] ) ) 
			stop("Unsymmetric priors with G selection is not available for collapsed sampler... rerun with beta set to a scalar")
		beta <- matrix(beta,G.max,M)
		if( G.sel ) beta <- matrix( beta[1], G.max, M )
	}else{
		if(length(beta)==M){
			if( G.sel && any( beta != beta[1] ) ) 
				stop("Unsymmetric priors with G selection is not available for collapsed sampler... rerun with beta set to a scalar") 
			beta <- matrix(beta,G.max,M, byrow=TRUE)
			} else {
			 	if( length(beta) == G.max*sum(ncat) && alpha==1 ){
			 		stop("Item probability prior improperly specified. For varying numbers of categories, use alpha to specify the prior.")
			 	}else stop("beta: Item probability prior improperly specified.")
			}
		}
	} 
	
	if( !any( ncat > 2 ) )
	{
		if( !all(dim(alpha) == c(G.max,M)) ) stop("alpha has not been specified with the correct dimensions")
		if( !all( dim(beta) == c(G.max,M) ) ) stop("beta has not been specified with the correct dimensions")
		# if the priors are unsymmetric in all instances, disable moves 1, 2 and 3
		if( any( alpha != alpha[1] ) || any( beta != beta[1] ) )
		{
			cat("Note: disabling label update moves other than Gibbs sampling due to identifiability of priors.")
			only.gibbs <- TRUE
		}
		# now restack the alpha and beta matrices into compatible format
		gamma <- matrix( nrow=2*G.max, ncol=M )	
		for( k in 1:G.max ) 
		{
			gamma[ 2*(k-1) + (1:2) ,  ] <- rbind( alpha[k,], beta[k,] )
		}
		prior.init.type <- 2
	}else{
		if( length(alpha) == 1 ){
			prior.init.type <- 1
			gamma <- rep( alpha, G.max*sum(ncat) )
		}else{
			if( length(alpha) == sum(ncat) ) 
			{
				gamma <- rep( alpha, G.max )
			}else if( length(alpha) == G.max*sum(ncat) ){
				gamma <- alpha
			}else{
				stop("alpha provided is not of a compatible length. Please check and rerun.")
			}
		}
	}

	
	if( G.sel & G.max < G ){
		stop("The maximum number of groups cannot be less than the initial number of groups. Please rerun with G.max (default 30) appropriately specified.")
	}
	
	# thinning parameter
	if( thin > 1 || thin == 0 ) stop("Argument thin gives a thinning rate and should be less than or equal to 1.")
	
	Thinby <- floor(1/thin)

	#if((iter-burn.in)%%thin != 0)
	#	stop("\t Please thin by an amount that divides perfectly into iter - burn.in. ")

	stored <- floor( iter / Thinby )	

	memberships <- numeric(stored*N)
	variable.inclusion.indicator <- numeric(stored*M)
	num.groups <- numeric(stored)
	log.post <- numeric(stored)
	prior.include <- numeric(stored)
	
	if( is.null(model.indicator) ) model.indicator <- rep(1,M) 
	
	hparam <- c( delta[1], beta[1] )
	accrts <- numeric(6)
	
	# working on passing the priors--
	
	w <- .C(	"BLCA_VS",																	as.integer(X),
				as.integer(N),																as.integer(M),
				as.integer(ncat),															as.double(hparam),
				as.integer(prior.init.type),											as.double(delta),
				as.double(gamma),
				as.integer(!G.sel),														as.integer(only.gibbs),
				as.integer(G),																as.integer(G.max),
				as.integer(iter),															as.integer(burn.in),
				as.integer(Thinby),														as.integer(n.gibbs),
				memberships = as.integer(memberships),								ngroups = as.integer(num.groups),
				as.double(G.prior),														as.integer(var.sel),
				variable.inclusion.indicator = as.integer(variable.inclusion.indicator),	as.double(prob.inc),
				log.posterior = as.double(log.post),								as.integer(hprior.model),
				prior.include = as.double(prior.include),							as.integer(model.indicator),
				as.integer(verbose),														as.integer(verbose.update),
				accrts=as.double(accrts),
				PACKAGE = "BayesLCA" )
				
	
	if( verbose ) cat("\nFinished sampling...\n")

	membership.mat <- matrix( w$memberships, nrow = stored, ncol=N, byrow=FALSE ) + 1
	vindicator.mat <- matrix( w$variable.inclusion.indicator, nrow=stored, ncol=M, byrow=FALSE )
	
	
	if( relabel ) 
	{
		if( G.sel ) grps <- w$ngroups else grps <- rep( G , stored )
		relabelled <- undo.label.switching( membership.mat, grps )
	}
	
	x = list()
	
	x$call <- match.call()
	
	x$classprob <- NULL
	x$classprob.sd <- NULL
	x$itemprob <- NULL
	x$itemprob.sd <- NULL
	
	x$samples <- list( logpost=w$log.posterior, G=w$ngroups )
	acc <- w$accrts
	x$accrates <- list()
	if( !only.gibbs ) x$accrates$labels <- list( m1 = acc[1] , m2 = acc[2] , m3 = acc[3] )
	if( G.sel ) x$accrates$G <- list( eject = acc[4], absorb = acc[5] )
	if( var.sel ) x$accrates$var.sel <- acc[6]
	
	if( var.sel )
	{
		x$samples$var.ind <- vindicator.mat
		if(!is.null(colnames(X))) colnames(x$samples$var.ind) <- colnames(X)
	}

	if( hprior.model ) x$samples$prob.inc <- w$prior.include
	
	if( relabel )
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
	x$G.sel <- G.sel
	#x$just.gibbs.updates <- just.gibbs.updates
	x$n.gibbs <- n.gibbs
	x$var.sel <- var.sel
	x$prob.inc <- prob.inc
	x$hprior.model <- hprior.model
	x$relabel <- relabel
	x$prior <- list( alpha=alpha, beta=beta, delta=delta, prob.inc=prob.inc, G.prior=G.prior, G.max=G.max )
	x$ncat <- ncat
	
	t2 = proc.time()
	
	ttime = t2-t1
  
  	x$time = ttime[3]
	
	#if( verbose ) cat("\nTotal time taken (approx): ",ttime[3]," seconds")
	
	#if(relabel && matchClasses(t(Z)%*%Z1, method="exact", verbose=FALSE)) warning("Label-switching (provisionally) corrected for - proceed with caution")
	#if(relabel && label.swap) warning("Label-switching (provisionally) corrected for - diagnostic plots are recommended. Use '?plot.blca' for details.")
	#if(!relabel && label.swap) warning("Label-switching may have occurred - diagnostic plots are recommended. Use '?plot.blca' for details.")
	
	if( post.hoc.run ){
	  
	  Gstar <- sort(unique(x$samples$G ), decreasing = FALSE)[which.max(table(x$samples$G))]
	  if( var.sel )
	  {
	 	 Mstar <- which(colMeans(x$samples$var.ind)>0.5)
   	}else{
	  	 Mstar <- 1:M
	  }
	  gamma.star <- rep(0,M)
	  gamma.star[Mstar] <- 1
	  #cat("\n Gamma star is ", gamma.star )
	    #ph.est <- blca.collapsed.post.hoc.estimates(x, Mstar, Gstar)
	  ph.est <- blca.gibbs( X, G=Gstar, model.indicator=gamma.star, ncat=ncat, verbose=FALSE )
	    
	  x$classprob <- ph.est$classprob
	  x$classprob.sd <- ph.est$classprob.sd
	  x$itemprob <- ph.est$itemprob #sapply(ph.est$itemprob, function(t1) t1[, 2])
	  x$itemprob.sd <- ph.est$itemprob.sd 
	  x$itemprob.tidy <- ph.est$itemprob.tidy#sapply(ph.est$itemprob.sd, function(t1) t1[, 2])

	  ## Name item probabilities
	  if( FALSE )
	  {
	  if(is.matrix(x$itemprob)){
	    if(is.null(colnames(X))) {
	      colnames( x$itemprob)<- paste("Col", Mstar)
	      } else{
	        colnames( x$itemprob)<- colnames(X)[Mstar]
	        } ## else is.null(colnames(X))
	    } else{
	      if(is.null(colnames(X))){
	        names(x$itemprob)<- paste("Col", Mstar)
	        } else{
	          names( x$itemprob)<- colnames(X)[Mstar]
	          } ## else is.null(colnames(X))
	      } ## else is.matrix(x$itemprob)
	    }
	} ## if post.hoc.est
	
	
	if( any(ncat > 2 ) ) class(x) <- c("blca.collapsed", "blca.multicat", "blca") else class(x)<-c("blca.collapsed", "blca")
	
	return(x)
	
}
