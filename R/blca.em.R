blca.em <-
function( X, G, formula = NULL, ncat=NULL, alpha=1, beta=1, delta=1, 
          start.vals = c("single","across"), counts.n=NULL, iter=2000, restarts=5, verbose=TRUE, sd=FALSE, 
          sd.method=c("delta","boot"), conv=1e-6, small=1e-10, MAP=TRUE, pars.init=NULL, for.boot=FALSE )
{
	# check if data is simulated 
  #if( class(X)[1] == "blca.rand" & !is.matrix(X) ) X <- X$X 
  
  args.passed <- as.list( environment() )
  #list of returns
  x <- list()
  x$call <- match.call()
  x$args <- args.passed

	# convert X into a numeric matrix for passage to C (index from 0) and check inputs
	D <- blca.prep.data( X, formula, counts.n, ncat )
	
	X <- D$X
	x$args$formula <- D$formula
	ncat <- D$ncat
	x$args$ncat <- ncat
	levnames <- D$levnames
	x$G <- G
	
	
	if( !is.null(D$missing.idx) )
	{
	  warning("Missing values encountered in X: rows with NA have been removed", call.=FALSE)
	  X <- na.omit(X)	  
	}

	N<-nrow(X) 
	M<-ncol(X)
	
	model.indicator <- rep(1,M)
	M.in <- sum( model.indicator ) 
	if( M.in == 0 ) stop("there are no variables included: check model.indicator")
	if( prod( ncat[model.indicator==1] ) <= (M.in+1)*G) stop(paste("maximum numer of classes that should be run for this data is ", floor(prod(ncat[model.indicator==1])/(M.in+1)) ))

	out.prior <- blca.check.prior( alpha, beta, delta, G, M, ncat )
	prior.init.type <- out.prior$prior.init.type
	gamma <- out.prior$gamma
	delta <- out.prior$delta
	
	if(is.numeric(restarts)){
	  if(length(restarts)>1){
	    restarts<- restarts[1]
	    warning("Restarts improperly specified: first value will be used, other values will be ignored", call.=FALSE)
	    }# else {stop("restarts improperly specified. Must be an integer of length 1.")}
	} else {stop("restarts improperly specified: must be an integer valued")}

	multistart.lp.store<- rep(0, restarts)

	#storage 
	weights <- numeric(G)
	variable.probs <- numeric(G*sum(ncat))
	group.probs <- numeric(N*G)
	log.post <- numeric( iter )	
	iters <- 0
	eps <- 0
	converged <- 0

	hparam <- c( delta[1], beta[1] )
	
	log.object.max <- 0.0
	log.post.max <- -Inf
	w.max <- NULL
	
	plus.plus <- FALSE
	init.vals <- numeric(N*G)
	if( is.character(start.vals[1]) ) 
	{
	  init <- match( start.vals[1], c("single","across")) - 1
	}else{
		init <- 2
		if(is.matrix(start.vals) & all(dim(as.matrix(start.vals)) == c(N,G))){
		  init.vals<- as.vector(start.vals) 
		  }else{
		  if(is.numeric(start.vals) & length(as.numeric(start.vals))==N)
		  { 
		    init.vals<- as.vector(unMAP(start.vals))
		    }else{ 
		    stop("start.vals improperly specified: see help")
		   }
		}
	}
	
	lpstarts <- numeric(restarts)
	cnv.warn <- FALSE
	#X <-  t(X) #row major
	x.vec <- as.vector(t(X))
	
	for( r in 1:restarts )
	{	
		# call the EM algorithm
	  
		w <- .C( 	"BLCA_EM_FIT", 							as.integer(x.vec), 			
							as.integer(N),								as.integer(M), 
							as.integer(ncat),	 						as.double(hparam),
							as.integer(prior.init.type), 			as.double(delta),
							as.double(gamma),
							as.integer(G),								as.integer(init),
							as.double(init.vals), 					as.integer(iter),
							iters = as.integer(iters),				group.probs = as.double(group.probs),
							weights = as.double(weights),			variable.probs = as.double(variable.probs),
							as.integer(model.indicator),			log.post = as.double(log.post),
							as.integer(MAP),							as.double(conv),
							eps = as.double(eps),
							converged = as.integer(converged),	log.object.max = as.double(log.object.max),
							as.integer(FALSE),
							PACKAGE = "BayesLCA" )
	
		log.post.this <- w$log.post[ w$iters ]
		
		lpstarts[r] <- log.post.this
		
		#store the run that gives the highest value  of the log posterior
		#	for the runs that have converged
		
		new.max <- FALSE
		if( log.post.this > log.post.max )
		{
			w.max <- w 
			log.post.max <- log.post.this
			new.max <- TRUE
		}
		
		if( verbose ) 
		{
		  if( MAP ) str.obj <- ", logpost = " else str.obj <- ", loglik = "
			if( new.max && r>1 )
			{
				cat( "\nNew maximum found... Restart number ",r,", logpost = ", round(log.post.this,2),"...", sep = "" )
			}else{ 
				cat( "\nRestart number ",r, str.obj, round(log.post.this,2),"...", sep = "" )
			}
		}
		
		if( as.logical(w$converged) == FALSE ) cnv.warn <- TRUE
	
	}	
	
	if(cnv.warn) warning("Some restarts failed to converge: rerun with a higher iter value or less stringent tolerance", call.=FALSE)
	if( verbose ) cat("\n")
	
	x$G <- G
	x$classprob <- w.max$weights
	
	o <- order(x$classprob, decreasing=TRUE)
	x$classprob <- x$classprob[o]
	
	var.probs.l <- list()
	
	l <- 1
	for( j in 1:M )
	{
		if(j == 1)
		{
			gap <- 0
		}else{
			gap <- G * sum( ncat[1:(j-1)] )
		}
		#variable probabilities stacked by group and then iteration
		if( model.indicator[j] )
		{
			var.probs.l[[l]] <- matrix( w.max$variable.probs[(gap+1):(gap + G*ncat[j])] , nrow =  G, ncol=ncat[j], byrow=TRUE )
			if(G>1){var.probs.l[[l]] <- var.probs.l[[l]][o,]}
			rownames( var.probs.l[[l]] ) <-  paste( "Group", 1:G )
			colnames( var.probs.l[[l]] ) <- levnames[[j]] 
			l <- l+1
		}
	}
	
	if( is.null(colnames(X)) )
	{
		names( var.probs.l ) <- paste("Variable", which( model.indicator==1 ) )
	}else{
		names( var.probs.l ) <- colnames(X)[  which( model.indicator == 1 ) ]
	}
	
	if( any( unlist(var.probs.l) < 0 ) |  any(x$classprob < 0) ) stop("MAP estimation returned negative probabilities (this can happen for some priors): use of default priors is recommended here")
	
	x$itemprob <- var.probs.l	
	
	vec.itemprobs <- unlist(x$itemprob)
	
	itemprobs.group.ind <- rep(paste("Group", 1:G), times = sum(ncat[ which(model.indicator == 1) ]))
	itemprobs.var.ind <- rep(names(x$itemprob), times = G * ncat[ which(model.indicator == 1) ])
	itemprobs.cat.ind <- paste("Cat", rep(as.numeric(unlist(apply(t(ncat[ which(model.indicator == 1) ]), 2, function(x) 0:(x-1)))), each = G))
	
	if(any(ncat > 2)){
		x$itemprob.tidy <- data.frame(itemprob = vec.itemprobs, group = itemprobs.group.ind, 
		                              variable = itemprobs.var.ind, category = itemprobs.cat.ind)
		# need to do some extra work here on itemprob.tidy to match with levnames
	}else{ 
	  x$itemprob <- matrix(vec.itemprobs[itemprobs.cat.ind == "Cat 1"], nrow = G, 
	                       ncol = sum(model.indicator), dimnames = list(paste("Group", 1:G), names(x$itemprob)))
	}

	x$Z <- matrix( w.max$group.probs, nrow=N, ncol=G )
	x$Z <- x$Z[,o]
	if(G>1){							      
	colnames(x$Z) <- paste( "Group", 1:G )	
	}	
  
  likl <- w.max$log.object.max
  
  # print warning message for BIC/AIC if priors different from default
  npars <- G * sum( model.indicator * ( ncat - 1 ) ) + G - 1
  likl <-  blca.compute.log.post( X, ncat, x$classprob, x$itemprob, model.indicator, prior=NULL )
  if( MAP ) x$logpost <- log.post.max else x$loglik <- likl
	
	x$BIC<- 2*likl- npars*log(N) 
	x$AIC<- 2*likl - 2*npars
	x$iter<- w.max$iters
	if( MAP ) x$poststore <- w.max$log.post[ 1:w.max$iters ] else x$loglikstore <- w.max$log.post[ 1:w.max$iters ]
	x$converged <- as.logical(w.max$converged)
	x$eps <- w.max$eps
	x$conv <- conv
	x$lpstarts <- lpstarts
	
	
	if( MAP ) x$prior <- list( alpha=alpha, beta=beta, delta=delta )
	x$model.indicator <- model.indicator
	x$MAP <- MAP
	
	x$ncat <- ncat
	
	x$for.boot <- for.boot
	if( sd ) x$for.boot <- TRUE
	if( x$for.boot )
	{
	  if(any(ncat>2))
	  {
	    vp <- lapply( x$itemprob, t )
	  }else{
	    vp <- vector( ncol(x$itemprob), mode="list")
	    for( k in seq_along(vp) ) vp[[k]] <- t( cbind(1-x$itemprob[,k],x$itemprob[,k]) )
	  }
	  # have to do this after the reordering
	  vp <- unlist( vp )
		x$boot.init = list( weights=x$classprob, var.probs=vp )
	}
	
	if( sd )
	{
		if( verbose ) cat("\nObtaining parameter uncertainty estimates...")
		if( sd.method[1] == "boot" ) y <- blca.boot( X, G=G, fit=x, iter=500, B=100, verbose=FALSE, conv=1e-4  )
		# still need to incorporate  the prior here
	  if( sd.method[1] == "delta" ) y <- blca.em.sd( X, ncat, x$classprob, var.probs.l, model.indicator, NULL, 1e-10 )
		x$itemprob.sd <-y$itemprob.sd
		names(x$itemprob.sd) <- names(x$itemprob)
		x$classprob.sd <- y$classprob.sd
		if( verbose ) cat("\nParameter uncertainty estimates complete...\n")
	}
	
	x$method <- "em"
	
	if( MAP & !for.boot ) warning("AIC and BIC will not be correct for model comparison if priors changed from default: set MAP to FALSE for MLE", call. = FALSE)
	
	x <- blca.return.order( x )

	if(any(ncat > 2)) class(x)<- c("blca.em", "blca.multicat","blca") else class(x)<-c("blca.em", "blca")

		return(x)
}
