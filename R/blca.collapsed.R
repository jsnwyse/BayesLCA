blca.collapsed <- function( X,  G, ncat=NULL , alpha=1, beta=1, delta=1, start.vals=c("single"), counts.n=NULL, iter=5000, burn.in=1000, thin=1, G.sel=FALSE, var.sel=FALSE, post.hoc.run=TRUE, n.gibbs=nrow(X), only.gibbs=TRUE, G.max=30, G.prior=dpois(1:G.max, lambda=1), model.indicator=NULL, prob.inc=.5, hprior.model=FALSE, relabel=TRUE, verbose=TRUE, verbose.update=1000 )
{

	t1 <- proc.time()
	
	# convert X into a numeric matrix and check inputs
	D <- blca.check.data( X, counts.n, ncat )
	
	X <- D$X
	ncat <- D$ncat

	N <- nrow(X)
	M <- ncol(X)
	
	if( G.sel == FALSE ) G.max <- G
	
	out.prior <- blca.check.prior.collapsed( alpha, beta, delta, G.sel, G, G.max, M, ncat, only.gibbs )
	prior.init.type <- out.prior$prior.init.type
	gamma <- out.prior$gamma
	delta <- out.prior$delta
	only.gibbs <- out.prior$only.gibbs
	
	if( G.sel & G.max < G ){
		stop("The maximum number of groups cannot be less than the initial number of groups. Please rerun with G.max (default 30) appropriately specified.")
	}
	
	# thinning parameter
	if( thin > 1 || thin == 0 ) stop("Argument thin gives a thinning rate and should be less than or equal to 1.")
	
	Thinby <- floor(1/thin)

	stored <- floor( iter / Thinby )	

	memberships <- numeric(stored*N)
	variable.inclusion.indicator <- numeric(stored*M)
	num.groups <- numeric(stored)
	log.post <- numeric(stored)
	prior.include <- numeric(stored)
	
	if( is.null(model.indicator) ) model.indicator <- rep(1,M) 
	
	hparam <- c( delta[1], beta[1] )
	accrts <- numeric(6)
	
	w <- .C(	"BLCA_VS",																as.integer(X),
				as.integer(N),																as.integer(M),
				as.integer(ncat),															as.double(hparam),
				as.integer(prior.init.type),									as.double(delta),
				as.double(gamma),
				as.integer(!G.sel),														as.integer(only.gibbs),
				as.integer(G),																as.integer(G.max),
				as.integer(iter),															as.integer(burn.in),
				as.integer(Thinby),														as.integer(n.gibbs),
				memberships = as.integer(memberships),				ngroups = as.integer(num.groups),
				as.double(G.prior),														as.integer(var.sel),
				variable.inclusion.indicator = as.integer(variable.inclusion.indicator),	as.double(prob.inc),
				log.posterior = as.double(log.post),					as.integer(hprior.model),
				prior.include = as.double(prior.include),			as.integer(model.indicator),
				as.integer(verbose),													as.integer(verbose.update),
				accrts=as.double(accrts),
				PACKAGE = "BayesLCA" )
				
	
	if( verbose ) cat("\nFinished sampling...\n")

	membership.mat <- matrix( w$memberships, nrow = stored, ncol=N, byrow=FALSE ) + 1
	vindicator.mat <- matrix( w$variable.inclusion.indicator, nrow=stored, ncol=M, byrow=FALSE )
	
	
	if( relabel ) 
	{
		if( G.sel ) grps <- w$ngroups else grps <- rep( G , stored )
		relabelled <- undo.label.switching( membership.mat, grps, w$log.posterior )
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
	  ph.est <- blca.gibbs( X, G=Gstar, model.indicator=gamma.star, ncat=ncat, verbose=FALSE )
	    
	  x$classprob <- ph.est$classprob
	  x$classprob.sd <- ph.est$classprob.sd
	  x$itemprob <- ph.est$itemprob 
	  x$itemprob.sd <- ph.est$itemprob.sd 
	  x$itemprob.tidy <- ph.est$itemprob.tid

	  ## Name item probabilities
	  if( FALSE )
	  {
	  if(is.matrix(x$itemprob)){
	    if(is.null(colnames(X))) {
	      colnames( x$itemprob)<- paste("Col", Mstar)
	      } else{
	        colnames( x$itemprob)<- colnames(X)[Mstar]
	        } 
	    } else{
	      if(is.null(colnames(X))){
	        names(x$itemprob)<- paste("Col", Mstar)
	        } else{
	          names( x$itemprob)<- colnames(X)[Mstar]
	          } 
	      }
	    }
	} 
	
	
	if( any(ncat > 2 ) ) class(x) <- c("blca.collapsed", "blca.multicat", "blca") else class(x)<-c("blca.collapsed", "blca")
	
	return(x)
	
}
