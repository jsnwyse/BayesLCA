#Functions for the fitting of Latent Class Analysis models using MCMC
#	methods. 
#	
#	Author:	Jason Wyse,
#			School of Computer Science and Statistics,
#			Lloyd Institute,
#			Trinity College,
#			Dublin 2,
#			Ireland.
#			mailto: wyseja@tcd.ie
#
#	Last modification of this code: Wed 11 May 2016 03:15:02 IST  			
#

blca.collapsed.post.hoc.estimates <- function( x, variables, G )
{

	#takes the variables to include and a given number of groups, runs a fixed sampler
	#	and computes the post-hoc parameter estimates in a call to C function... 
	
	if( !all.equal( class( x ) ,c("blca.collapsed","blca") ) ) stop("\t post-hoc parameter estimates may only be computed for objects of type blca.collapsed")
	
	X <- x$X
	N <- nrow(X)
	
	#make a sub-matrix of observations for the vars...
	
	if( max( variables ) > ncol( x$X ) )  stop("\t please check your variables argument- variables should be indexed 1,..., ncol(X)") 
	
	M <- length(variables)
	
	X <- X[ , variables ]
	
	num.categories <- x$num.categories[ variables ]
	
	alpha <- x$prior$alpha
	delta <- x$prior$delta
	beta <- x$prior$beta
	prob.inc <- x$prior$prob.inc
	iter <- x$iter
	burn.in <- x$burn.in
	thin <- x$thin
	
	hparam <- c( delta, beta )
	
	m <- blca.collapsed( X , G, alpha, beta, delta, num.categories, iter, burn.in, thin, fixed.G = TRUE, just.gibbs.updates = TRUE, post.hoc.est = FALSE, n.gibbs = nrow(X), prior.G = 1, G.max =  30, variable.selection = FALSE, prob.inc = .5, hprior.model = FALSE, relabel = TRUE, silent = TRUE )
	

	  #mcmc = collapsed.lca( Ynew, ncat, alpha,beta,fixed.groups = TRUE,just.gibbs.updates = TRUE,initial.groups = G,max.groups=G,n.iterations = n.iterations,n.burn = n.burn, thin.by = n.thin,prior.G=1,variable.selection = 0,prior.prob.include=.5)
	
	xx <- list()
	
	xx$G = G
	
	xx$variables = variables
	
	for(i in 1:length(variables)){	
	
		estimate = numeric( G * num.categories[i] )
		sd.estimate = numeric( G * num.categories[i] )
		classprob.estimate = numeric( G )
		sd.classprob.estimate = numeric( G )
		
		w = .C(		"R_LCA_VS_COMPUTE_POST_HOC_PARAMETER_ESTIMATES",					as.integer(X),
					as.integer(N),														as.integer(M),
					as.integer(num.categories),													as.double(hparam),
					as.integer(G),														as.integer((iter-burn.in)/thin),
					as.integer(m$labelstore-1),											as.integer(m$samples$var.ind),
					as.integer(i-1),													est = as.double(estimate),
					sd.est = as.double(sd.estimate),									classest = as.double(classprob.estimate),
					sd.classest = as.double(sd.classprob.estimate),
					PACKAGE = "BayesLCA" )
		
		
		xx$itemprob[[i]] = matrix(w$est,nrow=G,ncol=num.categories[i],byrow=FALSE)
		xx$sd.itemprob[[i]] = matrix(w$sd.est,nrow=G,ncol=num.categories[i],byrow=FALSE)
		xx$classprob = w$classest
		xx$sd.classprob = w$sd.classest
		
		#do row and col names for clarity
		
		gr = rep("Class",G)
		no = 1:G
		rownames(xx$itemprob[[i]]) = paste(gr,no)
		rownames(xx$sd.itemprob[[i]]) = paste(gr,no)
		
		ca = rep("Category",num.categories[i])
		nos = 1:num.categories[i]
		colnames(xx$itemprob[[i]]) = paste(ca,nos)
		colnames(xx$sd.itemprob[[i]]) = paste(ca,nos)
		
	}
	xx$itemprob.sd <- xx$sd.itemprob
	xx$classprob.sd <- xx$sd.classprob
	return(xx)
}

