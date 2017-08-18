summary.blca.collapsed <-
function(object, ...){
	sum1<- c( object$iter, object$burn.in, object$thin, mean(object$samples$logpost), object$G.sel, object$var.sel )
	names(sum1)<- c("IterNumber", "Burn-in", "ThiningRate", "LogPosterior", "GSelection","VSelection" )
	
	object$method<- "Collapsed Gibbs Sampling"
	object$printnames<- c("Chain Lengths:", "Burn-In:", "Thinning Rate:", "Mean Log-Posterior", "Latent class selection:", "Variable selection:")
	object$sum1<- sum1

	NextMethod("summary")
	}
