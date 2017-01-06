summary.blca.collapsed <-
function(object, ...){
	sum1<- c( object$iter, object$burn.in, object$thin, mean(object$samples$logpost), !object$fixed.G, object$variable.selection)
	names(sum1)<- c("IterNumber", "Burn-in", "ThiningRate", "LogPosterior", "GSearch","VSelect" )
	
	object$method<- "Collapsed Gibbs Sampling"
	object$printnames<- c("Chain Lengths:", "Burn-In:", "Thinning Rate:", "Mean Log-Posterior", "Latent class selection:", "Variable selection:")
	object$sum1<- sum1

	NextMethod("summary")
	}
