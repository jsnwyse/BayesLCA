# this is giving an error

summary.blca.collapsed <-
function(object, ...){
	sum1<- c( object$iter, object$burn.in, object$thin, object$G.sel, object$var.sel )
	names(sum1)<- c("IterNumber", "Burn-in", "ThiningRate", "GSelection","VSelection" )
	
	object$method.str <- "Collapsed Gibbs Sampling"
	object$printnames<- c("Chain Lengths:", "Burn-In:", "Thinning Rate:", "Latent class selection:", "Variable selection:")
	object$sum1<- sum1

	NextMethod("summary")
	}
