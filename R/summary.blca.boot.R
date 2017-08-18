summary.blca.boot <-
function(object, ...){
	sum1<- c( object$B, object$logpost, object$AIC, object$BIC)
	names(sum1)<- c("IterNumber", "Log-Posterior", "AIC", "BIC")
	
	object$method<- "Bootstrap"
	object$printnames<- c("Number of Samples:", "Log-Posterior:", "AIC:", "BIC:")
	object$sum1<- sum1

	NextMethod("summary")
	}
