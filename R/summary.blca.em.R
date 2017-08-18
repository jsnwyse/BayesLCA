summary.blca.em <-
function(object, ...){	
	sum1<- c(object$iter, object$eps, object$logpost,  object$AIC, object$BIC )
	if( object$MAP ){
		names(sum1)<- c("IterNumber", "ConvergenceDiff", "Log-Posterior", "AIC", "BIC")
	}else names(sum1)<- c("IterNumber", "ConvergenceDiff", "Log-Likelihood", "AIC", "BIC")
	object$method<- "EM algorithm"
	if( object$MAP ){
		object$printnames<- c("Number of iterations:","Log-Posterior Increase at Convergence:", "Log-Posterior:", "AIC:", "BIC:")
	}else object$printnames<- c("Number of iterations:","Log-Likelihood Increase at Convergence:", "Log-Likelihood:", "AIC:", "BIC:")
	object$sum1<- sum1

	NextMethod("summary")
	}
