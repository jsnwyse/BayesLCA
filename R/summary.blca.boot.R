summary.blca.boot <-
function(object, ...){
	if( object$MAP) 
	  sum1<- c( object$B, mean(object$samples$logpost), object$logpost, object$AIC, object$BIC)
	else 
	  sum1<- c( object$B, mean(object$samples$loglik), object$loglik, object$AIC, object$BIC)
	
	names(sum1)<- c("IterNumber", "Log-Posterior-mean", "Log-Posterior-max", "AIC", "BIC")
	object$method.str<- "Bootstrap (MAP)"
	object$printnames<- c("Number of Samples:", "Log-Posterior Mean (over all samples):", "Log-Posterior at estimated values:", "AIC:", "BIC:")
	if(!object$MAP)
	{
	  object$method.str <- "Bootstrap"
	  names(sum1)[2] <- "Log-Likelihood-mean"
	  names(sum1)[3] <- "Log-Likelihood-max"
	  object$printnames[2] <- "Log-Likelihood Mean:"
	  object$printnames[3] <- "Log-Likelihood Max:"
	  sum1[2] <- mean(object$samples$loglik)
	  sum1[3] <- max(object$samples$loglik)
	}
	object$sum1<- sum1

	NextMethod("summary")
	}
