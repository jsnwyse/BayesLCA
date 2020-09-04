summary.blca.em <-
function(object, ...){
  
  object$method.str<- "EM algorithm"
  
  if( object$MAP ) a <- object$logpost else a <- object$loglik
  
	sum1<- c(object$iter, object$eps, a,  object$AIC, object$BIC )
	names(sum1)<- c("IterNumber", "ConvergenceDiff", "Log-Likelihood", "AIC", "BIC")
	object$printnames<- c("Number of iterations:","Log-Likelihood Increase at Convergence:", "Log-Likelihood:", "AIC:", "BIC:")
	
	if( object$MAP ){
	  object$method.str<- "EM (MAP) algorithm"
	  names(sum1)[3] <- "Log-Posterior"
		object$printnames[2:3]<- c("Log-Posterior Increase at Convergence:", "Log-Posterior:")
		sum1[3] <- object$logpost
	}
	object$sum1 <- sum1
	NextMethod("summary")
	}
