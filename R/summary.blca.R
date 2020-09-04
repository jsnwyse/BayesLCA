# needs work for the collapsed

summary.blca <- function(object, ...){
	sum0<-sum2<-NULL
	
	pr <- FALSE
	if( object$method %in% c("vb","gibbs","collapsed") ){
		pr <- TRUE
	}else{
		if( object$MAP ) pr <- TRUE
	}
	
	if( pr ){
		sum2$ItemProb <- list(alpha=object$prior$alpha)
		sum2$ItemProb$beta<- object$prior$beta 
		sum2$ClassProb<- object$prior$delta
	}

	sum0$method <- object$method
	sum0$method.str <- object$method.str
	sum0$MAP <- object$MAP
	sum0$printnames<- object$printnames
	sum0$sum1<- object$sum1
	sum0$sum2<- sum2
	class(sum0)<-"summary.blca"
	return(sum0)
	}
