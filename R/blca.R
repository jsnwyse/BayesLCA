blca <-
function(X, G, ncat=NULL, method=c("em", "gibbs", "boot", "vb", "collapsed"),...){
	# if(length(method)>1){
		# method<- method[1]
		# warning(paste("More than one method specified", method, "method used only."))
		# } 
	method<- match.arg(method)
	
	if(method == "em") return(blca.em(X, G, ...))
	if(method == "vb") return(blca.vb(X, G, ...))
	if(method == "gibbs") return(blca.gibbs(X, G, ...))
	if(method == "boot") return(blca.boot(X, G, ...))
	if(method == "collapsed") return(blca.collapsed(X, G, ...))
	}
