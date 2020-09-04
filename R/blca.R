blca <-
function(X, G, ncat=NULL, method=c("em", "gibbs", "boot", "vb", "collapsed"),...){

	method<- match.arg(method)
	
	# check for missing values
	miss <- blca.check.missing(X)
	if( miss$missing & method != "gibbs" ) warning("X contains missing values: method = 'gibbs' will carry out data imputation.", call.=FALSE) 
	rm(miss)
	
	if(method == "em") return(blca.em(X, G, ...))
	if(method == "vb") return(blca.vb(X, G, ...))
	if(method == "gibbs") return(blca.gibbs(X, G, ...))
	if(method == "boot") return(blca.boot(X, G, ...))
	if(method == "collapsed") return(blca.collapsed(X, G, ...))
	}
