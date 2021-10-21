blca <-
function(X, G, formula=NULL, ncat=NULL, method=c("em", "gibbs", "boot", "vb", "collapsed"),...){

	method<- match.arg(method)
	
	# check for missing values
	miss <- blca.check.missing(X)
	if( miss$missing & method != "gibbs" ) warning("X contains missing values: method = 'gibbs' will carry out data imputation.", call.=FALSE) 
	rm(miss)
	
	if(method == "em") return(blca.em(X, G, formula, ...))
	if(method == "vb") return(blca.vb(X, G, formula, ...))
	if(method == "gibbs") return(blca.gibbs(X, G, formula, ...))
	if(method == "boot") return(blca.boot(X, G, formula, ...))
	if(method == "collapsed") return(blca.collapsed(X, G, formula, ...))
	}
