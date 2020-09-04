blca.return.order <- function( x )
{
  nm <- names(x)
  nm.ord <- c( "method", "call", "args", "G", "classprob", "classprob.sd", "itemprob", "itemprob.sd", "itemprob.tidy", "Z", "parameters", "logpost", "loglik", "LB", 
               "BIC", "AIC", "DIC", "BICM", "AICM", "samples", "accrates", "B", "iter", "burn.in", "thin", "G.sel", "var.sel", "prob.inc", "hprior.model", "relabel", "poststore", "lbstore", "converged", "eps", "conv", "lpstarts", "convergence", "converged", "prior", "classprob.initial", "itemprob.initial", 
               "model.indicator", "MAP", "ncat", "for.boot", "boot.init")
  pos <- match( nm, nm.ord)
  nm.ror <- nm.ord[ sort(pos) ]
  x <- x[ nm.ror ]
  return(x)
}
