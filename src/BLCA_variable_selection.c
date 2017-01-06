/*Functions for the fitting of Latent Class Analysis models using MCMC
	methods. Two implementations of the model are included in the Bayesian
	formulation: collapsed and not collapsed.
	
	Author:	Jason Wyse,
			School of Computer Science and Statistics,
			Lloyd Institute,
			Trinity College,
			Dublin 2,
			Ireland.
			mailto: wyseja@tcd.ie
			
	Last modification of this code: Mon 09 May 2016 15:38:03 IST   */
	
#include "BLCA_variable_selection.h"

void BLCA_update_model_by_variable_include_exclude(struct mix_mod *mixmod, int *accepted,int *proposed,int var)
{
  
	int i,j,S=0,D=0,Sn,Dn;
	double logP=0.,logPn,cstr;
	double *logprobs,logprob_undiscriminating,logratio,lp_current,lp_proposed,l_ratio_prior_variable_inclusion;
	
	for(j=0;j<mixmod->d;j++){
	  if(mixmod->varindicator[j]){
	    S += mixmod->ncat[j];
	    logP += log(mixmod->ncat[j]);
	    D += 1;
	  }
	}
	
	
	if(mixmod->varindicator[var]){
	  
	  /*this is the removal case-- check constraint and return if not possible*/
	  Sn = S - mixmod->ncat[var];
	  logPn = logP - log(mixmod->ncat[var]);
	  Dn = D -1;
	  
	  cstr = logPn - log(( Sn - Dn + 1 )*mixmod->G);
	  
	 // Rprintf("\nIn the removal case the cstr was %d",cstr);
	  
	  if(!(cstr>0)){
	    return;
	  }
	  
	  mixmod->varindicator[var] = 0;
	  
	  l_ratio_prior_variable_inclusion = BLCA_l_prior_variable_include(D,mixmod) - BLCA_l_prior_variable_include(D+1,mixmod);
	  
	}else{
	  
	  /*the is the inclusion case-- check constraint and return if not possible*/
	  Sn = S + mixmod->ncat[var];
	  logPn = logP + log(mixmod->ncat[var]);
	  Dn = D + 1;
	  
	  cstr = logPn - log(( Sn - Dn + 1 )*mixmod->G);
	  
	 // Rprintf("\nIn the add case the cstr was %d",cstr);
	  
	  if(!(cstr > 0)){
	    return;
	  }
	  
	  mixmod->varindicator[var] = 1;
	  
	  l_ratio_prior_variable_inclusion = BLCA_l_prior_variable_include(D,mixmod) - BLCA_l_prior_variable_include(D-1,mixmod);
	  
	}
	
	*proposed += 1;
  
	logprobs = calloc(mixmod->G,sizeof(double));
	
	/*temporary storage of the current log-probs*/
	lp_current = 0.;
	lp_proposed = 0.;
	for(i=0;i<mixmod->G;i++){
		logprobs[i] = mixmod->components[ mixmod->whereis[i] ]->log_prob;
		lp_current += logprobs[i];
		BLCA_recompute_marginal_likelihood_component(mixmod->components[ mixmod->whereis[i] ],mixmod);
		lp_proposed += mixmod->components[ mixmod->whereis[i] ]->log_prob;
	}	
	
	logprob_undiscriminating = mixmod->undiscriminating->log_prob;
	
	BLCA_recompute_marginal_likelihood_undiscriminating_variables(mixmod->undiscriminating,mixmod);
	
	logratio = lp_proposed - lp_current + mixmod->undiscriminating->log_prob - logprob_undiscriminating + l_ratio_prior_variable_inclusion;
	
	if(log(runif(0.0,1.0))<logratio){
	  
	  *accepted += 1;
	  
	}else{
	  
	  if(mixmod->varindicator[var]){
	    mixmod->varindicator[var] = 0;
	  }else{
	    mixmod->varindicator[var]=1;
	  }
	  
	  for(i=0;i<mixmod->G;i++){
	    mixmod->components[ mixmod->whereis[i] ]->log_prob = logprobs[i];
	  }
	  
	  mixmod->undiscriminating->log_prob = logprob_undiscriminating;
	  
	}
	
	free(logprobs);
  
  
	return;
  
}

