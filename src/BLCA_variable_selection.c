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
 
 Last modification of this code: Jul 15 2020   */

#include "BLCA_variable_selection.h"

void BLCA_update_model_by_variable_include_exclude_old(struct mix_mod *mixmod, int *accepted,int *proposed,int var)
{
  
  int i, j, d_in = 0, G = mixmod->G, d = mixmod->d, *ncat = mixmod->ncat, *vind = mixmod->varindicator, *wis = mixmod->whereis;
  //S=0,D=0,Sn,Dn;
  double cstr; //logP=0.,logPn,cstr;
  double *logprobs,logprob_undiscriminating,logratio,lp_current,lp_proposed,l_ratio_prior_variable_inclusion;
  struct component *comp;
  /*for(j=0;j<mixmod->d;j++){
   if(mixmod->varindicator[j]){
   S += mixmod->ncat[j];
   logP += log(mixmod->ncat[j]);
   D += 1;
   }
  }*/
  
  for( j=0; j<d; j++ ) d_in += vind[j];
  
  // check edge cases
  if( d_in == 1 & vind[ var ] == 1 ) return;
  
  
  if( vind[ var ] == 1 )
  {
    
    /*this is the removal case-- check constraint and return if not possible*/
    /*Sn = S - mixmod->ncat[var];
     logPn = logP - log(mixmod->ncat[var]);
     Dn = D -1;
     
     cstr = logPn - log(( Sn - Dn + 1 )*mixmod->G);*/
    
    // Rprintf("\nIn the removal case the cstr was %d",cstr);
    
    /*if(!(cstr>0)){
     return;
     }*/
    
    vind[ var ] = 0;
    
    l_ratio_prior_variable_inclusion = BLCA_l_prior_variable_include( d_in - 1, mixmod ) - BLCA_l_prior_variable_include( d_in , mixmod );
    
  }else{
    
    /*the is the inclusion case-- check constraint and return if not possible*/
    /* Sn = S + mixmod->ncat[var];
     logPn = logP + log(mixmod->ncat[var]);
     Dn = D + 1;
     
     cstr = logPn - log(( Sn - Dn + 1 )*mixmod->G);*/
    
    // Rprintf("\nIn the add case the cstr was %d",cstr);
    
    /*if(!(cstr > 0)){
     return;
     }*/
    
    vind[var] = 1;
    
    l_ratio_prior_variable_inclusion = BLCA_l_prior_variable_include( d_in + 1, mixmod ) - BLCA_l_prior_variable_include( d_in, mixmod );
    
  }
  
  *proposed += 1;
  
  logprobs = (double *)calloc( G, sizeof(double) );
  
  /*temporary storage of the current log-probs*/
  lp_current = 0.;
  lp_proposed = 0.;
  for( i=0; i<G; i++ )
  {
    comp = mixmod->components[ wis[i] ] ;
    logprobs[i] = comp->log_prob;
    lp_current += comp->log_prob;
    mixmod->component_compute = i ;
    BLCA_recompute_marginal_likelihood_component( comp, mixmod );
    lp_proposed += comp->log_prob;
  }	
  
  logprob_undiscriminating = mixmod->undiscriminating->log_prob;
  
  BLCA_recompute_marginal_likelihood_undiscriminating_variables( mixmod->undiscriminating, mixmod );
  
  logratio = lp_proposed - lp_current + mixmod->undiscriminating->log_prob - logprob_undiscriminating + l_ratio_prior_variable_inclusion;
  
  if(log(runif(0.0,1.0))<logratio)
  {
    
    *accepted += 1;
    
  }else{
    
    if( vind[var] == 1 ) vind[var] = 0; else vind[var]=1;
    
    
    for( i=0; i<G; i++ )
      mixmod->components[ wis[i] ]->log_prob = logprobs[i];
    
    mixmod->undiscriminating->log_prob = logprob_undiscriminating;
    
  }
  
  free(logprobs);
  
  
  return;
  
}

