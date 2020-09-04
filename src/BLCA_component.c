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


#include "BLCA_component.h"

struct component * BLCA_create_component( struct mix_mod *mixmod )
{
  struct component *comp = (struct component *)malloc(sizeof(struct component));
  BLCA_allocate_component( comp, mixmod );
  return( comp );
}

void BLCA_allocate_component(struct component *component,struct mix_mod *mixmod)
{
  /*allocate memory for a component*/
  int j, d=mixmod->d, *ncat=mixmod->ncat;
  
  component->n_g = 0;
  
  component->N = (int **)calloc(d,sizeof(int *));
  
  for( j=0; j<d; j++ ) component->N[j] = (int *)calloc(ncat[j],sizeof(int));
  
  if( !mixmod->collapsed )
  {
    component->prob_variables = (double **)calloc(d,sizeof(double *));
    for( j=0; j<d; j++ ) component->prob_variables[j] = (double *)calloc(ncat[j],sizeof(double));
  }
  
  if( mixmod->VB )
  {
    component->beta_ud = (double **)calloc( d, sizeof(double *));
    component->di_beta_ud = (double **)calloc( d, sizeof(double *));
    component->di_sum_beta_ud = (double *)calloc( d, sizeof(double));
    for( j=0; j<d; j++ )
    {
      component->beta_ud[j] = (double *)calloc(ncat[j],sizeof(double));
      component->di_beta_ud[j] = (double *)calloc(ncat[j],sizeof(double));
    }
  }
  
  return;
}

void BLCA_free_component(struct component *component,struct mix_mod *mixmod)
{
  
  int j;
  
  for(j=0;j<mixmod->d;j++){
    free(component->N[j]);
  }
  free(component->N);
  
  if( !mixmod->collapsed )
  {
    for(j=0;j<mixmod->d;j++){
      free(component->prob_variables[j]);
    }
    free(component->prob_variables);
  }
  
  if( mixmod->VB )
  {
    for(j=0;j<mixmod->d;j++){
      free(component->beta_ud[j]);
      free(component->di_beta_ud[j]);
    }	
    free(component->beta_ud);
    free(component->di_beta_ud);
    free(component->di_sum_beta_ud);
  }
  
  return;
  
}

void BLCA_copy_component(struct component *component_original,struct component *component_target,struct mix_mod *mixmod)
  /*copy the contents of the first argument into the second component argument*/
{
  
  int i, j, d=mixmod->d, *ncat=mixmod->ncat, *N_t, *N_o;
  double *a_t, *a_o;
  
  component_target->n_g = component_original->n_g;
  for( i=0; i<d; i++ )
  {
    N_t = component_target->N[i]; N_o = component_original->N[i];
    for( j=0; j<ncat[i]; j++ ) N_t[j] = N_o[j];
  }
  
  if( !mixmod->collapsed )
  {
    for( i=0; i<d; i++ )
    {
      a_t = component_target->prob_variables[i];
      a_o = component_original->prob_variables[i];
      for( j=0; j<ncat[i]; j++ ) a_t[j] = a_o[j];
    }
  }
  
  if( mixmod->VB )
  {
  
    for( i=0; i<d; i++ )
    {
      a_t = component_target->beta_ud[i];
      a_o = component_original->beta_ud[i];
      for( j=0; j<ncat[i]; j++ ) a_t[j] = a_o[j];
    }
    
    for( i=0; i<d; i++ )
    {
      a_t = component_target->di_beta_ud[i];
      a_o = component_original->di_beta_ud[i];
      for( j=0; j<ncat[i]; j++ ) a_t[j] = a_o[j];
    }
     
    for( i=0; i<d; i++ ) component_target->di_sum_beta_ud[i] = component_original->di_sum_beta_ud[i]; 
    
  }
  
  component_target->log_prob = component_original->log_prob;
  
  return;
  
}


void BLCA_add_to_component( struct component *component, int *x, struct mix_mod *mixmod, int sgn )
{
  int j, *N ;
  
  component->n_g += sgn ;
  
  for( j=0; j<mixmod->d; j++ ) 
  {
    N = component->N[j] ;
    N[ x[j] ] += sgn ; 
  }
  
  return;
}

void BLCA_add_components( struct component *component_original, struct component *component_target, struct mix_mod *mixmod )
{
  // add two components together- result is returned  in target
  int c, j, *N_t, *N_o, *ncat = mixmod->ncat;
  
  component_target->n_g += component_original->n_g;
  
  for( j=0; j<mixmod->d; j++ )
  {
    N_t = component_target->N[j];
    N_o = component_original->N[j];
    for( c=0; c<ncat[j]; c++ ) N_t[c] += N_o[c]; 
  }
  
  return;
}


void BLCA_recompute_sufficient_statistics_for_components(struct mix_mod *mixmod)
{
  
  int i,j,k;
  
  /*undiscriminating variables*/
  
  
  
  /*discriminating variables*/
  
  for(k=0;k<mixmod->G;k++){	
    /*reset*/
    for(j=0;j<mixmod->d;j++){
      for(i=0;i<mixmod->ncat[j];i++){
        mixmod->components[k]->N[j][i] = 0;
        mixmod->undiscriminating->N[j][i] = 0;
      }
    }
    mixmod->components[k]->n_g = 0;
  }
  
  int z;
  
  for(i=0;i<mixmod->n;i++){
    
    z = mixmod->z[i];	
    
    mixmod->components[ z ]->n_g += 1;
    
    for(j=0;j<mixmod->d;j++){
      mixmod->components[ z ]->N[j][ mixmod->Y[j][i] ] += 1;
      mixmod->undiscriminating->N[j][ mixmod->Y[j][i] ] += 1;
    }
    
    
  }
  
  return;
}


double BLCA_compute_log_data_probability_with_inclusion_in_component(int *x,struct component *component,struct mix_mod *mixmod)
{
  
  int i, j, c, k=mixmod->component_compute, d=mixmod->d, *vind=mixmod->varindicator, *ncat=mixmod->ncat, *N;
  double log_prob = 0., s, ind, *beta;
  
  // fix beta issue here...
  
  /*if( mixmod->prior_type == 0 )
   {
   
   for( j=0; j<d; j++ )
   {
   
   if( vind[j] )
   {
   
   log_prob += lgamma( ncat[j] * mixmod->beta ) - ncat[j] * lgamma( mixmod->beta ) - lgamma( component->n_g + 1.0 + ncat[j] * mixmod->beta );
   
   for( i=0; i<ncat[j]; i++ )
   {
   
   I = (x[j] == i) ?  1 : 0 ;
   
   log_prob += lgamma(component->N[j][i] + I + mixmod->beta);
   
   }
   
   }
   
   }
   
   }*/
  
  if( mixmod->prior_type > 0 )
  {
    for( j=0; j<d; j++ )
    {
      if( vind[j] )
      {
        s =  0.; 
        beta = mixmod->beta_prior[k][j] ;
        N = component->N[j];
        for( c=0; c<ncat[j]; c++ )
        { 
          ind = (x[j] == c) ? 1.0 : 0.0 ;
          log_prob += lgamma( N[c] + ind + beta[c] ) /*- lgamma( beta[c] )*/ ;
          s += beta[c] ;
        } 
        log_prob += mixmod->lg_sum_beta[k][j] - mixmod->lg_beta_sum[k][j] - lgamma( component->n_g + 1. + s ) ; 
      }
    }
  }
  
  return(log_prob);
}

double BLCA_compute_log_marginal_likelihood_with_inclusion_in_component( int *x, struct component *component, struct mix_mod *mixmod )
{
  // function used for class prediction
  
  int i, j, d=mixmod->d, *ncat=mixmod->ncat, *vind=mixmod->varindicator, *N;
  double ind, log_prob = 0.0, alpha = mixmod->alpha_prior[0], beta = mixmod->beta_prior[0][0][0] ;
  
  log_prob = lgamma( component->n_g + 1.0 + alpha );
  
  for( j=0; j<d; j++ )
  {
    
    if( vind[j] )
    {
      
      log_prob += lgamma( ncat[j] * beta ) - ncat[j] * lgamma(beta) - lgamma( component->n_g + 1.0 + ncat[j]*beta );
      N = component->N[j];
      for( i=0; i<ncat[j]; i++ )
      {
        ind = (x[j] == i) ?  1.0 : 0.0 ;
        log_prob += lgamma(N[i] + ind + beta);
      }
    }
    
  }
  
  return(log_prob);
}


double BLCA_compute_log_data_probability_component(struct component *component,struct mix_mod *mixmod)
{
  
  int i, j, c, k=mixmod->component_compute, d=mixmod->d, *ncat=mixmod->ncat, *vind=mixmod->varindicator, *N;
  double log_prob = 0., s, *beta;
  
  if( mixmod->prior_type ==  0 )
  {
    
    for( j=0; j<d; j++ )
    {
      
      if( vind[j] )
      {
        
        log_prob += lgamma( ncat[j] * mixmod->beta ) - ncat[j]*lgamma(mixmod->beta) - lgamma(component->n_g + ncat[j]*mixmod->beta );
        N = component->N[j];
        for( i=0; i<ncat[j]; i++ ) log_prob += lgamma(N[i] + mixmod->beta);
      }
      
    }
    
  }
  
  if( mixmod->prior_type > 0 )
  {
    for( j=0; j<d; j++ )
    {
      if( vind[j] )
      { 
        s = 0.;
        beta = mixmod->beta_prior[k][j];
        N = component->N[j];
        for( c=0; c<ncat[j]; c++ )
        { 
          log_prob += lgamma( N[c] + beta[c] ) - lgamma( beta[c] ) ;
          s += beta[c] ;
        } 
        log_prob += lgamma(s) - lgamma( component->n_g + 1.0 + s ) ; 
      } 
    }
  }
  
  return(log_prob);	
}

double BLCA_return_log_marginal_likelihood_component( struct component *component, struct mix_mod *mixmod )
{ 
  // these quantities will save some time for the Gibbs update of the labels
  int i,j,c, k, wis_k, d = mixmod->d, *vind = mixmod->varindicator, *ncat = mixmod->ncat, *N ;
  double log_prob, s, r, a, *beta;
  
  if( mixmod->prior_type == 0 )
  { 
    log_prob = lgamma( component->n_g + mixmod->alpha ); 
    if(component->n_g > 0)
    { 
      for( j=0; j<d; j++ )
      { 
        if( vind[j] )
        { 
          log_prob += lgamma( ncat[j] * mixmod->beta ) - ncat[j] * lgamma( mixmod->beta ) - lgamma( component->n_g + ncat[j] * mixmod->beta );
          N = component->N[j];
          for( i=0; i<ncat[j]; i++ )
          { 
            log_prob += lgamma( N[i] + mixmod->beta );
          }
        }
      }
    }	 
  }
  
  if( mixmod->prior_type > 0 )
  {
    
    k = mixmod->component_compute;
    log_prob = lgamma( component->n_g + mixmod->alpha_prior[k] ) ;
    if(component->n_g > 0)
    { 
      for( j=0; j<d ; j++ )
      { 
        if( vind[j] )
        {  
          s = 0.; //r = 0.;
          beta = mixmod->beta_prior[k][j] ;
          N = component->N[j];
          for( c=0; c<ncat[j]; c++ )
          {  
            a = N[c] + beta[c] ;
            // compile these into sums?
            log_prob += lgamma(a);// - lgamma(beta[c]) ;
            s += a ;
            //r += beta[c] ; 
          }  
          log_prob -= mixmod->lg_beta_sum[k][j];
          log_prob += mixmod->lg_sum_beta[k][j] - lgamma(s) ;
        }	
      }
    }
  }  
  
  return( log_prob );
}


void BLCA_recompute_marginal_likelihood_component(struct component *component,struct mix_mod *mixmod)
{
  
  double log_prob = BLCA_return_log_marginal_likelihood_component( component, mixmod );
  
  component->log_prob = log_prob;
  
  return;
}


double BLCA_return_marginal_likelihood_undiscriminating_variables(struct component *undiscriminating,struct mix_mod *mixmod)
{
  
  int i, j, d=mixmod->d, *vind=mixmod->varindicator, *ncat=mixmod->ncat, *N;
  double log_prob, s=0., r=0., *beta;
  
  log_prob = 0.;
  
  /*if( mixmod->prior_type == 0 )
   {
   for( j=0; j<d; j++ )
   {
   if( !vind[j] )
   { //variables not in mixture model 
   log_prob += lgamma( ncat[j] * mixmod->beta ) - ncat[j] * lgamma( mixmod->beta ) - lgamma( undiscriminating->n_g + mixmod->ncat[j] * mixmod->beta );
   for( i=0; i<ncat[j]; i++ )
   log_prob += lgamma(undiscriminating->N[j][i] + mixmod->beta);
   }
   } 
   }*/
  
  if( mixmod->prior_type > 0 )
  { 
    for( j=0; j<d; j++ )
    { 
      if( !vind[j] )
      {
        r=0.; s=0.;
        beta = mixmod->beta_prior[0][j];
        N = undiscriminating->N[j];
        for( i=0; i<ncat[j]; i++ )
        {
          r += N[i] + beta[i] ;
          log_prob += lgamma( N[i] + beta[i] ) ;
          //s +=  beta[i] ;
        }	 
        log_prob += /*lgamma(s)*/ mixmod->lg_sum_beta[0][j] - mixmod->lg_beta_sum[0][j] - lgamma(r) ; 
      } 
    }
  }
  
  return(log_prob);
  
}

void BLCA_recompute_marginal_likelihood_undiscriminating_variables(struct component *undiscriminating,struct mix_mod *mixmod)
{
  undiscriminating->log_prob = BLCA_return_marginal_likelihood_undiscriminating_variables( undiscriminating, mixmod );
  return;
}




