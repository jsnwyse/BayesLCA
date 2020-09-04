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

#include "BLCA_label_updates.h"

int BLCA_update_allocations_with_gibbs(struct mix_mod *mixmod)
{
  
  //printf("\nwithin Gibbs...\n");
  
  int i, ii, k, g, g_prime, n=mixmod->n, d=mixmod->d, position, G=mixmod->G, *order, *x, *y, *z, *wis = mixmod->whereis,  prior=mixmod->prior_type ;
  double *probs, *lpp_store, max, nrmcnst, lc, lcm, lp, lpp  ;
  
  probs = (double *)calloc(G,sizeof(double));
  lpp_store = (double *)calloc( G, sizeof(double) );
  order = (int *)calloc(mixmod->n,sizeof(int));
  
  struct component *comp_c, *comp_p ; 
  
  //cycle through the elements in turn... maybe randomize order?
  
  //randomize the order...
  for( i=0; i<n; i++ ) order[i] = i;
  
  //use gsl to shuffle
  BLCA_random_ranshuffle( order, n );
  
  y = mixmod->y;
  z = mixmod->z;
  
  for( i=0; i<mixmod->n_gibbs; i++ )
  {
    
    ii = order[i];
    g = z[ii];
    //point to the entry
    x = y + ii * d; //mixmod->Yobs[ind];
    
    comp_c = mixmod->components[ wis[g] ] ;
    
    lc = comp_c->log_prob;
    
    BLCA_add_to_component( comp_c, x, mixmod, -1 );
    
    if( prior == 1 ) mixmod->component_compute = g; // why is this here?
    
    lcm = BLCA_return_log_marginal_likelihood_component( comp_c, mixmod );
    
    //cycle through remaining groups
    for( k=0; k<G; k++ )
    {
      
      if( k!=g )
      {
        
        comp_p = mixmod->components[ wis[k] ] ;
        
        lp = comp_p->log_prob ;
        
        BLCA_add_to_component( comp_p, x, mixmod, 1 ) ;
        
        if( prior == 1 ) mixmod->component_compute = k;
        
        lpp = BLCA_return_log_marginal_likelihood_component( comp_p, mixmod );
        
        lpp_store[k] = lpp;
        
        probs[k] = lpp + lcm - ( lc + lp ) ;
        
        BLCA_add_to_component( comp_p, x, mixmod, -1 );
        
      }
      else
      {
        //if k==g this is easy
        probs[k] = 0.0;
      }
      
    }
    
    max = BLCA_get_max(probs,G);
    
    nrmcnst = 0.;
    
    for( k=0; k<G; k++  )
    {
      probs[k] -= max;
      probs[k] = exp(probs[k]);
      nrmcnst += probs[k];
    }
    
    for( k=0; k<G;  k++ ) probs[k] /= nrmcnst;
    
    
    /*sample allocation from the vector of weights*/
    g_prime = BLCA_sample_discrete(probs,G);
    
    
    if( g_prime != g )
    {
      
      comp_p = mixmod->components[ wis[g_prime] ] ;
      
      z[ii] = g_prime;
      
      BLCA_add_to_component( comp_p, x, mixmod, 1 );
      
      comp_p->log_prob = lpp_store[g_prime] ;
      
      comp_c->log_prob = lcm ;
      
    }
    else
    {
      BLCA_add_to_component( comp_c, x, mixmod, 1 );
    }
    
  }
  
  free(probs);
  free(order);
  free(lpp_store);
  
  return(TRUE);
}

/*double BLCA_get_log_relative_probability_for_gibbs(int *x,struct component *component_k,struct component *component_g,struct mix_mod *mixmod)
 //compute log of the relative probability for assigning entry where x is a pointer to the datapoint
 {
 int i,j; 
 double t1,t2,t3,t4,ll,I;
 
 
 t1 = 0.;
 t2 = 0.;
 
 for(j=0;j<mixmod->d;j++){
 
 if(mixmod->varindicator[j]){
 
 t1 += lgamma(mixmod->ncat[j]*mixmod->beta) - mixmod->ncat[j]*lgamma(mixmod->beta) -lgamma(component_k->n_g + 1. + mixmod->ncat[j]*mixmod->beta);
 t2 += lgamma(mixmod->ncat[j]*mixmod->beta) - mixmod->ncat[j]*lgamma(mixmod->beta) -lgamma(component_g->n_g - 1. + mixmod->ncat[j]*mixmod->beta);
 
 for(i=0;i<mixmod->ncat[j];i++){
 
 I = (x[j] == i) ? 1 : 0 ;
 
 t1 += lgamma(component_k->N[j][i] + I + mixmod->beta);
 t2 += lgamma(component_g->N[j][i] - I + mixmod->beta);
 
 }
 
 }
 
 }
 
 //add all up
 
 t3 = lgamma(component_k->n_g+1.+mixmod->alpha)+lgamma(component_g->n_g-1.+mixmod->alpha) 
 + t1 + t2;
 
 ll = t3 - component_k->log_prob - component_g->log_prob;
 
 return(ll);
 
 }*/

int BLCA_update_allocations_with_metropolis_move_1(struct mix_mod *mixmod,int *accepted,int *proposed)
  /*this performs the move M1 taken from Nobile and Fearnside (2007) Stats and Computing 17: p147-162*/
{
  
  //printf("\nwithin Move 1...\n");
  
  /*DO NOT PERFORM THIS MOVE UNLESS THERE ARE AT LEAST TWO COMPONENTS*/
  if(mixmod->G < 2) return(TRUE);
  
  *proposed += 1;
  
  int i, ii, kk, j, k, g1, g2, ig1, ig2, ntot, G=mixmod->G, d=mixmod->d, n=mixmod->n, 
    *indexes, *proposed_alloc, *x, *wis=mixmod->whereis, *z=mixmod->z, *y = mixmod->y;
  double p, log_acceptance, ag1=1.0, ag2=1.0, alpha;
  struct component *component_g1, *component_g2, *curr_component_g1, *curr_component_g2;
  
  //the integers g1 and g2 give the components, and ig1 and ig2 their whereis value
  //the doubles ag1 and ag2 give the parameters to the Beta distribution for generating p_1
  // take default uniform
  
  //sample the two components
  g1 = (int) ( runif(0.0,1.0) * G ) ;
  g2 = g1;
  while(g2 == g1) g2 =  (int) ( runif(0.0,1.0) * G ) ;
  
  //find where in mixmod->components g1 and g2 are
  ig1 = wis[g1];
  ig2 = wis[g2];
  
  curr_component_g1 = mixmod->components[ig1];
  curr_component_g2 = mixmod->components[ig2];
  
  ntot = curr_component_g1->n_g + curr_component_g2->n_g;
  
  //DO NOT PERFORM THIS MOVE UNLESS ntot > 0
  if(ntot == 0) return(TRUE);
  
  //allocate space for proposal stuff...
  
  component_g1 = BLCA_create_component(mixmod); //(struct component *)malloc(sizeof(struct component));
  component_g2 = BLCA_create_component(mixmod); //(struct component *)malloc(sizeof(struct component));
  
  //BLCA_allocate_component(component_g1,mixmod);
  //BLCA_allocate_component(component_g2,mixmod);
  
  /*allocate a vector to keep track of indexes*/
  indexes = (int *)calloc(ntot,sizeof(int));
  proposed_alloc = (int *)calloc(ntot,sizeof(int));
  
  k=0;
  for( i=0; i<n; i++ )
  {
    if( z[i] == g1 || z[i] == g2)
    {
      indexes[k] = i;
      k+=1;
    }
  }
  
  // generate p and begin reallocation
  p = rbeta(ag1,ag2);
  
  for( i=0; i<ntot; i++ )
  {
    
    ii = indexes[i];
    
    x = y + ii * d;  
    
    if( runif(0.0,1.0) < p )
    {
      //reallocate to g1
      proposed_alloc[i] = g1;
      BLCA_add_to_component( component_g1, x, mixmod, 1 );
    }else{
      //reallocate to g2
      proposed_alloc[i] = g2;
      BLCA_add_to_component( component_g2, x, mixmod, 1 );
    }		
  }
  
  //evaluate log of acceptance probability
  
  if( mixmod->prior_type == 1 )
  {
    mixmod->component_compute = 0; 
    alpha = mixmod->alpha_prior[0];
  }else{
    //alpha = mixmod->alpha;
  }
  //doesn't matter what this is set to as these moves are only for symmetric priors
  BLCA_recompute_marginal_likelihood_component( component_g1, mixmod );
  BLCA_recompute_marginal_likelihood_component( component_g2, mixmod );
  
  log_acceptance = component_g1->log_prob + component_g2->log_prob - curr_component_g1->log_prob - curr_component_g2->log_prob
    + lgamma( alpha + curr_component_g1->n_g) + lgamma( alpha + curr_component_g2->n_g) 
    - lgamma( alpha + component_g1->n_g) - lgamma( alpha + component_g2->n_g);
    
    Rprintf("\n m1 log acc %lf", log_acceptance);
    
    //printf("\nThe value of log acceptance is %.10f,",log_acceptance);
    
    if( log(runif(0.0,1.0)) < log_acceptance )
    {
      
      *accepted += 1;
      
      for( i=0; i<ntot; i++ ) z[ indexes[i] ] = proposed_alloc[i];
      
      //copy over the accepted components
      BLCA_copy_component(component_g1, curr_component_g1, mixmod);
      BLCA_copy_component(component_g2, curr_component_g2, mixmod);
      
    }
    
    
    BLCA_free_component(component_g1,mixmod);
    BLCA_free_component(component_g2,mixmod);
    
    free(component_g1);
    free(component_g2);
    
    free(indexes);
    free(proposed_alloc);
    
    return(TRUE);
}

// need to finish moves 2 and 3 and update the storage of Y for the collapsed and VB cases

int BLCA_update_allocations_with_metropolis_move_2(struct mix_mod *mixmod,int *accepted,int *proposed)
  //this performs the move M2 taken from Nobile and Fearnside (2007) Stats and Computing 17: p147-162
{
  
  /*DO NOT PERFORM THIS MOVE UNLESS THERE ARE AT LEAST TWO COMPONENTS*/
  if(mixmod->G < 2) return(TRUE);
  
  
  int i, ii, j, k, g1, g2, ig1, ig2, G=mixmod->G, d=mixmod->d, n=mixmod->n, curr_n_g1, curr_n_g2, m, ntot, c=0, 
    *indexes, *order, *z=mixmod->z, *y=mixmod->y, *wis=mixmod->whereis, *x;
  double log_acceptance;
  struct component *component_g1, *component_g2, *curr_component_g1, *curr_component_g2;
  
  /*sample the two components*/
  g1 =  (int) ( runif(0.0,1.0) * G );
  g2 = g1;
  while(g2 == g1) g2 =  (int) ( runif(0.0,1.0) * G );
  
  
  // find where in mixmod->components g1 and g2 are
  ig1 = wis[g1];
  ig2 = wis[g2];
  
  curr_component_g1 = mixmod->components[ ig1 ];
  curr_component_g2 = mixmod->components[ ig2 ];
  
  // check if move reasonable
  if(curr_component_g1->n_g == 0) return(TRUE);
  
  
  //current component sizes
  curr_n_g1 = curr_component_g1->n_g;
  curr_n_g2 = curr_component_g2->n_g;
  
  //allocate the candidate components
  
  component_g1 = BLCA_create_component( mixmod );
  component_g2 = BLCA_create_component( mixmod );
  
  /*component_g1 =(struct component *)malloc(sizeof(struct component));
   component_g2 =(struct component *)malloc(sizeof(struct component));
   
   BLCA_allocate_component(component_g1,mixmod);
   BLCA_allocate_component(component_g2,mixmod);*/
  
  *proposed += 1;
  
  order = (int *)calloc(curr_n_g1,sizeof(int));
  for(i=0;i<curr_n_g1;i++) order[i] = i;
  
  //shuffle the order
  BLCA_random_ranshuffle(order,curr_n_g1);
  
  indexes = (int *)calloc(curr_n_g1,sizeof(int));
  for( i=0; i<n; i++ )
  {
    if( z[i] == g1 )
    {
      indexes[c] = i;
      c += 1;
    }
  }
  
  m = (int)( runif(0.0,1.0) * curr_n_g1 ) ;  //gsl_rng_uniform_int(r,curr_n_g1);
  
  BLCA_copy_component( curr_component_g1, component_g1, mixmod );
  BLCA_copy_component( curr_component_g2, component_g2, mixmod );
  
  // component_g1->n_g -= m;
  //  component_g2->n_g += m;
  
  for( i=0; i<m; i++ )
  {
    
    k = order[i];
    ii = indexes[k];
    
    x = y + ii * d;
    
    BLCA_add_to_component( component_g1, x, mixmod, -1 );
    BLCA_add_to_component( component_g2, x, mixmod, 1 );
    
  }
  
  //compute the log acceptance probability
  
  if( mixmod->prior_type == 1 ) mixmod->component_compute =  0;
  
  BLCA_recompute_marginal_likelihood_component( component_g1, mixmod );
  BLCA_recompute_marginal_likelihood_component( component_g2, mixmod );
  
  log_acceptance = component_g1->log_prob + component_g2->log_prob - curr_component_g1->log_prob - curr_component_g2->log_prob
    + log(component_g1->n_g + m) - log(component_g2->n_g) + lgamma(component_g1->n_g + m + 1.0) + lgamma( component_g2->n_g - m + 1.0)
    -lgamma(component_g1->n_g + 1.0) - lgamma(component_g2->n_g + 1.0);
    
    if( log( runif(0.0,1.0) ) < log_acceptance )
    {
      //do the swap
      *accepted += 1;
      
      for( i=0; i<m; i++ )
      {
        k = order[i];
        ii = indexes[k];
        z[ii] = g2;		
      }
      
      //copy over the accepted components
      BLCA_copy_component( component_g1, curr_component_g1, mixmod );
      BLCA_copy_component( component_g2, curr_component_g2, mixmod );
    }
    
    
    BLCA_free_component(component_g1,mixmod);
    BLCA_free_component(component_g2,mixmod);
    
    free(component_g1);
    free(component_g2);
    
    free(order);
    free(indexes);
    
    return(TRUE);
}



int BLCA_update_allocations_with_metropolis_move_3(struct mix_mod *mixmod,int *accepted,int *proposed)
  // this performs the move M3 taken from Nobile and Fearnside (2007) Stats and Computing 17: p147-162
{
  
  // DO NOT PERFORM THIS MOVE UNLESS THERE ARE AT LEAST TWO COMPONENTS
  if( mixmod->G < 2) return(TRUE);
  
  int i, ii, j, k, g1, g2, ig1, ig2, curr_n_g1, curr_n_g2, m, ntot, c=0, G=mixmod->G, n=mixmod->n, d=mixmod->d, identify_g1, identify_g2, id;
  int *indexes, *order, *proposed_allocation, *z=mixmod->z, *wis=mixmod->whereis, *y=mixmod->y, *x ;
  double w, log_acceptance, log_transition_z_to_zprime=0., log_transition_zprime_to_z=0., 
    l1, l2, lp_curr_1, lp_curr_2, lp_prop_1, lp_prop_2, p1, p2, max, alpha = mixmod->alpha_prior[0], eps=1E-8;
  struct component *component_g1, *component_g2, *curr_component_g1, *curr_component_g2;
  
  *proposed += 1;
  
  // sample the two components
  g1 =  (int) ( runif(0.0,1.0) * G );
  g2 = g1;
  while(g2 == g1) g2 = (int) ( runif(0.0,1.0) * G );
  
  
  //find where in mixmod->components g1 and g2 are
  ig1 = wis[g1];
  ig2 = wis[g2];
  
  curr_component_g1 =  mixmod->components[ig1];
  curr_component_g2 = mixmod->components[ig2];
  
  ntot = curr_component_g1->n_g + curr_component_g2->n_g;
  
  // only attempt if possible
  if(ntot == 0) return(TRUE);
  
  indexes = (int *)calloc( ntot, sizeof(int) );
  order = (int *)calloc( ntot, sizeof(int) );
  proposed_allocation = (int *)calloc( ntot, sizeof(int) );
  
  //this move can still be done if either component empty
  for( i=0; i<n; i++ )
  {	
    if( z[i] == g1 ||  z[i] == g2)
    {
      indexes[c] = i;
      c += 1;
    }
  }
  
  for( i=0; i<ntot; i++ ) order[i] = i;
  
  //shuffle the order
  BLCA_random_ranshuffle( order, ntot );
  
  //allocate the candidate components
  component_g1 = BLCA_create_component( mixmod );
  component_g2 = BLCA_create_component( mixmod );
  
  //randomizing this part should give a 50% acceptance
  
  k = order[0];
  ii = indexes[k];
  x = y + ii * d; 
  
  if( z[ii] == g1 )
  {
    //put it into component_g1
    BLCA_add_to_component( component_g1, x, mixmod, 1  );
    proposed_allocation[0] = g1;
    log_transition_z_to_zprime = log(0.5) - log(1.);
  }
  else
  {
    //put it into component_g2
    BLCA_add_to_component( component_g2, x, mixmod, 1 );
    proposed_allocation[0] = g2;
    log_transition_z_to_zprime = log(0.5) - log(1.);
  }
  
  /*identify this as being in the same component as original
   ...this allows the computing of the backwards probability*/
  
  log_transition_zprime_to_z += log(0.5) - log(1.);
  
  for( i=1; i<ntot; i++ )
  {
    
    k = order[i];
    ii = indexes[k];
    x = y + ii * d; 
    
    // compute probability generated from g1 etc..
    
    l1 = BLCA_compute_log_data_probability_with_inclusion_in_component( x , component_g1, mixmod )
      + BLCA_compute_log_data_probability_component( component_g2, mixmod );
    
    l2 = BLCA_compute_log_data_probability_component( component_g1, mixmod )
      + BLCA_compute_log_data_probability_with_inclusion_in_component( x, component_g2, mixmod );
    
    w = ( (alpha + component_g1->n_g) / (alpha + component_g2->n_g) ) * exp( l1 - l2 );
    
    p1 = w/(1.0+w);
    
    // make a draw
    if( runif(0.0,1.0) < p1 )
    {
      // put it in g1
      BLCA_add_to_component( component_g1, x, mixmod, 1 );
      if( z[ii] != g1 )
      {
        log_transition_z_to_zprime += log(p1);
        if( 1.0 - p1 < eps )
          log_transition_zprime_to_z += log( eps );
        else  
          log_transition_zprime_to_z += log(1.-p1);
      }
      proposed_allocation[i] = g1;
    }
    else
    {
      //put it in g2
      BLCA_add_to_component( component_g2, x, mixmod, 1 );
      if( z[ii] != g2 )
      {
        log_transition_z_to_zprime += log(1.-p1);
        if( p1 < eps )
          log_transition_zprime_to_z += log(eps);
        else
          log_transition_zprime_to_z += log(p1);
      }
      proposed_allocation[i] = g2;
    }
    
  }
  
  
  //compute the acceptance probability
  BLCA_recompute_marginal_likelihood_component(component_g1,mixmod);
  BLCA_recompute_marginal_likelihood_component(component_g2,mixmod);
  
  log_acceptance = component_g1->log_prob + component_g2->log_prob 
    - curr_component_g1->log_prob - curr_component_g2->log_prob
    + log_transition_zprime_to_z - log_transition_z_to_zprime;
    
    Rprintf("\n m3 log acc %lf", log_acceptance);
    
    if( log( runif(0.0,1.0) ) < log_acceptance )
    {
      *accepted += 1;
      //accept the move and update all quantities
      //allocations first
      BLCA_copy_component(component_g1, curr_component_g1, mixmod);
      BLCA_copy_component(component_g2, curr_component_g2, mixmod);
      
      for( i=0; i<ntot; i++ )
      {
        k = order[i];
        ii = indexes[k];
        z[ii] = proposed_allocation[i];
      }
      
    }
    
    //free up all memory
    BLCA_free_component(component_g1,mixmod);
    BLCA_free_component(component_g2,mixmod);
    
    free(component_g1);
    free(component_g2);
    
    free(indexes);
    free(order);
    free(proposed_allocation);
    
    return(TRUE);
}


