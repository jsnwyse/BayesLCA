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

#include "BLCA_eject_absorb.h"

/*eject and absorb moves*/

int BLCA_update_allocations_with_ejection_move(struct mix_mod *mixmod,int *accepted,int *proposed,double pr_ej_G,double pr_ej_Gp1)
  /*this is the ejection move for one comonent ejecting another*/
{
  
  int i, ii, j, k, g1, g2, ig1, ig2, ntot, c=0, d = mixmod->d,
  *indexes, *order, *proposed_allocation, G = mixmod->G, n = mixmod->n, *x, *z = mixmod->z, *wis = mixmod->whereis, *ncat = mixmod->ncat, *y = mixmod->y;
  double w, a, prob_put_in_g2, log_acceptance,  log_transition_z_to_zprime=0.,  log_transition_zprime_to_z=0.;
  struct component *component_g1, *component_g2, *ej_component;
  
  // do a check for the constraint
  
  //int S=0,D=0;
  //double logP=0.,cstr;
  
  /*for(j=0;j<mixmod->d;j++){
   if(mixmod->varindicator[j]){
   S += mixmod->ncat[j];
   logP += log(mixmod->ncat[j]);
   D += 1;
   }
  }*/
  
  /*check constraint*/
  /*cstr = logP - log((S - D + 1)*(mixmod->G+1));
   if(!(cstr>0.)){
   //model will not be identifiable
   // Rprintf("\nConstraint not satisfied");
   return -1;
   }*/
  
  
  *proposed += 1;
  
  //sample the ejecting component
  g1 = (int) ( runif(0.0,1.0) * G );
  g2 = G;
  
  //find where in mixmod->components g1 is
  ig1 = wis[g1];
  ej_component = mixmod->components[ ig1 ];
  
  //if this component is empty we need a special case
  
  component_g1 = BLCA_create_component( mixmod );
  component_g2 = BLCA_create_component( mixmod );
  
  /*component_g1 =(struct component *)malloc(sizeof(struct component));
   component_g2 =(struct component *)malloc(sizeof(struct component));
   
   BLCA_allocate_component(component_g1,mixmod);
   BLCA_allocate_component(component_g2,mixmod);
   
   component_g1->n_g = 0;
   component_g2->n_g = 0;*/
  
  ntot = ej_component->n_g;	
  
  if(ntot > 0)
  { 
    //this is the case for ejecting from a non-empty component
    indexes = (int *)calloc( ntot, sizeof(int) );
    proposed_allocation = (int *)calloc( ntot, sizeof(int) );
    
    c = 0;
    for( i=0; i<n; i++ )
    { 
      if( z[i] == g1 )
      { 
        indexes[c] = i; 
        c += 1;
      } 
    }	 
    
    //copy the contents 
    BLCA_copy_component( ej_component, component_g1, mixmod );	
    
    // generate the probability of assignment to the new component
    if(ntot < 4)
    { 
      //just set a = 100
      a = 100.;
      prob_put_in_g2 = rbeta(a,a);
    } 
    else
    { 
      a = a_table[ntot-1];
      prob_put_in_g2 = rbeta(a,a);
    }  
    
    //now reassign or not
    for( i=0; i<ntot; i++ )
    { 
      ii = indexes[i];
      
      x = y + ii * d; 
      
      if( runif(0.0,1.0) < prob_put_in_g2)
      { 
        //then move this point to g2
        BLCA_add_to_component( component_g1, x, mixmod, -1 );
        BLCA_add_to_component( component_g2, x, mixmod, 1 );
        proposed_allocation[i] = g2;
      }else{ 
        proposed_allocation[i] = g1;
      }
    } 
  } 
  
  mixmod->component_compute = 0 ; //in this case the priors are symmetric
  BLCA_recompute_marginal_likelihood_component( component_g1, mixmod );
  BLCA_recompute_marginal_likelihood_component( component_g2, mixmod );	
  
  //compute the acceptance probability
  w = BLCA_log_normalizing_constant_model( G+1, mixmod ) - BLCA_log_normalizing_constant_model( G, mixmod );
  
  log_transition_z_to_zprime = log(pr_ej_G);
  
  if(ntot > 0)
  { 
    log_transition_z_to_zprime +=  lgamma(2.*a) - 2.*lgamma(a) + lgamma(a + component_g1->n_g) 
    + lgamma(a + component_g2->n_g) - lgamma(2.*a + ntot);
  } 
  
  log_transition_zprime_to_z = log(1.-pr_ej_Gp1);
  
  log_acceptance = w + component_g1->log_prob + component_g2->log_prob - ej_component->log_prob 
    - log_transition_z_to_zprime + log_transition_zprime_to_z + mixmod->log_prior_G[G+1] - mixmod->log_prior_G[G];
  
  if(log( runif(0.0,1.0) ) < log_acceptance)
  {  
    *accepted += 1;
    
    //update the model structure
    mixmod->G += 1;
    
    //relabel the appropriate indexes
    if(ntot>0)
    { 
      for( i=0; i<ntot; i++ )
      { 
        ii = indexes[i];
        z[ii] = proposed_allocation[i];
      }  
    }
    
    int new_whereis=0;
    
    //begin to encode an unused component by -1 in mixmod->whereis
    while(mixmod->components[new_whereis]->in_use == TRUE) new_whereis += 1;
     
    // put new componenet G in new_whereis
    // this will now be indexed as component G
    wis[G] = new_whereis;
    mixmod->components[new_whereis]->in_use = TRUE;
    
    BLCA_copy_component( component_g1, ej_component, mixmod );
    BLCA_copy_component( component_g2, mixmod->components[new_whereis], mixmod );
    
    //now do a swap between component G and one of the others...
    
    //generate component randomly
    g1 = (int)( runif(0.0,1.0) * ( G + 1 ) ) ;
    
    if( g1 != G )
    {
      //do a swap!
      ig1 = wis[g1];
      ig2 = wis[G];
      
      wis[g1] = ig2;
      wis[G] = ig1;
       
      //relabel the appropriate indexes
      for( i=0; i<n; i++ )
      { 
        if( z[i] == g1 ){
          z[i] = G;
        }else if( z[i] == G){ 
          z[i] = g1;
        }
      }
    }
  }
  
  BLCA_free_component( component_g1, mixmod );
  BLCA_free_component( component_g2, mixmod );
  
  free(component_g1);
  free(component_g2);
  
  if( ntot>0 )
  {
    free(indexes);
    free(proposed_allocation);
  } 
  
  return(TRUE);
}


int BLCA_update_allocations_with_absorb_move(struct mix_mod *mixmod,int *accepted,int *proposed,double pr_ej_G,double pr_ej_Gm1)
{ 
  int i, ii, j, k, g1, g2, ig1, ig2, ntot, n=mixmod->n, d = mixmod->d, G = mixmod->G;
  int n_g2, *wis = mixmod->whereis, *z=mixmod->z;
  double w, a, log_acceptance, log_transition_z_to_zprime=0.,log_transition_zprime_to_z=0.;
  struct component *component_g1, *component_g2, *ab_component;
  
  /*do a check for the constraint*/
  
  /*int S=0,D=0;
  double logP=0.,cstr;
  
  for(j=0;j<mixmod->d;j++){
    if(mixmod->varindicator[j]){
      S += mixmod->ncat[j];
      logP += log(mixmod->ncat[j]);
      D += 1;
    }
  }*/
  
  //check constraint
  /*cstr = logP - log((S - D + 1)*(mixmod->G-1));
  if(!(cstr>0)){
    //model will not be identifiable
    // Rprintf("\nConstraint not satisfied");
    return -1;
  }*/
  
  
  *proposed += 1;
  
  // choose component to absorb into and to absorb
  g1 =  (int) ( runif(0.0,1.0) * G ) ;
  g2 = g1;
  while(g2 == g1) g2 = (int) ( runif(0.0,1.0) * G ) ;
  
  
  // find where in mixmod->components g1 and g2 are
  ig1 = wis[g1];
  ig2 = wis[g2];
  
  component_g1 = mixmod->components[ig1];
  component_g2 = mixmod->components[ig2];
  
  //use a component to store the proposed
  ab_component = BLCA_create_component( mixmod );
  
  // form the proposed component by combining the two other components...
  //ntot = component_g1->n_g + component_g2->n_g;
  
  // copy first component and add second
  BLCA_copy_component( component_g1, ab_component, mixmod );
  BLCA_add_components( component_g2, ab_component, mixmod );
  
  /*n_g2 = component_g2->n_g;
  if(n_g2 > 0){
    
    //this is the case for non-empty components
    indexes = calloc(n_g2,sizeof(int));
    proposed_allocation = calloc(n_g2,sizeof(int));
    
    c = 0;
    
    for(i=0;i<mixmod->n;i++){
      if(mixmod->z[i] == g2){
        indexes[c] = i;
        c += 1;
      }
    }
    
    //now just place everything in component_g1
    
    for(i=0;i<n_g2;i++){
      
      ii = indexes[i];
      
      x = mixmod->Yobs[ii];
      
      BLCA_add_to_component( component_g1, x, mixmod, 1 );
      
    }	
    
  }*/
  
  mixmod->component_compute = 0;
  BLCA_recompute_marginal_likelihood_component( ab_component, mixmod );	
  
  /*compute the acceptance probability, remembering to add all necessary normalizing constants*/
  
  // w = log of difference in normalizing constants
  
  w = BLCA_log_normalizing_constant_model( G-1, mixmod ) - BLCA_log_normalizing_constant_model( G, mixmod ); //need to check this!
  
  //printf("\nDeath: w = %.10f ",w);
  
  log_transition_zprime_to_z = log( pr_ej_Gm1 );
  
  if(ab_component->n_g > 0)
  {
    if( ab_component->n_g < 4) a = 100.; else a = a_table[ ab_component->n_g - 1 ];
    log_transition_zprime_to_z +=  lgamma(2.*a) - 2.*lgamma(a) + lgamma(a + component_g1->n_g) + lgamma(a + component_g2->n_g) - lgamma(2.*a + ab_component->n_g);
  }
  
  log_transition_z_to_zprime = log(1.-pr_ej_G);	
  
  log_acceptance = w + ab_component->log_prob - component_g1->log_prob - component_g2->log_prob
    - log_transition_z_to_zprime + log_transition_zprime_to_z + mixmod->log_prior_G[G-1] - mixmod->log_prior_G[G];	
  
  if(log( runif(0.0,1.0) ) < log_acceptance)
    {
    
    *accepted += 1;
    
    //update the model structure
    mixmod->G -= 1;
    
    //relabel the appropriate indexes
    for( i=0; i<n; i++ ) 
    {
      if( z[i] == g2 ) z[i] = g1;
    }

    BLCA_copy_component( ab_component, component_g1, mixmod);
    
    component_g2->in_use = FALSE;

    // should relabel everything from component g2 upwards
    
    for( k=g2+1; k<G; k++ )
    {
      for( i=0; i<n; i++ )
      {
        if( z[i] == k )  z[i] = k-1;
      }
      
      j = wis[k];
      wis[k-1] = j;
    }
    
    wis[G-1] = -1;
  }
  
  BLCA_free_component( ab_component, mixmod );
  free( ab_component );

  return(TRUE);
}

