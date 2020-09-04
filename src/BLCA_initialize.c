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
 
 Last modification of this code: Tue 08 Aug 2017 21:04:52 IST   */


#include "BLCA_initialize.h"

void BLCA_set_prior_on_number_of_components(struct mix_mod *mixmod,int type)
{
  
  int i;
  
  switch(type)
  {
    
  case RICHARDSON_AND_GREEN:
    
    for(i=0;i<mixmod->maxgroups+1;i++){
      mixmod->log_prior_G[i] = 0.0;
    }
    
    break;
    
  case NOBILE_AND_FEARNSIDE:
    
    for(i=1;i<mixmod->maxgroups+1;i++){
      mixmod->log_prior_G[i] = -lgamma((double)i + 1.0);
    }
    
    break;
    
  }
  
  
  return;
  
}

void BLCA_initialize_data( struct mix_mod *mixmod, int *Y )
{
  int i, j, x, n = mixmod->n, d = mixmod->d ;
  
  for(i=0;i<n;i++)
  {
    for(j=0;j<d;j++)
    {
      x = Y[j*n + i];
      mixmod->Y[j][i] = x;
      mixmod->Yobs[i][j] = x;
    }
  }
}


void BLCA_initialize_priors( struct mix_mod *mixmod, double *alpha_prior, double *beta_prior, int type )
{
  int k, j, c, p, G = mixmod->maxgroups, d = mixmod->d, *ncat = mixmod->ncat ;
  double *prior;
  
  // three different types- type 0 is the simple initialization
  if( type == 0 )
  {
    // one parameter for each; these values have already been set above
    //mixmod->alpha = alpha_prior[0];
    //mixmod->beta = beta_prior[0];
  }
  
  if( type == 1 )
  {
    //all vectors
    //all values take the same
    mixmod->lg_alpha_sum = 0.0;
    mixmod->lg_sum_alpha = 0.0;
    for( k=0; k<G; k++ ) 
    {
      mixmod->alpha_prior[k] = alpha_prior[k];
      mixmod->lg_alpha_sum += lgamma( alpha_prior[k] );
      mixmod->lg_sum_alpha += alpha_prior[k];
    }
    mixmod->lg_sum_alpha = lgamma( mixmod->lg_sum_alpha );
    
    p = 0;
    for( j=0; j<d; j++ )
    {
      for( k=0; k<G; k++ )
      {
        prior = mixmod->beta_prior[k][j];
        for( c=0; c<ncat[j]; c++ ) 
        {
          prior[c] = beta_prior[ p + c ];
          mixmod->lg_beta_sum[k][j] += lgamma(prior[c]);
          mixmod->lg_sum_beta[k][j] += prior[c];
        }
        mixmod->lg_sum_beta[k][j] = lgamma( mixmod->lg_sum_beta[k][j] );
      }
      p += ncat[j] ;
    }	
  }
  
  if( type == 2 )
  {
    //all take different values
    mixmod->lg_alpha_sum = 0.0;
    mixmod->lg_sum_alpha = 0.0;
    for( k=0; k<G; k++ ) 
    {
      mixmod->alpha_prior[k] = alpha_prior[k]; 
      mixmod->lg_alpha_sum += lgamma( alpha_prior[k] );
      mixmod->lg_sum_alpha += alpha_prior[k];
    }
    mixmod->lg_sum_alpha = lgamma( mixmod->lg_sum_alpha );
    
    p = 0;
    for( j=0;  j<d; j++ )
    {
      for( k=0; k<G; k++ )
      {
        prior = mixmod->beta_prior[k][j];
        for( c=0; c<mixmod->ncat[j]; c++ ) 
        {
          prior[c] = beta_prior[ p*G + k*ncat[j] + c ]; 
          mixmod->lg_beta_sum[k][j] += lgamma(prior[c]);
          mixmod->lg_sum_beta[k][j] += prior[c];
        }
        mixmod->lg_sum_beta[k][j] = lgamma( mixmod->lg_sum_beta[k][j] );
      }
      p += ncat[j] ;
    }	
  }
  
  //trying to fix problem with eject absorb moves when prior_type==1
  mixmod->prior_type = 1;
  
  return;
  
}

void BLCA_initialize_Gibbs_sampler( int init_type, struct mix_mod *mixmod  )
{
  int g, k, j, c, p=0, gap_, n = mixmod->n, G = mixmod->G, d = mixmod->d, 
    *vind = mixmod->varindicator, *z = mixmod->z, *y = mixmod->y, *ncat = mixmod->ncat;
  double s=0.0, *prior, *prob, *sk, *w = mixmod->weights;
  struct component *comp;
  
  // go through three possible initialization types making the s 
  // matrix in each case
  
  // option 'single' and membership initialization for init_type == 3
  
  if( mixmod->collapsed || init_type == 0 || init_type == 2 )
  {
    for( k=0; k<n; k++ )
    {
      // generate a random group 0:G-1
      g = (int)( G * runif(0.0, 1.0) ) ;
      z[k] = g; 
      comp = mixmod->components[g];
      // this is where the problem is with the collapsed sampler... 
      BLCA_add_to_component( comp, y , mixmod, 1 ); 
      // backward compatibility
      if( !mixmod->collapsed ) mixmod->s[k][g] = 1.0;
      y += d;
    }	
  }
  
  // option 'across'
  if( init_type == 1 )
  {
    // generate a portion of a random group
    for( k=0; k<n; k++ )
    {
      sk = mixmod->s[k];
      s = 0.;
      for( g=0; g<G; g++ ) 
      {
        sk[g] = rgamma( 1.0, 1.0);
        s += sk[g];
      }
      for( g=0; g<G; g++ ) sk[g] /= s; 
    }			
  }
  
  
  if( init_type == 0 || init_type == 1 )
  {
    // class weights
    for( k=0; k<n; k++ )
    {
      sk = mixmod->s[k];
      for( g=0; g<G; g++ ) w[g] += sk[g];
    }
    // item probabilities
    y = mixmod->y;
    for( k=0; k<n; k++ )
    {
      for( g=0; g<G; g++ )
      {
        comp = mixmod->components[g];
        for( j=0; j<d; j++ )
        {
          prob = comp->prob_variables[j];
          prob[ y[j] ] += mixmod->s[k][g];
          for( c=0; c<ncat[j]; c++ ) prob[c] /= w[g];
        }
      }
      y += d;
    }
    // get initial item probs
    for( g=0; g<G; g++ )
    {
      comp = mixmod->components[g];
      for( j=0; j<d; j++ )
      {
        if( vind[j] )
        {
          y = mixmod->y;
          prob = comp->prob_variables[j];
          for( k=0; k<n; k++ ){ prob[ y[j] ] += mixmod->s[k][g]; y+=d; }
          for( c=0; c<mixmod->ncat[j]; c++ ) prob[c] /= w[g];
        }
      }
    }
    // normalize weights
    for( g=0; g<G; g++ ) w[g] /= n; 
  }
  
  // option 'prior'
  
  if( init_type == 2 )
  {
    //generate the probabilities from the prior
    
    //component weights
    s = 0.0;
    for( g=0; g<G; g++)
    {
      w[g] = rgamma( mixmod->alpha_prior[g] , 1.0 ) ;
      s += w[g];
    }
    for( g=0; g<G; g++) w[g] /= s;
    
    //variable probabilities
    for( g=0; g<G; g++)
    {
      comp = mixmod->components[g];
      for( j=0; j<d; j++)
      {
        if( vind[j] )
        {
          //only sample if the variable is included in the model
          s = 0.0;
          prob = comp->prob_variables[j];
          for( c=0; c<ncat[j]; c++ )
          {
            prob[c] =  rgamma( mixmod->beta_prior[g][j][c] , 1.0 ) ;
            s += prob[c];
          }
          for( c=0; c<ncat[j]; c++ ) prob[c] /= s;
        }
      }
    }
    
  }
  
  if( mixmod->collapsed )
  {
    mixmod->component_compute = 0;
    
    //compute log marginal likelihood for each component
    for( g=0; g<G; g++ )
      BLCA_recompute_marginal_likelihood_component( mixmod->components[g], mixmod );
    
    //fill up the part for the undiscriminating variables
    y = mixmod->y;
    for( k=0; k<mixmod->n; k++ ) 
    {
      BLCA_add_to_component( mixmod->undiscriminating, y, mixmod, 1 ); 
      y += d;
    }
    
    /*for( j=0; j<mixmod->d; j++ )
     {
     for( k=0; k<mixmod->n; k++ )
     mixmod->undiscriminating->N[j][ mixmod->Y[j][k] ] += 1;
     } */
    
    BLCA_recompute_marginal_likelihood_undiscriminating_variables( mixmod->undiscriminating, mixmod );
  }
  
  for( g=0; g<mixmod->G; g++ )
    mixmod->whereis[g] = g;
  
  return;
}


void BLCA_initialize_EM( int init_type, double *init_vals, struct mix_mod *mixmod , double *group_weights, double *prob_variables)
{
  int g, k, j, c, p=0, gap_, n = mixmod->n, G = mixmod->G, d = mixmod->d, *vind = mixmod->varindicator, *ncat = mixmod->ncat;
  double s=0.0, *prob, *sk, *w = mixmod->weights;
  struct component *comp;
  
  //initialize the membership matrix using either across/single or prespecified
  
  if( init_type == 0 )
  {
    for( k=0; k<n; k++ )
    {
      //generate random group 0:G-1
      g = (int)( G * runif( 0.0,1.0) ) ;
      mixmod->s[k][g] = 1.;
    }
  }
  
  if( init_type == 1 )
  {
    for( k=0; k<n; k++ )
    {
      s = 0.;
      sk = mixmod->s[k];
      for( g=0; g<G; g++ ) 
      {
        sk[g] = runif( 0.0, 1.0);
        s += sk[g];
      }
      for( g=0; g<G; g++ ) sk[g] /= s; 
    }
  }
  
  if( init_type == 2 )
  {
    for( k=0; k<n; k++ )
    {
      sk = mixmod->s[k];
      for( g=0; g<G; g++ ) sk[g] = init_vals[k + g*n ] ; 
    }
    
  }
  
  // randomly initialise the item probabilities
  //component weights
  s = 0.0;
  for( g=0; g<G; g++ )
  {
    w[g] = rgamma(1.0, 1.0) ;
    s += w[g];
  }
  for( g=0; g<G; g++ ) w[g] /= s;
  
  //variable probabilities
  for( g=0; g<G; g++ )
  {
    comp = mixmod->components[g];
    for( j=0; j<d; j++ )
    {
      if( vind[j] )
      {
        //only sample if the variable is included in the model
        s = 0.0;
        prob = comp->prob_variables[j];
        for( c=0; c<ncat[j]; c++ )
        {
          prob[c] =  rgamma(1.0, 1.0) ;
          s += prob[c];
        }
        for( c=0; c<ncat[j]; c++ ) prob[c] /= s;
      }
    }
  }
  
  return;
}


void BLCA_initialize_VB( int init_type, double *init_vals, struct mix_mod *mixmod , double *alpha_ud, double *beta_ud )
{
  int k, g, j, c, p=0, gap_, G = mixmod->G, n = mixmod->n, d = mixmod->d, *ncat = mixmod->ncat, *vind = mixmod->varindicator;
  double s, *sk;
  struct component *comp;
  
  for( k=0; k<n; k++ )
  {
    s = 0.0;
    sk = mixmod->s[k];
    for( g=0; g<G; g++ )
    {
      sk[g] = rgamma(1.0,1.0);
      s += sk[g];
    }
    for( g=0; g<G; g++ ) sk[g] /= s;
  }
  
  
  //initialize everything randomly
  /*mixmod->di_sum_alpha_ud = 0.0;
  for( g=0; g<G; g++ ) 
  {
    mixmod->alpha_ud[g] = 1.0 * runif( 0.0 , 1.0 ) ;
    mixmod->di_sum_alpha_ud += mixmod->alpha_ud[g];
    mixmod->di_alpha_ud[g] = digammaRN( mixmod->alpha_ud[g] );
  }
  mixmod->di_sum_alpha_ud = digammaRN( mixmod->di_sum_alpha_ud );
  
  for( g=0; g<G; g++ )
  {
    comp = mixmod->components[g] ;
    for( j=0; j<d; j++ )
    {
      if( vind[j] )
      {
        comp->di_sum_beta_ud[j] = 0.0;
        for( c=0; c<ncat[j]; c++ )
        {
          comp->beta_ud[j][c]  = 1.0 * runif( 0.0 , 1.0 );
          comp->di_sum_beta_ud[j] += comp->beta_ud[j][c] ;
          comp->di_beta_ud[j][c] = digammaRN( comp->beta_ud[j][c] );
        }
        comp->di_sum_beta_ud[j] = digammaRN( comp->di_sum_beta_ud[j] ) ;
      }
    }
  }*/
  return;
}

void BLCA_initialize_Boot( struct mix_mod *mixmod, double *group_weights, double *prob_variables  )
{
  
  int j, g, c, p, gap_, G = mixmod->G, n = mixmod->n, d = mixmod->d, *ncat = mixmod->ncat;
  double *prob;
  
  for( g=0; g<G; g++ ) mixmod->weights[g] = group_weights[g]; 
  
  p = 0;
  
  for( j=0; j<d; j++)
  {
    gap_ = p * G;
    for( g=0; g<G; g++)
    {
      prob = mixmod->components[g]->prob_variables[j] ;
      for( c=0; c<ncat[j]; c++ ) prob[c] = prob_variables[ gap_ + g*ncat[j] + c  ] ;
    }
    p += ncat[j];
  }
  return;
}



