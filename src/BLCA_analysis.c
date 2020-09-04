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
 
 Last modification of this code: Thu 04 Aug 2016 21:02:24 IST */

#include "BLCA_analysis.h"

struct results *BLCA_analysis_MCMC_collapsed( struct mix_mod *mixmod, int num_iteration, int num_burnin, int thin_by, int fixedG, int onlyGibbs, 
                                              int selectVariables, int *group_memberships, int *variable_inclusion_indicator, int *n_groups, double *log_joint_posterior, 
                                              double *prior_include, int *var_pattern, int verbose, int verbose_update )
  /*fixedG takes the value either TRUE or FALSE as defined in the macros*/
{
  
  int i, j, k, l, itmod, d_in, v_in, ej_case, vs_case, maxgroups = mixmod->maxgroups, n = mixmod->n, d = mixmod->d, 
    *ncat = mixmod->ncat, *vind = mixmod->varindicator, *z = mixmod->z ;
  
  int gap = (int)( num_iteration / thin_by) ;
  
  //writeToFile = FALSE;
  
  //struct results res,*results;
  //results = &res;
  
  
  double pr_ej_G,pr_ej_Gm1,pr_ej_Gp1;
  
  struct results *results = (struct results *)malloc(sizeof(struct results));
  BLCA_allocate_results( results, num_iteration, num_burnin, thin_by, n, d  );
  
  /*initialize with all variables in model*/
  //pass in the pattern from R
  for( j=0; j<d; j++ ) vind[j] = var_pattern[j];
  
  
  if( verbose ) Rprintf("\nInitializing sampler... starting burn-in...");
  
  for( l=0; l<num_iteration + num_burnin; l++ )
  {
    
    R_CheckUserInterrupt();
    
    /*check for violation of identifiability constraint*/
    /*int s = 0;
     double log_p = 0;
     int v_in = 0;
     
     for( j=0; j<d; j++ )
     {
     
     if( vind[j] )
     {
     s += ncat[j];
     log_p += log( ncat[j] );
     v_in += 1;
     }
     }*/
    
    //if(!(log_p > log((s - v_in + 1)*mixmod->G)) && !(l>0)){
    //Rprintf("\n***Warning***: identifiability constraint not satisfied; log_p = %d, s=%d, in=%d, G=%d",log_p,s,in,mixmod->G);
    //}
    
    BLCA_update_allocations_with_gibbs(mixmod);
    
    if(!onlyGibbs)
    {
      BLCA_update_allocations_with_metropolis_move_1(mixmod,&(results->accepted_m1),&(results->proposed_m1));
      BLCA_update_allocations_with_metropolis_move_2(mixmod,&(results->accepted_m2),&(results->proposed_m2));
      BLCA_update_allocations_with_metropolis_move_3(mixmod,&(results->accepted_m3),&(results->proposed_m3));
    }
    
    if(!fixedG)
    {
      
      if(mixmod->G == 1){
        ej_case = 0;
      }else if(mixmod->G == maxgroups){
        ej_case = 1;
      }else if(mixmod->G == 2){
        ej_case = 2;
      }else if(mixmod->G == maxgroups-1){
        ej_case = 3;
      }else{
        ej_case = 4;
      }
      
      switch(ej_case)
      {
      case 0:
        pr_ej_G = 1.; pr_ej_Gp1 = .5; pr_ej_Gm1 = 0.;
        break;
        
      case 1:
        pr_ej_G = 0.; pr_ej_Gp1 = 0.; pr_ej_Gm1 = .5;
        break;
        
      case 2:
        pr_ej_G = .5; pr_ej_Gp1 = .5; pr_ej_Gm1 = 1.;
        break;
        
      case 3:
        pr_ej_G = 0.5; pr_ej_Gp1 = 0.; pr_ej_Gm1 = 0.5;
        break;
        
      case 4:
        pr_ej_G = 0.5; pr_ej_Gp1 = 0.5; pr_ej_Gm1 = 0.5;
        break;
      }
      
      if( runif(0.0,1.0) < pr_ej_G)
      {		
        BLCA_update_allocations_with_ejection_move(mixmod,&(results->accepted_eject),&(results->proposed_eject),pr_ej_G,pr_ej_Gp1);
      } 
      else
      {
        BLCA_update_allocations_with_absorb_move(mixmod,&(results->accepted_absorb),&(results->proposed_absorb),pr_ej_G,pr_ej_Gm1);
      }
    }
    
    if(selectVariables)
    {
      j = BLCA_random_integer( d );
      BLCA_update_model_by_variable_include_exclude_old(mixmod,&(results->accepted_include_exclude_variable),&(results->proposed_include_exclude_variable),j);
      if(mixmod->hprior_model)
      {
        //update the prior probability variable inclusion using hyperprior
        j = 0;
        for( i=0; i<d; i++ ) j += vind[i];
        mixmod->prior_prob_variable_include = rbeta( j + mixmod->hprior_model_a0, mixmod->d - j + mixmod->hprior_model_b0 );
      }
    }
    
    if( verbose && l == num_burnin-1 ) Rprintf("\nBurn-in completed...");
    
    if( l == num_burnin-1 ) BLCA_reset_results( results );
    
    if(l > num_burnin-1 && (l+1-num_burnin)%thin_by == 0)
    {
      
      itmod = ((l+1-num_burnin)/thin_by)-1;
      
      if(selectVariables)
      {
        for( j=0; j<d; j++ ) variable_inclusion_indicator[ itmod + j*gap ] = vind[j];
        
        /* for( i=0; i<d; i++)
         {
         if( vind[i] ) results->variable_prob_inclusion[i] += 1.;
         }*/
      }
      /*else
       {
       for( j=0; j<d; j++ ) variable_inclusion_indicator[ itmod + j*gap ] = 1;
       }*/
      
      //results->ngroups[itmod] = mixmod->G;
      n_groups[itmod] = mixmod->G;
      prior_include[itmod] = mixmod->prior_prob_variable_include;	
      for( j=0; j<n; j++ ) group_memberships[ itmod + j*gap ] = z[j];
      
      
      /*store the value of the log posterior here- full value (incl prior for no. components)*/
      log_joint_posterior[itmod] = BLCA_get_full_log_posterior(mixmod);
      
      if( verbose && ( (l+1-num_burnin)/thin_by )%verbose_update == 0 ) Rprintf("\n%d of %d samples completed....",(l+1-num_burnin)/thin_by,gap);
      
    }
    
  }
  
  if( verbose ) Rprintf("\n");
  
  
  /*if(selectVariables){
   for(i=0;i<mixmod->d;i++){
   //results->variable_prob_inclusion[i] /= ((num_iteration-num_burnin)/thin_by);
   }
  }*/
  
  return(results);
  
}


void BLCA_analysis_MCMC_Gibbs_sampler( struct mix_mod *mixmod, int num_iteration, int num_burnin, int thin_by, int *group_memberships, 
                                       double *group_weights, double *variable_probabilities, double *log_joint_posterior, double *log_like, int sample_missing_data, int n_missing, 
                                       int *imputed_missing_values, int *position_missing, int verbose, int verbose_update )
{
  //this is the non-collapsed form of the model which uses Gibbs updates for all unknowns...
  // in this model there is no search for the number of groups, and no variable selection moves
  
  int g, i, j, k, l, c, p, idx, g_new, y_old, y_new, ind, itmod, 
  n = mixmod->n, d = mixmod->d, G = mixmod->G, *y, *z = mixmod->z, *N, *pmiss, *vind = mixmod->varindicator, *ncat = mixmod->ncat, *grp_mem_samp = group_memberships;
  double a, s, m, *prz, *prob, *w = mixmod->weights, eps = 1E-8, *beta_prior, *wei_samp = group_weights;
  struct component *comp;
  
  prz = (double *)calloc(G,sizeof(double));

  
  int gap = (int)( num_iteration/thin_by ) , gap_;
  
  if( verbose ) Rprintf("Initializing sampler... starting burn-in.\n");
  
  for( l=0; l<num_iteration+num_burnin; l++ )
  {
    
    R_CheckUserInterrupt();
    
    y = mixmod->y;
    for( k=0; k<n; k++ )
    {
      for( g=0; g<G; g++ )
      {
        comp = mixmod->components[g];
        a = w[g];
        if( a < eps ) a = eps;
        prz[g] = log( a ) ;
        for( j=0; j<d; j++ )
        {
          if( vind[j] )
          {
            prob = comp->prob_variables[j];
            a = prob[ y[j] ];
            if( a < eps ) a = eps;
            prz[g] += log( a );	
          }	
        }
      }
      
      m = BLCA_get_max( prz, G );
      s = 0.0;
      for( g=0; g<G; g++ )
      {
        prz[g] -= m;
        prz[g] = exp( prz[g] );
        s += prz[g];
      }
      for( g=0; g<G; g++ ) prz[g] /= s;
      
      g_new = BLCA_sample_discrete( prz, G );
      
      if(g_new != z[k])
      {
        //take out of component counts in old and put into new
        BLCA_add_to_component( mixmod->components[ z[k] ], y, mixmod, -1 );
        BLCA_add_to_component( mixmod->components[ g_new ], y, mixmod, 1 );
        z[k] = g_new;
      }
      
      y += d;
    }		
    
    //sample the weights
    s = 0.0;
    for( g=0; g<G; g++ )
    {
      w[g] = rgamma( (double) mixmod->components[g]->n_g + mixmod->alpha_prior[g] , 1.0 ) ;
      s += w[g];
    }
    for( g=0; g<G; g++ ) w[g] /= s;
    
    
    //sample the variable probabilities
    for( g=0; g<G; g++ )
    {
      comp = mixmod->components[g];
      for( j=0; j<d; j++ )
      {
        if( vind[j] )
        {
          //only sample if the variable is included in the model
          prob = comp->prob_variables[j];
          N = comp->N[j];
          beta_prior = mixmod->beta_prior[g][j];
          s = 0.0;
          for( c=0; c<ncat[j]; c++ )
          {
            prob[c] =  rgamma( (double)N[c] + beta_prior[c] , 1.0 ) ;
            s += prob[c];
          }
          for( c=0; c<ncat[j]; c++ ) prob[c] /= s;
        }
      }
    }
    
    // impute the missing data before carrying out the update step
    if( sample_missing_data )
    {
      pmiss = position_missing;
      for( k=0; k<n_missing; k++ )
      {
        y = mixmod->y;
        
        i = pmiss[0]; 
        j = pmiss[1];
        
        if( vind[j] )
        {
          g = z[i];
          comp = mixmod->components[g];
          prob = comp->prob_variables[j];
          y += d*i ;
          // take old vector out of component
          BLCA_add_to_component( comp, y, mixmod, -1 );
          
          y_new = BLCA_sample_discrete( prob, ncat[j] );
          y[j] = y_new;
          //mixmod->Y[j][i] = y_new;
          
          // add new vector to the component
          BLCA_add_to_component( comp, y, mixmod, 1 );
        }
        pmiss += 2;
      }
    }
    
    //storage of results
    if(l > num_burnin-1 && (l+1-num_burnin)%thin_by == 0)
    {
      
      if( verbose && l == num_burnin ) 
        Rprintf("Burn-in completed...");
      
      itmod = ((l+1-num_burnin)/thin_by)-1;
      
      //for( k=0; k<n; k++ ) group_memberships[itmod * n + k ] = z[k];
      for( k=0; k<n; k++ ) grp_mem_samp[k] = z[k];
      grp_mem_samp += n;
      
      p = 0;
      for( j=0; j<d; j++)
      {
        gap_ = p * G * gap;
        for( g=0; g<G; g++ )
        {
          prob = mixmod->components[g]->prob_variables[j];
          for( c=0; c<ncat[j]; c++ )
            variable_probabilities[ gap_ + itmod*ncat[j]*G + g*ncat[j] + c ] = prob[c];
        }
        p += ncat[j];
      }
      
      //for( g=0; g<G; g++ ) group_weights[itmod*G + g] = w[g];
      for( g=0; g<G; g++ ) wei_samp[g] = w[g];
      wei_samp += G;
      
      // if any missing data 
      
      if( sample_missing_data )
      {
        for( k=0; k<n_missing; k++ )
        {
          y = mixmod->y;
          i = position_missing[ 2 * k ];
          j = position_missing[ 2 * k + 1 ];
          y += d * i;
          imputed_missing_values[ itmod*n_missing + k ] = y[j];
        }
      }			
      
      log_joint_posterior[ itmod ] = BLCA_get_full_log_posterior_x2(mixmod);
      log_like[ itmod ] = BLCA_get_log_likelihood(mixmod);
      
      if( verbose && (l+1-num_burnin)%verbose_update == 0 )
        Rprintf("\n%d of %d samples completed...",l+1-num_burnin,num_iteration);
      
    }	
  }
  
  //if(verbose) Rprintf("\n");
  
  free(prz);
  //for(i=0;i<mixmod->d;i++) free(v[i]);
  //free(v);
  //free(order);
  
  return;
  
}


void BLCA_analysis_EM( struct mix_mod *mixmod, int max_iterations, int *iterations, double *membership_probabilities, 
                       double *group_weights, double *variable_probabilities, double *log_likelihood, int MAP, double tol, double *eps, int *converged ) 
{
  
  int i, j, g, k, c, p, iter = 0, n = mixmod->n, d = mixmod->d, G = mixmod->G, *ncat = mixmod->ncat ;
  double cond = DBL_MAX, llike_new, llike_old = -DBL_MAX, llike_inf, c_inf, *prob ;
  
  //need to make a modification to this for when MAP == TRUE (not done yet)
  
  while( cond > tol  && iter < max_iterations )
  {
    
    R_CheckUserInterrupt();
    
    BLCA_E_step( mixmod ); // the log prior has to be computed in here!
    BLCA_M_step( mixmod );
    
    if( mixmod->EM_MAP ) mixmod->log_like += mixmod->log_prior;
    log_likelihood[ iter ] = mixmod->log_like;
    
    if( iter > 1 )
    {
      // use Aitken's acceleration to get convergence
      c_inf = ( log_likelihood[iter] - log_likelihood[iter-1] ) / ( log_likelihood[iter-1] - log_likelihood[iter-2] );
      llike_inf = log_likelihood[iter-2] + ( log_likelihood[iter-1] - log_likelihood[iter-2] ) / (1.0-c_inf) ;
      cond = fabs( llike_inf - llike_old ) ;
      llike_old = llike_inf;
    }
    
    iter++;
    
  }
  
  if( cond < tol ) *converged = 1; else *converged = 0;
  
  //number of iterations to converge (if converged)
  *iterations = iter;
  *eps = cond;
  
  //store the results before returning
  for( g=0; g<G ; g++ ) group_weights[g] = mixmod->weights[g] ;
  
  for( k=0; k<n; k++ ) 
  {
    for( g=0; g<G; g++ )
      membership_probabilities[ (mixmod->n)*g + k ] = mixmod->s[k][g] ;
  }
  
  p = 0;
  int gap_;
  
  for( j=0; j<d; j++ )
  {
    gap_ = p * G;
    for( g=0; g<G; g++ )
    {
      prob = mixmod->components[g]->prob_variables[j];
      for( c=0; c<ncat[j]; c++ ) 
        variable_probabilities[ gap_ + g * ncat[j] + c ] = prob[c] ;
    } 
    
    p += ncat[j] ;
  }
  
  return;
}



void BLCA_E_step( struct mix_mod *mixmod )
{
  int k, g, j, c, G = mixmod->G, n = mixmod->n, d = mixmod->d, 
    BOOT = mixmod->BOOT, EM_MAP=mixmod->EM_MAP, *boot_idx, *y, *vind = mixmod->varindicator, *ncat = mixmod->ncat;
  double a, z, r=0.0, llik=0.0, lprior=0.0, max, eps = 1E-10, *prob, *s, *beta_prior, *alpha_prior, *w = mixmod->weights;
  struct component *comp;
  
  double *sk = (double *)calloc( n*G, sizeof(double));
  s = sk;
  
  llik = 0.0;
  
  // the appropriate data is pointed to if BOOT==TRUE  
  y = mixmod->y;
  
  for( k=0; k<n; k++ )
  {
    for( g=0; g<G; g++ )
    {
      
      comp = mixmod->components[g];
      
      a = w[g] ;
      if( a < eps ) a = eps ;
      s[g] = log( a );
      
      for( j=0; j<d; j++ )
      {
        if( vind[j] )
        {
          prob = comp->prob_variables[j] ;
          c = y[j] ;
          a = prob[c] ;
          if( a < eps ) a = eps ; 
          s[g] += log( a ) ;
        }
      }	
    }
    s += G;
    y += d;
  }
  
  s = sk;
  for( k=0; k<n; k++ )
  {
    llik += BLCA_get_log_sum_exp( s, G );
    for( g=0; g<G; g++ ) mixmod->s[k][g] = s[g];
    s += G;
  }
  mixmod->log_like = llik; 
  
  // need to compute log prior here for prior to correspond to same par set as llik
  if( EM_MAP )
  {
    alpha_prior = mixmod->alpha_prior;
    for( g=0; g<G; g++ )
    {
      comp = mixmod->components[g];
      for( j=0; j<d; j++ )
      {
        if( vind[j] )
        {
          z = 0.0;
          beta_prior = mixmod->beta_prior[g][j];
          prob = comp->prob_variables[j];
          for( c=0; c<ncat[j]; c++ ) 
          {
            a = prob[c];
            if( a < eps ) a = eps;
            lprior += ( beta_prior[c] - 1.0 ) * log( a );
            lprior -= lgamma( beta_prior[c] );
            z += beta_prior[c];
          }
          lprior += lgamma( z );
        }
      }
      lprior -= lgamma( alpha_prior[g] );
      a = w[g];
      if( a < eps ) a = eps;
      lprior += ( alpha_prior[g] - 1.0 ) * log( a );
      r += alpha_prior[g];
    }
    lprior += lgamma( r );
    mixmod->log_prior = lprior;
  }
  
  free(sk);
  return;
}


void BLCA_M_step( struct mix_mod *mixmod )
{
  int k, g, j, c, *y, G = mixmod->G, n = mixmod->n, d = mixmod->d, EM_MAP = mixmod->EM_MAP, BOOT = mixmod->BOOT,  
    *boot_idx, *ncat = mixmod->ncat, *vind = mixmod->varindicator;
  double a, b, s = 0.0, r=0.0, eps=1E-10, lprior=0.0, *w, *sk, *prob, *beta_prior, *alpha_prior = mixmod->alpha_prior;
  
  struct component *comp;
  
  // if bootstrapping need indices
  y = mixmod->y;
  
  // get the sum of the probabilities
  w = mixmod->weights;
  for( g=0; g<G; g++ ) w[g] = 0.0;
  for( k=0; k<n; k++ )
  {
    sk = mixmod->s[k];
    for(g=0; g<G; g++ ) w[g] += sk[g];
  }
  
  // reset probabilities to zero- will accumulate values
  for( g=0; g<G; g++ )
  {
    comp = mixmod->components[g];
    
    // set all vars to zero
    for( j=0; j<d; j++ )
    { 
      if( vind[j] ) 
      {
        prob = comp->prob_variables[j];
        for( c=0; c<ncat[j]; c++ ) prob[c] = 0.0;
      }
    }
  }
  
  for( k=0; k<n; k++ )
  {
    sk = mixmod->s[k];
    // contributions from datum k to all item probs
    for( g=0; g<G; g++ )
    {
      comp = mixmod->components[g];
      for( j=0; j<d; j++ )
      {
        if( vind[j] )
        {
          prob = comp->prob_variables[j];
          c = y[j];
          prob[c] += sk[g];
        }
      }
    }
    // move pointer along
    y += d;
  }
  
  for( g=0; g<G; g++ )
  {
    
    comp = mixmod->components[g];
    
    // gather up values of log prior
    for( j=0; j<d; j++ )
    {
      if( vind[j] )
      {
        prob = comp->prob_variables[j];
        beta_prior = mixmod->beta_prior[g][j];
        
        s = 0.;
        
        for( c=0; c<ncat[j]; c++ )
        {
          if( EM_MAP )
          {
            prob[c] += beta_prior[c] - 1.0 ;
            s += beta_prior[c] - 1.0 ; 
          }
          else
          {
            b = w[g] ;
            if( b < eps ) b = eps; 
            prob[c] /= b ;
          }
        }
        
        if( EM_MAP )
        {
          b = w[g];
          if( b < eps ) b = eps; 
          for( c=0; c<ncat[j]; c++ ) 
          {
            prob[c] /= ( b + s );
          }
        }
        
      }
    }
    
    if( EM_MAP )
    {
      w[g] +=  alpha_prior[g] - 1.0 ;
      r += alpha_prior[g] - 1.0;
    }
    else
    {
      w[g] /= n ;
    }
    
  }
  
  if( EM_MAP )
  {
    for( g=0; g<G; g++ ) w[g] /= ( n + r );
  }
  
  return;
}

/*-------------------------------------- VB algorithm -----------------------------------*/

void BLCA_analysis_VB( struct mix_mod *mixmod, int max_iterations, int *iterations, double *group_probabilities, double *group_weights, 
                       double *sd_group_weights, double *vb_pars_group_weights, double *prob_variables, double *sd_prob_variables, double *vb_pars_prob_variables, double *lower_bound, double tol, int *converged, double *log_post )
{
  int i, j, g, k, c, p, iter = 0, n = mixmod->n, d = mixmod->d, G = mixmod->G, *vind = mixmod->varindicator, *ncat = mixmod->ncat ;
  
  double cond = DBL_MAX, bound_new, bound_old = -DBL_MAX, bound_inf, c_inf ;
  
  while( cond > tol  && iter < max_iterations )
  {
    
    R_CheckUserInterrupt();
    
    BLCA_VB_alpha_beta_step( mixmod );
    
    BLCA_VB_phi_step( mixmod );
    
    lower_bound[ iter ] = BLCA_get_VB_bound( mixmod ) ; 
    
    if( iter > 1 )
    {
      // use Aitken's acceleration to get convergence
      c_inf = ( lower_bound[iter] - lower_bound[iter-1] ) / ( lower_bound[iter-1] - lower_bound[iter-2] );
      bound_inf = lower_bound[iter-2] + ( lower_bound[iter-1] - lower_bound[iter-2] ) / (1.0-c_inf) ;
      cond = fabs( bound_inf - bound_old ) ;
      bound_old = bound_inf;
    }
    
    iter++;
  }
  
  if( cond < tol ) *converged = 1; else *converged = 0;
  
  //number of iterations to converge (if converged)
  *iterations = iter;
  
  //store the results before returning
  
  double sum_alpha_ud = 0.0, var, *alpha_ud = mixmod->alpha_ud ;
  for( g=0; g<G ; g++ ) sum_alpha_ud += alpha_ud[g] ;
  
  for( g=0; g<G ; g++ )  
  {
    group_weights[g] = alpha_ud[g] / sum_alpha_ud ;
    vb_pars_group_weights[g] = alpha_ud[g] ;
    mixmod->weights[g] = group_weights[g] ;
    var = alpha_ud[g] * ( sum_alpha_ud - alpha_ud[g] ) / ( sum_alpha_ud * sum_alpha_ud * ( sum_alpha_ud + 1.0 ) );
    sd_group_weights[g] = sqrt( var );
  }
  
  double *s;
  for( k=0; k<n; k++ ) 
  {
    s = mixmod->s[k] ;
    for( g=0; g<G; g++ )
      group_probabilities[ k * G + g ] = s[g] ;
  }
  
  p = 0;
  int gap_;
  
  //convert the beta_ud into expect values of the theta's
  
  struct component *comp;
  double sum_b = 0.0, *beta_ud, ***var_prob_variables ;
  
  var_prob_variables = (double ***)malloc( G * sizeof(double**) );
  for( g=0; g<G; g++ )
  {
    var_prob_variables[g] = (double **)malloc( d * sizeof(double *) );
    for( j=0; j<d; j++ ) 
      var_prob_variables[g][j] = (double *)calloc( ncat[j], sizeof(double) );
  }
  
  for( g=0; g<G; g++ )
  {
    comp = mixmod->components[g] ;
    
    for( j=0; j<d; j++ )
    {
      if( vind[j] )
      {
        beta_ud = comp->beta_ud[j];
        sum_b = 0.0;
        for( c=0; c<ncat[j]; c++ ) sum_b += beta_ud[c] ;
        for( c=0; c<ncat[j]; c++ ) 
        {
          comp->prob_variables[j][c] = beta_ud[c] / sum_b;
          var_prob_variables[g][j][c] = beta_ud[c] * ( sum_b - beta_ud[c] ) / ( sum_b * sum_b * ( sum_b + 1.0 )  ) ;
        }
      }
    }
    
  }
  
  for( j=0; j<d; j++ )
  {
    gap_ = p * G;
    for( g=0; g<G; g++ )
    {
      for( c=0; c<ncat[j]; c++ )
      {
        prob_variables[ gap_ + g * ncat[j] + c ] = mixmod->components[g]->prob_variables[j][c] ;
        sd_prob_variables[ gap_ + g * ncat[j] + c ] = sqrt( var_prob_variables[g][j][c] ) ;
        vb_pars_prob_variables[ gap_ + g * ncat[j] + c ] = mixmod->components[g]->beta_ud[j][c] ;
      }
    } 
    p += ncat[j] ;
  }
  
  for( g=0; g<G; g++ )
  {
    for( j=0; j<d; j++ ) free( var_prob_variables[g][j] );
    free( var_prob_variables[g] );
  }
  free( var_prob_variables );

  //compute the log of the posterior... 
  int EM_map = mixmod->EM_MAP;
  mixmod->EM_MAP = TRUE;
  *log_post = BLCA_get_log_likelihood( mixmod );
  mixmod->EM_MAP = EM_map;
  
  return;
}

void BLCA_VB_phi_step( struct mix_mod *mixmod )
{
  int k, g, j, c, n = mixmod->n, G = mixmod->G, d = mixmod->d, *vind = mixmod->varindicator, *y = mixmod->y;
  double lse, *s, **di_beta_ud, *di_sum_beta_ud;
  struct component *comp;
  
  for( k=0; k<n; k++ )
  {
    s = mixmod->s[k];
    for( g=0; g<G; g++ )
    {
      comp = mixmod->components[g] ;
      di_beta_ud = comp->di_beta_ud ;
      di_sum_beta_ud = comp->di_sum_beta_ud  ;
      
      s[g] = mixmod->di_alpha_ud[g] - mixmod->di_sum_alpha_ud ;
      
      for( j=0; j<d; j++ )
      {
        if( vind[j] )
        {
          c = y[j] ;
          s[g] += di_beta_ud[j][c] - di_sum_beta_ud[j] ; //check that beta_ud correspond
        }
      }	
    }
    
    lse = BLCA_get_log_sum_exp( s, G );
    //s += G;
    y += d;
  }
  
 /* s = sk;
  for( k=0; k<n; k++ )
  {
    for( g=0; g<G; g++ ) mixmod->s[k][g] = s[g];
    s += G;
  }*/
  
 // free(sk);
  
  return;	
}

void BLCA_VB_alpha_beta_step( struct mix_mod *mixmod )
{
  
  //update the alphas first
  int k, g, j, c, n = mixmod->n, d = mixmod->d, G = mixmod->G, *ncat = mixmod->ncat, *vind = mixmod->varindicator, *y = mixmod->y ;
  double *colsums = (double *)calloc( G, sizeof(double) ), *s, *alpha_ud = mixmod->alpha_ud, w, *beta_ud, *beta_prior, *di_beta_ud, *di_sum_beta_ud;
  struct component *comp;
  
  for( k=0; k<n; k++ )
  {
    s = mixmod->s[k];
    for( g=0; g<G; g++ ) colsums[g] += s[g] ;
  }
  
  //alpha update 
  w = 0.;
  for( g=0; g<G; g++ )
  {
    alpha_ud[g] = colsums[g] + mixmod->alpha_prior[g]; 
    mixmod->di_alpha_ud[g] = digammaRN( alpha_ud[g] );
    w += alpha_ud[g];
    //initialize the beta update -- is this so necessary?
    comp = mixmod->components[g];
    for( j=0; j<d; j++ )
    {
      if( vind[j] )
      {
        beta_ud = comp->beta_ud[j] ;
        beta_prior = mixmod->beta_prior[g][j];
        for( c=0; c<ncat[j]; c++ )
          beta_ud[c] = beta_prior[c] ;
      }
    }
  }
  mixmod->di_sum_alpha_ud = digammaRN( w );
  
  //beta update -- need to think about the nesting here... 
  for( k=0; k<n; k++ )
  {
    s = mixmod->s[k];
    for( g=0; g<G; g++ )
    {
      comp = mixmod->components[g];
      for( j=0; j<d; j++ )
      {
        if( vind[j] )
        {
          beta_ud = comp->beta_ud[j];
          c = y[j];
          beta_ud[c] += s[g];
        }
      }
    }
    //s += G;
    y += d;
  }

  
  for( g=0; g<G; g++ )
  {
    comp = mixmod->components[g];	
    //beta_ud = comp->beta_ud;
    //di_beta_ud = comp->di_beta_ud;
    di_sum_beta_ud = comp->di_sum_beta_ud;
    for( j=0; j<d; j++ )
    {
      if( vind[j] )
      {
        beta_ud = comp->beta_ud[j];
        di_beta_ud = comp->di_beta_ud[j];
        di_sum_beta_ud[j] = 0.0;
        for( c=0; c<ncat[j]; c++ )
        {
          di_beta_ud[c] = digammaRN( beta_ud[c] ) ;
          di_sum_beta_ud[j] += beta_ud[c] ;
        }
        di_sum_beta_ud[j] = digammaRN( di_sum_beta_ud[j] ) ;
      }
    }
  }
  
  free( colsums );
  return;
}

void BLCA_analysis_Boot( struct mix_mod *mixmod, int boot_samples, int max_iterations, int *boot_samp_idx, double *group_probs,
                         double *group_weights, double *prob_variables, double *log_posterior_boot, double *log_likelihood_boot, double tol, int verbose, int verbose_update )
{
  
  int b, iter, i, j, k, g, l, c, p, ch, idx, g_new, ind, gap = boot_samples, gap_,
    n = mixmod->n, G = mixmod->G, d = mixmod->d, *ncat = mixmod->ncat, *boot_y = (int *)calloc(n*d, sizeof(int)), *ptr, *ptr_ref, *ptr_boot_samp_idx = boot_samp_idx;
  double cond = DBL_MAX, llike_new, c_inf, llike_inf, llike_old = -DBL_MAX, *prob, *wei_samp = group_weights, *grp_probs = group_probs ;
  double *log_likelihood = (double *)calloc(max_iterations, sizeof(double));
  
  //store the probabilities & weights from original EM run
  struct mix_mod *mixmod_0 = BLCA_clone_mixmod( mixmod );
  mixmod_0->y = mixmod->y;
  
  for( b=0; b<boot_samples; b++ )
  {
    // generate indexes for a bootstrap sample
    ptr_ref = mixmod_0->y;
    ptr = boot_y;
    
    ptr_ref += ptr_boot_samp_idx[0] * d ;
    for( j=0; j<d; j++ ) ptr[j] = ptr_ref[j];
    ptr += d;
    
    for( k=1; k<n ; k++ ) 
    {
      ch =  ptr_boot_samp_idx[k] - ptr_boot_samp_idx[k-1];
      if( ch > 0 ) ptr_ref += ch * d ;
      for( j=0; j<d; j++ ) ptr[j] = ptr_ref[j];
      ptr += d;
    }
    // increment for next bootstrap sample
    ptr_boot_samp_idx += n;
    
    mixmod->y = boot_y;
    
    // reset the component weights and item probabilities from original fit
    for( g=0; g<G; g++ ) 
    {
      BLCA_copy_component( mixmod_0->components[g], mixmod->components[g], mixmod );
      mixmod->weights[g] = mixmod_0->weights[g];
    }
    
    cond = DBL_MAX; llike_old = -DBL_MAX; iter = 0;
    
    while( cond > tol  && iter < max_iterations )
    {
      
      R_CheckUserInterrupt();
      
      BLCA_E_step( mixmod );
      
      BLCA_M_step( mixmod );
      
      if( mixmod->EM_MAP ) mixmod->log_like += mixmod->log_prior;
      
      log_likelihood[ iter ] = mixmod->log_like;
      
      if( iter > 1 )
      {
        // use Aitken's acceleration to get convergence
        c_inf = ( log_likelihood[iter] - log_likelihood[iter-1] ) / ( log_likelihood[iter-1] - log_likelihood[iter-2] );
        llike_inf = log_likelihood[iter-2] + ( log_likelihood[iter-1] - log_likelihood[iter-2] ) / (1.0-c_inf) ;
        cond = fabs( llike_inf - llike_old ) ;
        llike_old = llike_inf;
      }
      
      iter++;
    }
    
    log_posterior_boot[b] = mixmod->log_like;
    log_likelihood_boot[b] = BLCA_get_log_likelihood( mixmod ); //mixmod->log_like;
    
    p = 0;
    //storage of results	
    for( j=0; j<d; j++ )
    {
      gap_ = p * G * gap ;
      for( g=0; g<G; g++ )
      {
        prob = mixmod->components[g]->prob_variables[j];
        for( c=0; c<ncat[j]; c++ ) prob_variables[ gap_ + b*ncat[j]*G + g*ncat[j] + c ] = prob[c];
      }
      p += ncat[j];
    }
    
    for( g=0; g<G; g++ ) wei_samp[g] = mixmod->weights[g];
    wei_samp += G;
    
    // storage of group probabilities
    for( k=0; k<n; k++ ) 
    {
      for( g=0; g<G; g++ ) grp_probs[ n*g + k ] = mixmod->s[k][g] ;
    }
    grp_probs += n * G;
    
    if( verbose && (b+1)%verbose_update == 0 )
      Rprintf("\n%d of %d bootstrap samples completed....",b+1,boot_samples);
    
  }
  
  free(boot_y);
  free(log_likelihood);
  BLCA_free_mixmod( mixmod_0 );
  
  return;
  
}


