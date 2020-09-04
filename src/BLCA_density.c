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
			
	Last modification of this code: Wed 29 July 2020    */

#include "BLCA_density.h"

double BLCA_log_normalizing_constant_model( int G, struct mix_mod *mixmod )
/*returns the log of the normalizing constant for a model with G components*/
{
	double z, s=0.0;
	int k;
	
	if( mixmod->prior_type == 0 ) s = lgamma(G*mixmod->alpha) - G*lgamma(mixmod->alpha) - lgamma(mixmod->n+G*mixmod->alpha);
	
	if( mixmod->prior_type == 1 )
	{
		z = 0.0;
		for( k=0; k<G; k++ ) 
		{
			z += mixmod->alpha_prior[k];
			s -= lgamma(mixmod->alpha_prior[k]);
		}
		s += lgamma( z ) - lgamma( mixmod->n + z ) ;
	}
	
	return(s); 	
}

double BLCA_l_prior_variable_include( int D, struct mix_mod *mixmod )
{
	
	double l;
	
	l = D * log(mixmod->prior_prob_variable_include) + (mixmod->d - D) * log(1.-mixmod->prior_prob_variable_include) ;
	
	return(l);

}

double BLCA_get_full_log_posterior(struct mix_mod *mixmod)
{

  int j, g, v = 0, G = mixmod->G, d = mixmod->d, *wis = mixmod->whereis;
	double log_full_posterior = 0.0;
	

	/*model normalizing constant*/
	
	log_full_posterior += BLCA_log_normalizing_constant_model( G, mixmod );
	
	/*components - discriminating*/
	
	for( g=0; g<G; g++ ) log_full_posterior += mixmod->components[ wis[g] ]->log_prob;
	
	/*undiscriminating*/
	
	log_full_posterior += mixmod->undiscriminating->log_prob;
	
	/*prior on variable inclusion*/
	for( j=0; j<d; j++ ) v += mixmod->varindicator[j];
	
	log_full_posterior += BLCA_l_prior_variable_include(d,mixmod);
	
	log_full_posterior += mixmod->log_prior_G[mixmod->G];
	
	return(log_full_posterior);
	
	
}

double BLCA_get_full_log_posterior_x2(struct mix_mod *mixmod)
{

  int j, g, c, G = mixmod->G, d = mixmod->d, *vind = mixmod->varindicator, *ncat = mixmod->ncat, *N;
	double log_full_posterior=0.0, a, s=0.0, r=0.0, eps=1E-10, *alpha_prior = mixmod->alpha_prior, *beta_prior, *prob, *w = mixmod->weights;
	struct component *comp;

	for( g=0; g<G; g++ )
	{
	  comp = mixmod->components[g];
		s += alpha_prior[g];
		log_full_posterior -= lgamma( alpha_prior[g] );
		a = w[g];
		if( a < eps ) a = eps;
		log_full_posterior += (comp->n_g + alpha_prior[g] - 1.0) * log( a );
		for( j=0; j<d; j++ )
		{
			if( vind[j] )
			{
				r = 0.0;
			  N =  comp->N[j];
			  beta_prior = mixmod->beta_prior[g][j];
			  prob = comp->prob_variables[j];
				for( c=0; c<ncat[j]; c++ )
				{
					log_full_posterior -= lgamma( beta_prior[c] ); 
				  a = prob[c];
				  if( a < eps ) a = eps;
					log_full_posterior += ( N[c] + beta_prior[c] - 1.0) * log( a );				
					r += beta_prior[c] ;
				}
				log_full_posterior += lgamma(r);
				//log_full_posterior += lgamma(mixmod->ncat[j]*mixmod->beta) - mixmod->ncat[j]*lgamma(mixmod->beta);
			}
		}
	}
	log_full_posterior += lgamma(s);
	
	return(log_full_posterior);
}

double BLCA_get_log_likelihood(struct mix_mod *mixmod)
{
 	/*get log likelihood using stable log-sum-exp evaluation-- only to be used for EM*/
 	int g, k, j, c, G = mixmod->G, n = mixmod->n, d = mixmod->d, 
 	  *ncat = mixmod->ncat, *vind = mixmod->varindicator, *y = mixmod->y;
	double a, llik = 0.0, s=0.0, r=0.0, eps=1E-10, *prob, *w = mixmod->weights, *beta_prior, *alpha_prior, *prob_variables ;
  double *ld = (double *)calloc(G,sizeof(double)) ;
	struct component *comp;
	
	for( k=0; k<n; k++ )
	{
		for( g=0; g<G; g++ )
		{
		  comp = mixmod->components[g];
			ld[g] = 0.0;
			a = w[g] ;
			if( a < eps ) a = eps;
			ld[g] += log( a );
			for( j=0; j<d; j++ )
			{
				if( vind[j] )
				{
				  prob = comp->prob_variables[j];
					a = prob[ y[j] ];
					if( a < eps ) a = eps;
					ld[g] += log( a ) ; 
				}
			}
		}
		llik += BLCA_get_log_sum_exp( ld, G ) ; 
	  y += d;
	}
	
	// additional terms for prior if looking for MAP 
	if( mixmod->EM_MAP )
	{
	  alpha_prior = mixmod->alpha_prior ; 
		for( g=0; g<G; g++ ) 
		{
			s += alpha_prior[g];
			llik -= lgamma( alpha_prior[g] );
			llik += ( alpha_prior[g] - 1.0 ) *  log( mixmod->weights[g] ) ;
			for( j=0; j<d; j++ )
			{
				if( vind[j] )
				{
				  beta_prior = mixmod->beta_prior[g][j] ;
				  prob_variables = mixmod->components[g]->prob_variables[j];
					r = 0.0;
					for( c=0; c<ncat[j]; c++ ) 
					{
						r += beta_prior[c];
						llik -= lgamma( beta_prior[c] );
						llik += ( beta_prior[c] - 1.0 ) * log( prob_variables[c] );
					}
					llik += lgamma(r);
				}
			}
		
		}
		
		llik += lgamma(s);
		
		//Rprintf("\n loglike 2 = %lf ", log_likelihood);
	}
	
	free(ld);
	return( llik ) ;
}


double BLCA_get_VB_bound( struct mix_mod *mixmod )
{
	
	int k, g, j, c, n = mixmod->n, d = mixmod->d, G = mixmod->G, *vind = mixmod->varindicator, *ncat = mixmod->ncat, *y = mixmod->y ;
	double ex_q_logp = 0.0, ex_q_logq = 0.0, sm, *s, *di_beta_ud, eps = 1E-10, *lg_sum_beta, *lg_beta_sum, *alpha_ud, *di_alpha_ud, *colsums = (double *)calloc(G, sizeof(double)) ;
	struct component *comp;
	
	//s = mixmod->s[0] ;
	
	ex_q_logp = mixmod->lg_sum_alpha - mixmod->lg_alpha_sum ; //prior terms
	
	for( k=0; k<n; k++ )
	{
	  s = mixmod->s[k] ; 
		for( g=0; g<G; g++ )
		{
		  colsums[g] += s[g] ; 
		  
			comp = mixmod->components[g] ;
			//ex_q_logp += s[g] * ( mixmod->di_alpha_ud[g] - mixmod->di_sum_alpha_ud ) ;
		
			for( j=0; j<d; j++ )
			{
				if( vind[j] )
				{
					c = y[j] ;
					ex_q_logp += s[g] * ( comp->di_beta_ud[j][c] - comp->di_sum_beta_ud[j] ) ;
				}	
			}
			
			if( s[g] > eps ) ex_q_logq += s[g] * log( s[g] ); //don't do anything if below threshold
		}
		y += d ;
	}
	
	for( g=0; g<G; g++ ) ex_q_logp += colsums[g] * ( mixmod->di_alpha_ud[g] - mixmod->di_sum_alpha_ud ) ;
	
	sm = 0.0;
	for( g=0; g<G; g++ )
	{
		ex_q_logp += (mixmod->alpha_prior[g] - 1.0) * ( mixmod->di_alpha_ud[g] - mixmod->di_sum_alpha_ud ) ;
		
		comp = mixmod->components[g];
		
		lg_sum_beta = mixmod->lg_sum_beta[g];
		lg_beta_sum = mixmod->lg_beta_sum[g];
	  
		for( j=0; j<d; j++ )
		{
			if( vind[j] )
			{
				sm = 0.0;
			  ex_q_logp += lg_sum_beta[j] - lg_beta_sum[j];
				for( c=0; c<ncat[j]; c++ )
				{
				  ex_q_logp += (mixmod->beta_prior[g][j][c] - 1.0) * ( comp->di_beta_ud[j][c] - comp->di_sum_beta_ud[j] ) ; 
					ex_q_logq += ( comp->beta_ud[j][c] - 1.0 ) * ( comp->di_beta_ud[j][c] - comp->di_sum_beta_ud[j] ) - lgamma( comp->beta_ud[j][c] ) ;
					sm  += comp->beta_ud[j][c]; 
				}
				ex_q_logq += lgamma( sm ) ;
			}
		}
	}
	
	sm = 0.0;
	alpha_ud = mixmod->alpha_ud ;
	di_alpha_ud = mixmod->di_alpha_ud;
	for( g=0; g<G; g++ )
	{
		ex_q_logq += ( alpha_ud[g] - 1.0 ) * ( di_alpha_ud[g] - mixmod->di_sum_alpha_ud  ) - lgamma( alpha_ud[g] )  ;
	  sm += alpha_ud[g] ;
	}
	ex_q_logq += lgamma( sm );
	
	free(colsums) ; 
	
	return( ex_q_logp - ex_q_logq );
}
