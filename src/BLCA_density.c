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
			
	Last modification of this code: Thu 04 Aug 2016 21:43:55 IST    */

#include "BLCA_density.h"

double BLCA_log_normalizing_constant_model(int G,struct mix_mod *mixmod)
/*returns the log of the normalizing constant for a model with G components*/
{
	double z;
	
	z =  lgamma(G*mixmod->alpha) - G*lgamma(mixmod->alpha) - lgamma(mixmod->n+G*mixmod->alpha);
	
	return(z); 	
}

double BLCA_l_prior_variable_include(int D,struct mix_mod *mixmod)
{
	
	double l;
	
	l = D*log(mixmod->prior_prob_variable_include) + (mixmod->d - D)*log(1.-mixmod->prior_prob_variable_include);
	
	return(l);

}

double BLCA_get_full_log_posterior(struct mix_mod *mixmod)
{

	double log_full_posterior = 0.;
	int i,d;

	/*model normalizing constant*/
	
	log_full_posterior += BLCA_log_normalizing_constant_model(mixmod->G,mixmod);
	
	/*components - discriminating*/
	
	for(i=0;i<mixmod->G;i++){
		log_full_posterior += mixmod->components[ mixmod->whereis[i] ]->log_prob;
	}
	
	/*undiscriminating*/
	
	log_full_posterior += mixmod->undiscriminating->log_prob;
	
	/*prior on variable inclusion*/
	d = 0;
	for(i=0;i<mixmod->d;i++) d += mixmod->varindicator[i];
	
	log_full_posterior += BLCA_l_prior_variable_include(d,mixmod);
	
	log_full_posterior += mixmod->log_prior_G[mixmod->G];
	
	return(log_full_posterior);
	
	
}

double BLCA_get_full_log_posterior_x2(struct mix_mod *mixmod)
{

	double log_full_posterior;
	int j,i,c;

	/*model normalizing constant*/
	
	log_full_posterior = lgamma(mixmod->G*mixmod->alpha) - mixmod->G*lgamma(mixmod->alpha);
	
	for(i=0;i<mixmod->G;i++){
		log_full_posterior += (mixmod->components[i]->n_g + mixmod->alpha -1.)*log(mixmod->weights[i]);
		for(j=0;j<mixmod->d;j++){
		
			if( mixmod->varindicator[j] )
			{
				for(c=0;c<mixmod->ncat[j];c++){
					log_full_posterior += (mixmod->components[i]->N[j][c] + mixmod->beta - 1.)*log(mixmod->components[i]->prob_variables[j][c]);
				}
				log_full_posterior += lgamma(mixmod->ncat[j]*mixmod->beta) - mixmod->ncat[j]*lgamma(mixmod->beta);
			}
		}
	}
	
	return(log_full_posterior);
}

double BLCA_get_log_likelihood(struct mix_mod *mixmod)
{
 	/*get log likelihood using stable log-sum-exp evaluation*/

	double log_likelihood = 0., *ld = calloc(mixmod->G , sizeof(double)) ;
	int g, k, j, c;
	
	for( k=0; k<mixmod->n; k++ ){
		for( g=0; g<mixmod->G; g++ ){
			ld[g] = 0.;
			ld[g] += log( mixmod->weights[g] );
			for( j=0; j<mixmod->d; j++ ){
				if( mixmod->varindicator[j] )
				{
					c = mixmod->Y[j][k] ;
					ld[g] += log( mixmod->components[g]->prob_variables[j][c] ) ; 
				}
			}
		
		}
		log_likelihood += BLCA_get_log_sum_exp( ld, mixmod->G ) ; 
	}
	
	// additional terms for prior if looking for MAP 
	if( mixmod->EM_MAP )
	{
		log_likelihood += lgamma( mixmod->G * mixmod->alpha ) - mixmod->G  * lgamma( mixmod->alpha ) ;
		for( g=0; g<mixmod->G; g++ ) 
		{
			log_likelihood += ( mixmod->alpha - 1. ) *  log( mixmod->weights[g] ) ;
			for( j=0; j<mixmod->d; j++ )
			{
				if( mixmod->varindicator[j] )
				{
					log_likelihood += lgamma( mixmod->ncat[j] * mixmod->beta ) - mixmod->ncat[j] * lgamma( mixmod->beta );
					for( c=0; c<mixmod->ncat[j]; c++ ) log_likelihood += ( mixmod->beta - 1. ) * log( mixmod->components[g]->prob_variables[j][c] );
				}
			}
		
		}
	
	}
	
	free(ld);
	return( log_likelihood ) ;
}

double BLCA_get_VB_bound( struct mix_mod *mixmod )
{
	
	int k, g, j, c;
	double Elogp = 0., Elogq = 0.;
	struct component *comp;
	
	for( k=0; k<mixmod->n; k++ )
	{
		for( g=0; g<mixmod->G; g++ )
		{
			comp = mixmod->components[g] ;
		
			Elogp += mixmod->s[k][g] * ( mixmod->di_alpha_ud[g] - mixmod->di_sum_alpha_ud ) ;
		
			for( j=0; j<mixmod->d; j++ )
			{
				if( mixmod->varindicator[j] )
				{
					c = mixmod->Y[j][k] ;
					Elogp += comp->di_beta_ud[j][c] - comp->di_sum_beta_ud[j] ;
				}	
			}
			
			if( mixmod->s[k][g] < 1E-10 )
				Elogq += 0.;
			else
				Elogq += mixmod->s[k][g] * log( mixmod->s[k][g] );
			
		}
	}
	
	double sum_j = 0;
	
	for( g=0; g<mixmod->G; g++ )
	{
		comp = mixmod->components[g];
		
		for( j=0; j<mixmod->d; j++ )
		{
			if( mixmod->varindicator[j] )
			{
				sum_j = 0.;
				Elogp += lgamma( mixmod->ncat[j] * mixmod->beta ) - mixmod->ncat[j] * lgamma( mixmod->beta ) ;
				for( c=0; c<mixmod->ncat[j]; c++ )
				{
					Elogp += ( mixmod->beta - 1. ) * ( comp->di_beta_ud[j][c] - comp->di_sum_beta_ud[j] ) ;
					Elogq += -lgamma( comp->beta_ud[j][c] ) + ( comp->beta_ud[j][c] - 1. ) * ( comp->di_beta_ud[j][c] - comp->di_sum_beta_ud[j] ) ;
					sum_j  += comp->beta_ud[j][c]; 
				}
				Elogq += lgamma( sum_j ) ;
			}
		}
	}
	
	Elogp += lgamma( mixmod->G * mixmod->alpha ) - mixmod->G * lgamma( mixmod->alpha ) ; 
	
	double sum_g = 0.;
	for( g=0; g<mixmod->G; g++ )
	{
		Elogp += ( mixmod->alpha - 1. ) * ( mixmod->di_alpha_ud[g] - mixmod->di_sum_alpha_ud ) ;
		Elogq += - mixmod->alpha_ud[g] + ( mixmod->alpha_ud[g] - 1. ) * ( mixmod->di_alpha_ud[g] - mixmod->di_sum_alpha_ud  ) ;
		sum_g += mixmod->alpha_ud[g] ;
	}
	Elogq += lgamma( sum_g ) ;
	
	return( Elogp - Elogq );
}
